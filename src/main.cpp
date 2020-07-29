#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp> 
#include <iostream>
#include "math.h"
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include "kseq.h"
#include "zlib.h"
#include "graph.h"
#include "parameters.h"
#include "store.h"

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////////
KSEQ_INIT2(, gzFile, gzread)

struct timeval start_time, end_time, start_time_complete, end_time_complete;
long useconds, seconds, mseconds;

Configuration cfg;
SeedPosTable *sa;

DRAM *seq_DRAM = nullptr;
DRAM *seq_rc_DRAM = nullptr;
size_t seq_len = 0;

std::vector<std::string> chr_name;
std::vector<size_t>      chr_start;
std::vector<uint32_t>    chr_len;

std::vector<size_t>   seq_block_start;
std::vector<uint32_t> seq_block_len;
std::vector<uint32_t> seq_num_intervals;

////////////////////////////////////////////////////////////////////////////////

void RevComp(size_t rc_start, size_t start, size_t len) {

    size_t r = rc_start;
    for (size_t i = start+len; i> start; i--) {
        
        switch (seq_DRAM->buffer[i-1]) {
            case 'a': seq_rc_DRAM->buffer[r++] = 't';
                      break;

            case 'A': seq_rc_DRAM->buffer[r++] = 'T';
                      break;

            case 'c': seq_rc_DRAM->buffer[r++] = 'g';
                      break;

            case 'C': seq_rc_DRAM->buffer[r++] = 'G';
                      break;

            case 'g': seq_rc_DRAM->buffer[r++] = 'c';
                      break;

            case 'G': seq_rc_DRAM->buffer[r++] = 'C';
                      break;

            case 't': seq_rc_DRAM->buffer[r++] = 'a';
                      break;

            case 'T': seq_rc_DRAM->buffer[r++] = 'A';
                      break;

            case 'n': seq_rc_DRAM->buffer[r++] = 'n';
                      break;

            case 'N': seq_rc_DRAM->buffer[r++] = 'N';
                      break;

            case '&': seq_rc_DRAM->buffer[r++] = '&';
                      break;

            default: printf("Bad Nt char! '%c' %lu\n", seq_DRAM->buffer[i], i);
        }
    }
}

int main(int argc, char** argv){

    po::options_description desc{"Options"};
    desc.add_options()
        ("strand", po::value<std::string>(&cfg.strand)->default_value("both"), "strand to search - plus/minus/both")
        ("scoring", po::value<std::string>(&cfg.scoring_file), "Scoring file in LASTZ format")
        ("ambiguous", po::value<std::string>(&cfg.ambiguous), "ambiguous nucleotides - n/iupac")
        ("seed", po::value<std::string>(&cfg.seed_shape)->default_value("12of19"), "seed pattern-12of19(1110100110010101111)/14of22(1110101100110010101111)/an arbitrary pattern of 1s, 0s, and Ts ")
        ("notransition", po::bool_switch(&cfg.transition)->default_value(false), "don't allow one transition in a seed hit")
        ("step", po::value<uint32_t>(&cfg.step)->default_value(1), "Offset between the starting positions of successive target words considered for generating seed table")
        ("xdrop", po::value<int>(&cfg.xdrop)->default_value(910), "x-drop value for gap-free extension")
        ("hspthresh", po::value<int>(&cfg.hspthresh)->default_value(3000), "segment score threshold for high scoring pairs")
        ("queryhsplimit", po::value<uint32_t>(&cfg.queryhsplimit)->default_value(-1), "segment score threshold for high scoring pairs")
        ("noentropy", po::bool_switch(&cfg.noentropy)->default_value(false), "don't adjust low score segment pair scores using entropy factor after filtering stage")
        ("nogapped", po::bool_switch(&cfg.gapped)->default_value(false), "don't perform gapped extension stage")
        ("ydrop", po::value<int>(&cfg.ydrop)->default_value(9430), "y-drop value for gapped extension")
        ("gappedthresh", po::value<int>(&cfg.gappedthresh), "score threshold for gapped alignments")
        ("notrivial", po::bool_switch(&cfg.notrivial)->default_value(false), "Don't output a trivial self-alignment block if the target and query sequences are identical")
        ("format", po::value<std::string>(&cfg.output_format)->default_value("maf-"), "format of output file (same formats as provided by LASTZ) - lav, lav+text, axt, axt+, maf, maf+, maf-, sam, softsam, sam-, softsam-, cigar, BLASTN, differences, rdotplot, text")
        ("output", po::value<std::string>(&cfg.output), "output filename")
        ("markend", po::bool_switch(&cfg.markend), "write a marker line just before completion")
        ("wga_chunk", po::value<uint32_t>(&cfg.wga_chunk_size)->default_value(DEFAULT_WGA_CHUNK), "chunk sizes for GPU calls for Xdrop - change only if you are a developer")
        ("lastz_interval", po::value<uint32_t>(&cfg.lastz_interval_size)->default_value(DEFAULT_LASTZ_INTERVAL), "LASTZ interval for ydrop - change only if you are a developer")
        ("num_gpu", po::value<int>(&cfg.num_gpu)->default_value(-1), "Specify number of GPUs to use - -1 if all the GPUs should be used")
        ("num_neighbor_intervals", po::value<uint32_t>(&cfg.num_neighbor_intervals)->default_value(5), "Specify number of neighboring chunks for repeat masking per query chunk")
        ("debug", po::bool_switch(&cfg.debug)->default_value(false), "print debug messages")
        ("help", "Print help messages");

    po::options_description hidden;
    hidden.add_options()
        ("target", po::value<std::string>(&cfg.seq_filename)->required(), "target sequence file in FASTA format");

    po::options_description all_options;
    all_options.add(desc);
    all_options.add(hidden);

    po::positional_options_description p;
    p.add("target", 1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        if(!vm.count("help")){
            if(!vm.count("target")){
                fprintf(stderr, "You must specify a target file\n"); 
            }
        }

        fprintf(stderr, "Usage: run_segalign target query [options]\n"); 
        std::cout << desc << std::endl;
        return 1;
    }

    cfg.transition = !cfg.transition;
    cfg.gapped = !cfg.gapped;
    if(cfg.seed_shape == "12of19"){
        cfg.seed = "TTT0T00TT00T0T0TTTT"; 
        cfg.seed_size = 19;
    }
    else if(cfg.seed_shape == "14of22"){
        cfg.seed = "TTT0T0TT00TT00T0T0TTTT";
        cfg.seed_size = 22;
    }
    else{
        int seed_len = cfg.seed_shape.size();
        cfg.seed = cfg.seed_shape;
        for(int i = 0; i< seed_len; i++){
            if(cfg.seed_shape[i] == '1')
                cfg.seed[i] = 'T';
            else
                cfg.seed[i] = '0';
        }
        cfg.seed_size = cfg.seed.size();
    }

    if(vm.count("gappedthresh") == 0)
        cfg.gappedthresh = cfg.hspthresh; 

    int ambiguous_reward = -100;
    int ambiguous_penalty = -100;
    int fill_score = -100;
    int bad_score = -1000;
    std::string ambiguous_field = "x";

    std::vector <std::string> fields;
    boost::split( fields, cfg.ambiguous, boost::is_any_of( "," ) );
    ambiguous_field = fields[0];
    if(fields.size() == 3){
        ambiguous_reward  = std::stoi(fields[1]);
        ambiguous_penalty = -1*std::stoi(fields[2]);
    }
    else if (cfg.ambiguous == "n" || cfg.ambiguous == "iupac"){
        ambiguous_reward  = 0;
        ambiguous_penalty = 0;
    }

    if(vm.count("scoring") == 0){

        //ACGT
        int tmp_sub_mat[L_NT][L_NT] = {{   91, -114,  -31, -123},
                                 { -114,  100, -125,  -31},
                                 {  -31, -125,  100, -114},
                                 { -123,  -31, -114,  91}};

        for(int i = 0; i < L_NT; i++){
            for(int j = 0; j < L_NT; j++){
                cfg.sub_mat[i*NUC+j] = tmp_sub_mat[i][j];
            }
        }

        //lower case characters
        for(int i = 0; i < L_NT; i++){
            cfg.sub_mat[i*NUC+L_NT] = bad_score;
            cfg.sub_mat[L_NT*NUC+i] = bad_score;
        }
        cfg.sub_mat[L_NT*NUC+L_NT] = bad_score;

        //N
        if(ambiguous_field == "n" || ambiguous_field == "iupac"){
            for(int i = 0; i < N_NT; i++){
                cfg.sub_mat[i*NUC+N_NT] = ambiguous_penalty;
                cfg.sub_mat[N_NT*NUC+i] = ambiguous_penalty;
            }
            cfg.sub_mat[N_NT*NUC+N_NT] = ambiguous_reward;
        }
        else{
            for(int i = 0; i < N_NT; i++){
                cfg.sub_mat[i*NUC+N_NT] = bad_score;
                cfg.sub_mat[N_NT*NUC+i] = bad_score;
            }
            cfg.sub_mat[N_NT*NUC+N_NT] = bad_score;
        }

        //other IUPAC
        if(ambiguous_field == "iupac"){
            for(int i = 0; i < X_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = ambiguous_penalty;
                cfg.sub_mat[X_NT*NUC+i] = ambiguous_penalty;
            }
            cfg.sub_mat[X_NT*NUC+X_NT] = ambiguous_reward;
        }
        else{
            for(int i = 0; i < L_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = fill_score;
                cfg.sub_mat[X_NT*NUC+i] = fill_score;
            }

            for(int i = L_NT; i < X_NT; i++){
                cfg.sub_mat[i*NUC+X_NT] = bad_score;
                cfg.sub_mat[X_NT*NUC+i] = bad_score;
            }
            cfg.sub_mat[X_NT*NUC+X_NT] = fill_score;
        }

        for(int i = 0; i < E_NT; i++){
            cfg.sub_mat[i*NUC+E_NT] = -10*cfg.xdrop;
            cfg.sub_mat[E_NT*NUC+i] = -10*cfg.xdrop;
        }
        cfg.sub_mat[E_NT*NUC+E_NT] = -10*cfg.xdrop;
    }

    cfg.num_threads = tbb::task_scheduler_init::default_num_threads();
    cfg.num_threads = (cfg.num_threads == 1) ? 2 : cfg.num_threads;
    tbb::task_scheduler_init init(cfg.num_threads);

    if(cfg.debug){
        fprintf(stderr, "Target %s\n", cfg.seq_filename.c_str());
        fprintf(stderr, "Seed %s\n", cfg.seed.c_str());
        fprintf(stderr, "Seed size %d\n", cfg.seed_size);
        fprintf(stderr, "Transition %d\n", cfg.transition);
        fprintf(stderr, "Gapped %d\n",cfg.gapped);
        fprintf(stderr, "xdrop %d\n", cfg.xdrop);
        fprintf(stderr, "ydrop %d\n", cfg.ydrop);
        fprintf(stderr, "HSP threshold %d\n", cfg.hspthresh);
        fprintf(stderr, "gapped threshold %d\n", cfg.gappedthresh);
        fprintf(stderr, "ambiguous %s\n", cfg.ambiguous.c_str());

        for(int i = 0; i < NUC; i++){
            for(int j = 0; j < NUC; j++){
                fprintf(stderr, "%d ", cfg.sub_mat[i*NUC+j]);
            }
            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "Using %d threads\n", cfg.num_threads);

    cfg.num_gpu = g_InitializeProcessor (cfg.num_gpu, cfg.transition, cfg.wga_chunk_size, cfg.seed_size, cfg.sub_mat, cfg.xdrop, cfg.hspthresh, cfg.noentropy);

    seq_DRAM = new DRAM;
    seq_rc_DRAM = new DRAM;
    
    // Read query file
    if(cfg.debug){
	    gettimeofday(&start_time, NULL);
    }

    fprintf(stderr, "\nReading genome file ...\n");

    gzFile f_rd = gzopen(cfg.seq_filename.c_str(), "r");
    if (!f_rd) { 
        fprintf(stderr, "cant open file: %s\n", cfg.seq_filename.c_str()); 
        exit(7); 
    }
        
    kseq_t *kseq_rd = kseq_init(f_rd);
    uint32_t total_chr = 0;

    memset(seq_DRAM->buffer + seq_DRAM->bufferPosition, '&', 1);
    seq_DRAM->bufferPosition += 1;
    seq_len += 1;

    while (kseq_read(kseq_rd) >= 0) {
        uint32_t curr_chr_len = kseq_rd->seq.l;
        std::string curr_chr_name = std::string(kseq_rd->name.s, kseq_rd->name.l);

        chr_name.push_back(curr_chr_name);
        chr_start.push_back(seq_DRAM->bufferPosition);
        chr_len.push_back(curr_chr_len);

        if (seq_DRAM->bufferPosition + curr_chr_len > seq_DRAM->size) {
            fprintf(stderr, "Not enough memory in query DRAM. Need to reduce input query sequence size to less than 6GB!\n"); 
            exit(9); 
        }
        
        memcpy(seq_DRAM->buffer + seq_DRAM->bufferPosition, kseq_rd->seq.s, curr_chr_len);
        seq_DRAM->bufferPosition += curr_chr_len;
        seq_len += curr_chr_len;

        total_chr++;

        memset(seq_DRAM->buffer + seq_DRAM->bufferPosition, '&', 1);
        seq_DRAM->bufferPosition += 1;
        seq_len += 1;
    }

    RevComp(seq_rc_DRAM->bufferPosition, 0, seq_len);

    cfg.num_chr = total_chr;
    gzclose(f_rd);

    cfg.num_blocks = ceil((float) seq_len/SEQ_BLOCK_SIZE);

    size_t block_start = 0;
    size_t ref_block_start = 0;
    size_t ref_block_end = 0;
    uint32_t ref_block_len = 0;

    uint32_t num_intervals = 0;

    uint32_t interval_start;
    uint32_t interval_end;
    uint32_t curr_pos;
    uint32_t end_pos;

    float overlap_blocks = cfg.num_neighbor_intervals/2;
    uint32_t overlap_len = floor(overlap_blocks)*cfg.lastz_interval_size;
    uint32_t max_ref_interval_len = cfg.num_neighbor_intervals*cfg.lastz_interval_size;
    std::vector<seed_interval> interval_list;
    interval_list.clear();
    uint32_t total_intervals = 0;

    while(block_start < seq_len){

        if(block_start > (cfg.num_neighbor_intervals-1)*cfg.lastz_interval_size){
            ref_block_start = block_start-overlap_len;
            curr_pos = overlap_len; 
            end_pos = overlap_len;
        }
        else{
            ref_block_start = 0;
            curr_pos = 0;
            end_pos = 0;
        }

        ref_block_end = block_start+SEQ_BLOCK_SIZE+overlap_len;
        end_pos += SEQ_BLOCK_SIZE;

        if(ref_block_end > seq_len-1){
            ref_block_end = seq_len-1;
            end_pos = seq_len-1-ref_block_start;
        }

        ref_block_len = ref_block_end - ref_block_start+1; 
        seq_block_start.push_back(ref_block_start);
        seq_block_len.push_back(ref_block_len);
        block_start += SEQ_BLOCK_SIZE;

        num_intervals = 0;
        while (curr_pos < end_pos) {
            seed_interval inter;
            inter.query_start = curr_pos;
            inter.query_end = std::min(end_pos, curr_pos + cfg.lastz_interval_size);
            if(inter.query_start < overlap_len){
                inter.ref_start = 0;
            }
            else{
                inter.ref_start = inter.query_start - overlap_len;
            }

            inter.ref_end = inter.ref_start + max_ref_interval_len;

            if(inter.ref_end > end_pos){
                inter.ref_end = end_pos;
            }

            inter.num_invoked = 0;
            inter.num_intervals = 0;

            interval_list.push_back(inter);
            curr_pos += cfg.lastz_interval_size;
            num_intervals++;
        }
        seq_num_intervals.push_back(num_intervals);
    }

    total_intervals = interval_list.size();

//    for(int i = 24; i < 43; i++){
//        printf("%c", seq_DRAM->buffer[i]);
//    }
//    printf("\n");
//
//    for(int i = 5491038; i < 5491057; i++){
//        printf("%c", seq_rc_DRAM->buffer[i]);
//    }
//    printf("\n");

    if(cfg.debug){
    	gettimeofday(&end_time, NULL);
    	useconds = end_time.tv_usec - start_time.tv_usec;
    	seconds = end_time.tv_sec - start_time.tv_sec;
    	mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    	fprintf(stderr, "Time elapsed (loading complete query from file): %ld msec \n\n", mseconds);
    }

    // start alignment
    fprintf(stderr, "\nStart alignment ...\n");
    tbb::flow::graph align_graph;

    printer_node printer(align_graph, tbb::flow::unlimited, segment_printer_body());

    tbb::flow::function_node<seeder_input, printer_input> seeder(align_graph, tbb::flow::unlimited, seeder_body());

    tbb::flow::make_edge(seeder, printer);

    tbb::flow::join_node<seeder_input> gatekeeper(align_graph);

    tbb::flow::make_edge(gatekeeper, seeder);

    tbb::flow::buffer_node<size_t> ticketer(align_graph);

    // Allocate tickets
    for (size_t t = 0ull; t < cfg.num_threads; t++)
        ticketer.try_put(t);

    tbb::flow::make_edge(tbb::flow::output_port<0>(printer), ticketer);

    tbb::flow::make_edge(ticketer, tbb::flow::input_port<1>(gatekeeper));

    uint32_t blocks_sent = 0;
    bool send_block = true;
    uint32_t send_block_len;
    uint32_t send_block_start;
    uint32_t send_block_intervals;

    uint32_t invoked_block_intervals = 0;
    uint32_t block_intervals_num = 0;

    uint32_t completed_intervals = 0;  

    gettimeofday(&start_time, NULL);
    tbb::flow::source_node<seeder_payload> reader(align_graph,
        [&](seeder_payload &op) -> bool {

        while(true){
            if (send_block) {

                send_block_start = seq_block_start[blocks_sent];
                send_block_len   = seq_block_len[blocks_sent];
                send_block_intervals = seq_num_intervals[blocks_sent];

                fprintf(stderr, "\nSending block %u ...\n", blocks_sent);
                if(blocks_sent > 0)
                    g_clearQuery(0);
                g_SendQueryWriteRequest (send_block_start, send_block_len, 0);

                if(cfg.debug){
                    gettimeofday(&start_time_complete, NULL);
                }

                sa = new SeedPosTable (seq_DRAM->buffer, send_block_start, send_block_len, cfg.seed, cfg.step);

                if(cfg.debug){
                    gettimeofday(&end_time_complete, NULL);
                    useconds = end_time_complete.tv_usec - start_time_complete.tv_usec;
                    seconds = end_time_complete.tv_sec - start_time_complete.tv_sec;
                    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
                    fprintf(stderr, "\nTime elapsed (seed position table create and copy to GPU) : %ld msec \n", mseconds);
                }

                seeder_body::total_xdrop = 0;
                
                completed_intervals += invoked_block_intervals;

                invoked_block_intervals = 0;
                send_block = false;
            }
            else if(send_block_intervals == seeder_body::total_xdrop.load()){
                if(blocks_sent < cfg.num_blocks){
                    blocks_sent++;
                    send_block = true; 
                }
            }

            if(blocks_sent < cfg.num_blocks){
                if (invoked_block_intervals < send_block_intervals) {
                    reader_output& chrom = get<0>(op);
                    chrom.block_start = send_block_start;
                    chrom.block_len = send_block_len;
                    chrom.block_index = blocks_sent; 

                    seed_interval& inter = get<1>(op);
                    seed_interval curr_inter = interval_list[completed_intervals + invoked_block_intervals++];
                    inter.ref_start = curr_inter.ref_start;
                    inter.ref_end   = curr_inter.ref_end;
                    inter.query_start = curr_inter.query_start;
                    inter.query_end   = curr_inter.query_end;
                    inter.num_invoked = invoked_block_intervals;
                    inter.num_intervals = send_block_intervals;

                    return true;
                }
            }
            else if(blocks_sent == cfg.num_blocks && invoked_block_intervals == send_block_intervals){
                return false;
            }
            }

        }, true);

    tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));
    
    align_graph.wait_for_all();

    g_ShutdownProcessor();

    if(cfg.debug){
        gettimeofday(&end_time, NULL);
        seconds = end_time.tv_sec - start_time.tv_sec;
        fprintf(stderr, "Time elapsed (complete pipeline): %ld sec \n\n", seconds);
    	fprintf(stderr, "#seeds: %lu \n", seeder_body::num_seeds.load());
    	fprintf(stderr, "#seed hits: %lu \n", seeder_body::num_seed_hits.load());
    	fprintf(stderr, "#HSPs: %lu \n", seeder_body::num_hsps.load());
    }

    return 0;
}
