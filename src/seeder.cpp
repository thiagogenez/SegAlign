#include "graph.h"
#include "store.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);

printer_input seeder_body::operator()(seeder_input input) {

    seeder_payload &payload = get<0>(input);

    auto &chrom = get<0>(payload);
    auto &interval_data = get<1>(payload);

    size_t token = get<1>(input);

    std::vector<segment> fw_segments;
    std::vector<segment> rc_segments;
    fw_segments.clear();
    rc_segments.clear();

    uint64_t index = 0;
    uint64_t transition_index = 0;
    uint64_t seed_offset;

    uint32_t inter_ref_start = interval_data.ref_start;
    uint32_t inter_ref_end   = interval_data.ref_end;
    uint32_t inter_query_start = interval_data.query_start;
    uint32_t inter_query_end   = interval_data.query_end;
    uint32_t num_invoked   = interval_data.num_invoked;
    uint32_t num_intervals = interval_data.num_intervals;

    size_t block_start   = chrom.block_start;
    uint32_t block_len   = chrom.block_len;
    uint32_t block_index = chrom.block_index;

    uint32_t rc_query_start = seq_len - inter_query_end;
    uint32_t rc_query_end   = seq_len - inter_query_start;

    fprintf (stderr, "Chromosome block %u interval %u/%u (%u:%u)\n", block_index, num_invoked, num_intervals, inter_query_start, inter_query_end);

    if(cfg.strand == "plus" || cfg.strand == "both"){
        for (uint32_t i = inter_query_start; i < inter_query_end; i += cfg.wga_chunk_size) {

            //end position
            uint32_t e = std::min(i + cfg.wga_chunk_size, inter_query_end);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = i; j < e; j++) {

                index = GetKmerIndexAtPos(seq_DRAM->buffer, block_start+j, cfg.seed_size);
                if (index != ((uint32_t) 1 << 31)) {
                    seed_offset = (index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 

                    if (cfg.transition) {
                        for (int t=0; t < sa->GetKmerSize(); t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segment> anchors = g_SeedAndFilter(seed_offset_vector, false, 0);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    fw_segments.insert(fw_segments.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    if(cfg.strand == "minus" || cfg.strand == "both"){
        for (uint32_t i = rc_query_start; i < rc_query_end; i += cfg.wga_chunk_size) {
            uint32_t e = std::min(i + cfg.wga_chunk_size, rc_query_end);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();
            for (uint32_t j = i; j < e; j++) {
                index = GetKmerIndexAtPos(seq_rc_DRAM->buffer, block_start+j, cfg.seed_size);
                if (index != ((uint32_t) 1 << 31)) {
                    seed_offset = (index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 
                    if (cfg.transition) {
                        for (int t=0; t < sa->GetKmerSize(); t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segment> anchors = g_SeedAndFilter(seed_offset_vector, true, 0);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    rc_segments.insert(rc_segments.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(num_invoked, fw_segments, rc_segments, block_start, block_len, block_index, inter_ref_start, inter_ref_end, inter_query_start, inter_query_end), token);
}
