#include <algorithm>
#include <iterator>
#include <iostream>
#include <map>
#include "graph.h"

std::mutex io_lock;

void segment_printer_body::operator()(printer_input input, printer_node::output_ports_type & op){

    auto &payload = get<0>(input); 
    size_t token  = get<1>(input);

    auto &index = get<0>(payload);
    auto &fw_segments = get<1>(payload);
    auto &rc_segments = get<2>(payload);
    auto &block_start = get<3>(payload);
    auto &block_len   = get<4>(payload);
    auto &block_index = get<5>(payload);
    auto &inter_ref_start = get<6>(payload);
    auto &inter_ref_end   = get<7>(payload);
    auto &inter_query_start = get<8>(payload);
    auto &inter_query_end   = get<9>(payload);
    size_t block_end = block_start + block_len;

    std::string base_filename;
    std::string segment_filename;

    uint32_t num_hsps = fw_segments.size() + rc_segments.size();

    if(num_hsps > 0){

        uint32_t start_r_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_ref_start) - chr_start.cbegin(); 
        uint32_t end_r_chr   = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_ref_end) - chr_start.cbegin(); 
        uint32_t start_q_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_query_start) - chr_start.cbegin(); 
        uint32_t end_q_chr   = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_query_end) - chr_start.cbegin(); 

        uint32_t curr_q_chr_index = start_q_chr;
        std::string curr_q_chr = chr_name[start_q_chr];

        size_t curr_q_chr_start = chr_start[start_q_chr]; 
        size_t curr_q_chr_end = curr_q_chr_start + chr_len[start_q_chr];

        std::string out_str;
        FILE* segmentFile;

        if(fw_segments.size() > 0){

            base_filename = "tmp"+std::to_string(index)+".block"+std::to_string(block_index)+".r"+std::to_string(block_start)+".plus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for (auto e: fw_segments) {
                size_t seg_r_start = block_start + e.ref_start;
                size_t seg_q_start = block_start + e.query_start;
                size_t r_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_r_start) - chr_start.cbegin() - 1;

                if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                    out_str = chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n";
                }
                else{
                    size_t q_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_q_start) - chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = chr_name[curr_q_chr_index];
                    curr_q_chr_start = chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + chr_len[curr_q_chr_index];

                    out_str = chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t+\t" + std::to_string(e.score) + "\n";
                }
                fprintf(segmentFile, "%s", out_str.c_str());
            }

            fclose(segmentFile);
        }

        start_q_chr = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_query_start) - chr_start.cbegin();
        end_q_chr   = std::upper_bound(chr_start.cbegin(), chr_start.cend(), inter_query_end) - chr_start.cbegin();

        curr_q_chr_index = start_q_chr;
        curr_q_chr = chr_name[start_q_chr];
        curr_q_chr_start = chr_start[start_q_chr];
        curr_q_chr_end = curr_q_chr_start + chr_len[start_q_chr];

        if(rc_segments.size() > 0){
            base_filename = "tmp"+std::to_string(index)+".block"+std::to_string(block_index)+".r"+std::to_string(block_start)+".minus"; 
            segment_filename = base_filename+".segments";
            segmentFile = fopen(segment_filename.c_str(), "w");

            for(int r = rc_segments.size()-1; r >= 0; r--){
                auto e =  rc_segments[r];
                size_t seg_r_start = e.ref_start + block_start;
                size_t seg_q_start = (e.query_start + block_start);
                size_t r_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_r_start) - chr_start.cbegin() - 1;

                if(seg_q_start >= curr_q_chr_start && seg_q_start < curr_q_chr_end){
                    out_str = chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n";
                }
                else{

                    size_t q_index = std::upper_bound(chr_start.cbegin(), chr_start.cend(), seg_q_start) - chr_start.cbegin() - 1;
                    curr_q_chr_index = q_index;
                    curr_q_chr = chr_name[curr_q_chr_index];
                    curr_q_chr_start = chr_start[curr_q_chr_index];
                    curr_q_chr_end = curr_q_chr_start + chr_len[curr_q_chr_index];

                    out_str = chr_name[r_index] + '\t' + std::to_string(seg_r_start+1-chr_start[r_index]) + '\t' + std::to_string(seg_r_start+e.len+1-chr_start[r_index]) + '\t' + curr_q_chr + '\t' +  std::to_string(seg_q_start+1-curr_q_chr_start) + '\t' + std::to_string(seg_q_start+e.len+1-curr_q_chr_start) + "\t-\t" + std::to_string(e.score) + "\n";
                }
                fprintf(segmentFile, "%s", out_str.c_str());
            }
            fclose(segmentFile);
        }
    }

    get<0>(op).try_put(token);
};
