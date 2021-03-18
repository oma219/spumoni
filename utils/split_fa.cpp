//
//  split_fa.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <zlib.h>
#include <stdio.h>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread)

void split_file(std::string& path, std::size_t n_seqs, std::size_t n_blocks)
{
    std::size_t seqs_per_block = 0, residual = 0;
    if (n_seqs % n_blocks == 0) { seqs_per_block = n_seqs / n_blocks; residual = n_seqs / n_blocks; }
    else if (n_blocks == 2) {seqs_per_block = n_seqs / 2; residual = (n_seqs % 2) + (n_seqs / 2);}
    else { seqs_per_block = n_seqs / (n_blocks - 1); residual = n_seqs % (n_blocks - 1); }
    
    std::size_t last_index = path.find_last_of('.'); std::string remove_gz = path.substr(0, last_index);
    last_index = remove_gz.find_last_of('.', last_index);
    std::string base_name = path.substr(0, last_index);
    
    int l;
    gzFile fp;
    kseq_t *seq;
    fp = gzopen(path.c_str(), "r");
    seq = kseq_init(fp);
    for (std::size_t i = 0; i < n_blocks - 1; i++)
    {
        std::cout << "\rSplitting sequences... " << std::to_string(i + 1) << "/" << n_blocks << "  "
                  << std::to_string((double(i + 1) / double(n_blocks)) * 100) << "%" << std::flush;
        
        std::string out_path = base_name + "_" + std::to_string(i + 1) + ".fa";
        std::ofstream out_file(out_path);
        std::size_t it = 0;
        while (it < seqs_per_block)
        {
            l = kseq_read(seq);
            if (seq->seq.l > 0)
            {
                out_file.put('>'); out_file.write(seq->name.s, seq->name.l); out_file.put('\n');
                out_file.write(seq->seq.s, seq->seq.l); out_file.put('\n');
            }
            it++;
        }
        out_file.close();
    }
    // remember to write residual to last file
    std::string out_path = base_name + "_" + std::to_string(n_blocks) + ".fasta";
    std::ofstream out_file(out_path);
    std::size_t it = 0;
    std::cout << "\rSplitting sequences... " << std::to_string(n_blocks) << "/" << n_blocks << "  "
              << std::to_string((double(n_blocks) / double(n_blocks)) * 100) << "%" << std::flush;
    while (it < residual)
    {
        l = kseq_read(seq);
        if (seq->seq.l > 0)
        {
            out_file.put('>'); out_file.write(seq->name.s, seq->name.l); out_file.put('\n');
            out_file.write(seq->seq.s, seq->seq.l); out_file.put('\n');
        }
        it++;
    }
    out_file.close();
    
    
    // free
    kseq_destroy(seq);
    gzclose(fp);
}

// Count sequences
std::size_t count_seqs(std::string& path)
{
    std::size_t sequences_count = 0;
    int l;
    gzFile fp;
    kseq_t *seq;
    fp = gzopen(path.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        sequences_count++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    
    return sequences_count;
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <in.seq> <n blocks>\n", argv[0]);
        return 1;
    }
    
    std::string path = argv[1];
    std::cout << "In path: " << path << std::endl;
    std::size_t n_blocks = std::stoi(argv[2]);
    std::cout << "Blocks: " << n_blocks << std::endl;
    
    std::cout << "Reading sequences...";
    std::size_t n_seq = count_seqs(path);
    std::cout << " done. N: " <<  n_seq << std::endl;
    
    std::cout << "Splitting sequences...";
    split_file(path, n_seq, n_blocks);
    std::cout << " done." << std::endl;
    
    return 0;
}
