 /*
  * File: batch_loader.hpp 
  * Description: Header file for batch_loader.cpp
  *             
  * Start Date: April 21, 2022
  *
  * Note: The BatchLoader class is inspired by the implementation
  *       in the Kraken2 repository at this link 
  *       (https://github.com/DerrickWood/kraken2/blob/master/src/seqreader.h). It uses
  *       different logic to determine the batch size to load balance between
  *       multiple reader threads.
  */

#ifndef BATCH_LOADER_H
#define BATCH_LOADER_H

#include <spumoni_main.hpp>
#include <sstream> 

struct Read {
    std::string header_line; // name of read
    std::string id; // only header until first whitespace
    std::string seq; // sequence
    std::string qual; // qualities if needed

    query_input_type read_format;
};

class BatchLoader {
public:
    BatchLoader();
    ~BatchLoader() {};

    bool loadBatch(std::ifstream& input, size_t num_bases);
    bool grabNextRead(Read& curr_read);
    
private:
    query_input_type input_format;
    std::stringstream batch_stream; // will hold entire batch
    std::string batch_buffer; // will be used for reading
};

#endif /* End of BATCH_LOADER_H */