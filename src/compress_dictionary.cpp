/* compress_dictionary - Computes the compressed dictionary from prefix-free parse dictionary
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file compress_dictionary.cpp
   \brief compress_dictionary.cpp Computes the compressed dictionary from prefix-free parse dictionary.
   \author Massimiliano Rossi
   \date 16/09/2020
*/

#define VERBOSE
#include <iostream>
#include <common.hpp>
#include <spumoni_main.hpp>

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";
#if defined(_WIN32) || defined(_WIN64)
#define popen _popen
#define pclose _pclose
#endif
    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
    if (!pipe) {THROW_EXCEPTION(std::runtime_error("popen() failed!"));}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        THROW_EXCEPTION(std::runtime_error("Error occurred while reading popen() stream."));
    }
    pclose(pipe);
    return output;
}

int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  FORCE_LOG("compress_dict", "starting to compress dictionary.");
  auto overall_start = std::chrono::high_resolution_clock::now();

  // open output files (*.dicz and *.dicz.len)
  std::string dicz_filename = args.filename + ".dicz";
  std::string dicz_len_filename = args.filename + ".dicz.len";

  FILE *dicz;
  FILE *dicz_len;

  if ((dicz = fopen(dicz_filename.c_str(), "w")) == nullptr)
    FATAL_ERROR(("open() file " + std::string(dicz_filename) + " failed").data());

  if ((dicz_len = fopen(dicz_len_filename.c_str(), "w")) == nullptr)
    FATAL_ERROR(("open() file " + std::string(dicz_len_filename) + " failed").data());

  // open the dictionary and read it into an array of uint8_t
  std::string dict_filename = args.filename + ".dict";
  std::vector<uint8_t> dict;
  read_file(dict_filename.c_str(), dict);
    
  // counting the number of dollar at beginning of parse
  size_t i = 0, j = 0;
  while(dict[i++] == Dollar) {j++;}
  dict.erase(dict.begin(), dict.begin() + j);

  // computing the length of each phrase
  FORCE_LOG("compress_dict", "computing the length of each phrase");
  std::vector<size_t> lengths(1,0);

  for(auto chr: dict) {
    // skip the \x00
    if(chr == EndOfDict) continue;

    // hit end of phrase \x01
    if(chr == EndOfWord) lengths.push_back(0);
    else lengths.back()++;
  }
  if (lengths.back()==0)
    lengths.pop_back();

  FORCE_LOG("compress_dict", "found %d phrases in the dictionary", lengths.size());


  // start writing out to *.dicz and *.dicz.len files
  uint8_t* ptr = dict.data();
  bool empty_first_phrase = false;

  for (size_t i = 0; i < lengths.size(); i++) 
  {
    // we only want to write the portion of phrase minus trigger string
    size_t compressed_length = lengths[i] - args.w;

    // special case: starts with a trigger string
    if (i==0 && compressed_length == 0) {
      ptr += lengths[i] + 1; 
      empty_first_phrase = true;
      continue;
    } else if (i > 0 && compressed_length == 0) {
      error("encountered a length=0 phrase after removing trigger string, which should not occur.");
    }

    // write length as 4-byte integer to *.dicz.len
    if ((fwrite(&compressed_length, 4, 1, dicz_len)) != 1)
      FATAL_ERROR(("fwrite() file " + std::string(dicz_len_filename) + " failed").data());
    // write character as 1-byte char to *.dicz (ignoring last w chars)
    if ((fwrite(ptr, sizeof(uint8_t), compressed_length, dicz)) != compressed_length)
      FATAL_ERROR(("fwrite() file " + std::string(dicz_filename) + " failed").data());

    ptr += lengths[i] + 1;
  }
  fclose(dicz);
  fclose(dicz_len);

  // re-writes parse file to shift down all the phrase ids by 1 
  // since we removed the empty beginning phrase
  if (empty_first_phrase) {
    FORCE_LOG("compress_dict", "alert: found that the first phrase length is 0"
              " so we will rewrite *.parse file to generated correct SLP.");

    // read in all the phrase ids in parse
    std::string parse_filename = args.filename + ".parse";
    std::vector<uint32_t> parse_arr;
    read_file(parse_filename.c_str(), parse_arr);

    // make sure first phrase is lowest lexicographically and then remove it
    ASSERT((parse_arr[0] == 1),
           "parse should begin with lowest lexicographic phrase.");
    parse_arr.erase(parse_arr.begin());

    // rename the old parse file as *.parse_with_empty_phrase
    std::ostringstream command_stream;
    command_stream << "mv " << parse_filename << " " << (args.filename + ".parse_with_empty_phrase");
    auto mv_log = execute_cmd(command_stream.str().c_str());

    FORCE_LOG("compress_dict", "executed this command: %s", command_stream.str().data());

    // open new parse file for writing
    FILE* new_parse_file;
    if ((new_parse_file = fopen((args.filename + ".parse").c_str(), "w")) == nullptr)
      FATAL_ERROR(("open() file " + std::string(args.filename + ".parse") + " failed").data());

    // iterate through each element of parse and decrement by 1
    for (size_t i = 0; i < parse_arr.size(); i++) {
      ASSERT((parse_arr[i] > 1), "issue occurred when creating new parse file.");
      parse_arr[i]--;
      
      // write it out 
      if ((fwrite(&parse_arr[i], 4, 1, new_parse_file)) != 1)
        FATAL_ERROR(("fwrite() file " + std::string(args.filename + ".parse") + " failed").data()); 
    }
    fclose(new_parse_file);
  }

  auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - overall_start);
  FORCE_LOG("compress_dict", "elapsed time (s): %.5f", time.count());
  return 0;
}