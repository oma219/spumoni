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

#include <iostream>

#define VERBOSE

#include <common.hpp>

//#include <malloc_count.h>

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index

  verbose("Compressing the dictionary");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Open output files
  std::string dicz_filename = args.filename + ".dicz";
  std::string dicz_len_filename = args.filename + ".dicz.len";

  FILE *dicz;
  FILE *dicz_len;

  if ((dicz = fopen(dicz_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_filename) + " failed");

  if ((dicz_len = fopen(dicz_len_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_len_filename) + " failed");

  // Open the dictionary
  std::string dict_filename = args.filename + ".dict";
  std::vector<uint8_t> dict;
  read_file(dict_filename.c_str(), dict);

  // Start processing

  
  // Generating phrase lengths
  verbose("Generating phrase lengths");
  std::vector<size_t> lengths(1,0);
  
  // Counting the number of Dollars at the beginning
  size_t i = 0, j = 0;
  while(dict[i++] == Dollar)
    j++;
  dict.erase(dict.begin(), dict.begin() + j);

  for(auto chr: dict)
  {
    // Skip the Dollars
    if(chr == EndOfDict)
      continue;

    // Hit end of phrase
    if(chr == EndOfWord)
      lengths.push_back(0);
    else
      lengths.back()++;
  }

  if (lengths.back()==0)
    lengths.pop_back();

  verbose("Found", lengths.size(), " phrases ");

  verbose("Generating phrases");
  uint8_t* ptr = dict.data(); // Beginning of the current phrase
  for(auto length: lengths)
  {
    size_t compressed_length = length - args.w;

    if ((fwrite(&compressed_length, 4, 1, dicz_len)) != 1)
      error("fwrite() file " + std::string(dicz_len_filename) + " failed");

    if ((fwrite(ptr, sizeof(uint8_t), compressed_length, dicz)) != compressed_length)
      error("fwrite() file " + std::string(dicz_filename) + " failed");

    ptr += length + 1;
  }


  fclose(dicz);
  fclose(dicz_len);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  //verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  //auto mem_peak = malloc_count_peak();
  //verbose("Memory peak: ", malloc_count_peak());
  
  /*
  size_t space = 0;
  if (args.memo)
  {
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;
  */
  
  return 0;
}