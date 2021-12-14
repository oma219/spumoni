/* build_spumoni - Build the SPUMONI index
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
   \file build_spumoni.cpp
   \brief build_spumoni.cpp Build the SPUMONI index.
   \author Massimiliano Rossi
   \date 03/03/2021
*/

#include <iostream>
#define VERBOSE
#include <common.hpp>
#include <sdsl/io.hpp>
#include <spumoni.hpp>
//#include <malloc_count.h>

int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index
  verbose("Building the SPUMONI index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
  ms_pointers<> ms(args.filename, true);
  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  //verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  std::string outfile = args.filename + ms.get_file_extension();
  std::ofstream out(outfile);
  ms.serialize(out);


  t_insert_end = std::chrono::high_resolution_clock::now();

  //verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  //auto mem_peak = malloc_count_peak();
  //verbose("Memory peak: ", malloc_count_peak());
  
  /*
  size_t space = 0;
  if (args.memo)
  {
    sdsl::nullstream ns;

    size_t ms_size = ms.serialize(ns);
    verbose("MS size (bytes): ", ms_size);
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;
  */
  
  return 0;
}