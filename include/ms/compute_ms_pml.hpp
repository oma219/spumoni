/* 
 *  compute_ms_pml.hpp
 *  Copyright (C) 2020 Omar Ahmed
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 */
 /*
  * File: compute_ms_pml.hpp
  * Description: Header file for compute_ms_pml.cpp
  *
  * Authors: Omar Ahmed
  * Start Date: October 13, 2021
  */

#ifndef COMPUTE_MS_PML_H
#define COMPUTE_MS_PML_H

/* Function Declarations */
int run_spumoni_ms_main(SpumoniRunOptions* run_opts);
int run_spumoni_main(SpumoniRunOptions* run_opts);

#endif /* End of include of COMPUTE_MS_PML_H */