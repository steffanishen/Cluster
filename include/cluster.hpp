/*
    read_dcd : c++ class + main file example for reading a CHARMM dcd file
    Copyright (C) 2013  Florent Hedin
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//***************** Fully contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#ifndef CLUSTER_HPP_INCLUDED
#define CLUSTER_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include "atom.hpp"

using namespace std;

class CLUSTER
{
public:
    
    /*content of ICNTRL : non detailed ones are 0 */
    
    int parent_cluster_ind; // Number of atoms
    int cluster_ind; // Number of atoms
   
    vector<vector<int>> residue_members; 

    vector<CLUSTER*> kid_clusters;

    int head_C1 = -1;
    int tail_C2 = -1;

    vector<int> head_cap;
    vector<int> tail_cap;
   

    CLUSTER(int cluster_ind);
    //coordinates stored in simple precision (IMPORTANT)
    
    //protected methods

    //write_cluster_file();   

    ~CLUSTER();
};

#endif // DCD_HPP_INCLUDED
