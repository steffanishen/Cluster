/*
 *  read_dcd : c++ class + main file example for reading a CHARMM dcd file
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//***************** Partially contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

#include "cluster.hpp"

#define PI 3.14159265

using namespace std;

CLUSTER::CLUSTER(int cluster_ind)
{
    this->cluster_ind = cluster_ind;
    this->parent_cluster_ind = -1;
    kid_clusters.clear();
    residue_members.clear();
    head_cap.resize(2,-1);
    tail_cap.resize(2,-1);
}


CLUSTER::~CLUSTER()
{
    residue_members.clear();
    kid_clusters.clear();
}
