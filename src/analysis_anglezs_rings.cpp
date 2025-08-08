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

//***************** Partially contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include "analysis_anglezs_rings.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_ANGLEZS_RINGS::ANALYSIS_ANGLEZS_RINGS(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename) : ANALYSIS()
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
}

void ANALYSIS_ANGLEZS_RINGS::init() {
}

vector<vector<float>>ANALYSIS_ANGLEZS_RINGS::compute_2D_vector() {
  //  clock_t t;  //debug

//    vector<vector<float>> angles;

   // t = clock(); //debug
    sel1->anglezs.clear();
    for (auto &segment:sel1->segments_ind) {
	vector<float> com = {0.0,0.0,0.0};
	vector<float> angles_seg;
	vector<float> dist ;//= {0.0,0.0,0.0};
// calculate the COM of each ring
	for (int ind : segment) {
	    com[0] = com[0] + system->x[ind];
	    com[1] = com[1] + system->y[ind];
	    com[2] = com[2] + system->z[ind];
	}
        for (int dim = 0; dim < 3; dim++) {
	    com[dim] = com[dim]/float(segment.size());
	}
// Calculate the zangle of all N's
	for (int ind : segment) {
	    vector<float> coords;
	    coords.push_back(system->x[ind]);
	    coords.push_back(system->y[ind]);
	    coords.push_back(system->z[ind]);
//	        cout << "frame: " << i << " atom:" << ind << endl;
	    dist = getDistPoints(coords, com);
	    float anglez = angle_z(dist);
   	    angles_seg.push_back(anglez);
 	}
	sel1->anglezs.push_back(angles_seg);
	
    }
//    t = clock() - t;
//    cout << "angle calculation cost time: " << ((double) t) / CLOCKS_PER_SEC << " seconds" << endl;
    return sel1->anglezs;
}

ANALYSIS_ANGLEZS_RINGS::~ANALYSIS_ANGLEZS_RINGS()
{
    system = NULL;
    sel1 = NULL;
}

