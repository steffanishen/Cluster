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
#include "analysis_scale.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_SCALE ::ANALYSIS_SCALE (PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, float scaling) : ANALYSIS()
{
    this->system = system;
    this->sel1 = sel1;
    this->scaling = scaling;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
}

void ANALYSIS_SCALE ::init() {
}

void ANALYSIS_SCALE::compute_void() {

    FILE * fpdb;
    const char *filechar = filename.c_str();
    fpdb = fopen(filechar,"w");

    vector<float> shift = {0.0,0.0,0.0};
    vector<float> com_res = {0.0,0.0,0.0};
    vector<float> com_res_scaled = {0.0,0.0,0.0};

    for (int restempid = 0; restempid < sel1->segments_ind.size(); restempid++) {
	
	for (int dim = 0; dim < 3; dim++) {
	    shift[dim] = 0.0;
	    com_res[dim] = 0.0;
	    com_res_scaled[dim] = 0.0;
	}
        for (int ind : sel1->segments_ind[restempid]) {
	    com_res[0] = com_res[0] + system->x[ind];
	    com_res[1] = com_res[1] + system->y[ind];
	    com_res[2] = com_res[2] + system->z[ind];
	}
	int Natoms_res = sel1->segments_ind[restempid].size();
	for (int dim = 0; dim < 3; dim++) {
	    com_res[dim] = com_res[dim]/float(Natoms_res);
	    com_res_scaled[dim] = com_res[dim] * scaling;
	    shift[dim] = com_res_scaled[dim] - com_res[dim];
	}

        for (int ind : sel1->segments_ind[restempid]) {
 	    system->x[ind] = system->x[ind] + shift[0];
 	    system->y[ind] = system->y[ind] + shift[1];
 	    system->z[ind] = system->z[ind] + shift[2];
        }
    }

    for (int restempid = 0; restempid < sel1->residues.size(); restempid++) {
        for (int ind : sel1->residues[restempid]) {
            fprintf(fpdb,"%4s%7d  %-4s%-4s%1d%4d%12.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM", ind + 1, system->atomname[ind].c_str(), system->resname[ind].c_str(), 1, system->resid[ind], system->x[ind], system->y[ind], system->z[ind], 0,0);
        }
    }

}

ANALYSIS_SCALE::~ANALYSIS_SCALE()
{
    system = NULL;
    sel1 = NULL;
}

