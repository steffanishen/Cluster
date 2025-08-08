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

//***************** Contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unordered_map>
#include "analysis_packing.hpp"


#define PI 3.14159265

using namespace std;

//ANALYSIS_RANDOM_WALK::ANALYSIS_RANDOM_WALK(PSF *system, GROUP *sel1, int vector1d, int vector2d,int voidf, float bl, string name1, string name2, string name3, string name4, string name_ref1, string name_ref2, string filename, vector<PSF*> monomers)
ANALYSIS_PACKING::ANALYSIS_PACKING(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, float lx, float ly, float lz) //constructor
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->lx = lx;
    this->ly = ly;
    this->lz = lz;
}

void ANALYSIS_PACKING::init() {
}


void ANALYSIS_PACKING::compute_void() {
    FILE * fpdb;
    const char *filechar = filename.c_str();
    fpdb = fopen(filechar,"w");
  //  vector<float> next_bond(3);

    cout << "Nresidues: " << sel1->residues.size() << endl;


///////////determine the initial orientation of the C-C bond////////////////////////////////////////

/////////////////////////////end determining the C-C bond///////////////////////

    srand(time(NULL));

    for (int segtempid = 0; segtempid < sel1->segments_ind.size(); segtempid++) {
// calculate the geometrical center of the segment
	float center_x = 0.0;
	float center_y = 0.0;
	float center_z = 0.0;
        for (int ind : sel1->segments_ind[segtempid]) {
	    center_x += system->x[ind];
	    center_y += system->y[ind];
	    center_z += system->z[ind];
	}
	float segment_norm = float (sel1->segments_ind[segtempid].size());
	center_x = center_x / segment_norm;
	center_y = center_y / segment_norm;
	center_z = center_z / segment_norm;

// place the chains at random coordinates

	float x_random = lx * (float) rand()/RAND_MAX - lx * 0.5;
	float y_random = ly * (float) rand()/RAND_MAX - ly * 0.5;
	float z_random = lz * (float) rand()/RAND_MAX - lz * 0.5;

	float x_shift = x_random - center_x;
	float y_shift = y_random - center_y;
	float z_shift = z_random - center_z;

        for (int ind : sel1->segments_ind[segtempid]) {
	    system->x[ind] += x_shift;
	    system->y[ind] += y_shift;
	    system->z[ind] += z_shift;
	}

// randomly rotate the chains by azimuth and polar angles
	float theta = 2 * PI * (float) rand()/RAND_MAX;

	vector<float> sp_azi = {0.0, 0.0, 1.0};
	vector<float> sp_pol = {1.0, 0.0, 0.0};

        for (int ind : sel1->segments_ind[segtempid]) {
            vector<float> dist = {system->x[ind] - x_random,system->y[ind] - y_random,system->z[ind] - z_random};
	    vector<float> dist_rot = rotate(dist,sin(theta),cos(theta),sp_azi);
	    system->x[ind] = x_random + dist_rot[0];
	    system->y[ind] = y_random + dist_rot[1];
	    system->z[ind] = z_random + dist_rot[2];
    	}

        for (int ind : sel1->segments_ind[segtempid]) {
            vector<float> dist = {system->x[ind] - x_random,system->y[ind] - y_random,system->z[ind] - z_random};
	    vector<float> dist_rot = rotate(dist,sin(theta),cos(theta),sp_pol);
	    system->x[ind] = x_random + dist_rot[0];
	    system->y[ind] = y_random + dist_rot[1];
	    system->z[ind] = z_random + dist_rot[2];
    	}

    }

    for (int restempid = 0; restempid < sel1->residues.size(); restempid++) {
        for (int ind : sel1->residues[restempid]) {
            fprintf(fpdb,"%4s%7d  %-4s%-4s%1d%4d%12.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM", ind + 1, system->atomname[ind].c_str(), system->resname[ind].c_str(), 1, system->resid[ind], system->x[ind], system->y[ind], system->z[ind], 0,0);
        }
    }

    fclose(fpdb);    
}


ANALYSIS_PACKING::~ANALYSIS_PACKING()
{
    system = NULL;
    sel1 = NULL;
}

