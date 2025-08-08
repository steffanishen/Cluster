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
#include <stdio.h>
#include <limits>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

#include <math.h>
#include "analysis_density_profile_droplet.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_DENSITY_PROFILE_DROPLET::ANALYSIS_DENSITY_PROFILE_DROPLET(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, int nbins,string which_density_profile,float zshift,float dr)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->nbins = nbins;
    this->which_density_profile = which_density_profile;
    this->zshift = zshift;
    this->dr = dr;
}

void ANALYSIS_DENSITY_PROFILE_DROPLET::init() {
}

vector<float> ANALYSIS_DENSITY_PROFILE_DROPLET::compute_vector() {
    sel1->anglezs.clear();
    float x,y,z;
    int izbin;
    int iybin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float R_dense = 3.0;
    vector<float> ytemp;
    vector<float> ztemp;

    vector<float> profile(this->nbins,0.0);
    vector<double> box;
    box.resize(3);
    box[0] = system->box_first_frame[0];
    box[1] = system->box_first_frame[1]; 
    box[2] = system->box_first_frame[2];
    //box[0] = system->pbc[0];
    //box[1] = system->pbc[2];
    //box[2] = system->pbc[5];
    this->yshift = box[1] * 0.5;
    float dz = box[2]/float(this->nbins);
    float z_shift = 3.4;

    for (auto &segment:sel1->segments_ind) {
      for (int ind : segment) {
	      x = system->x[ind];
	      y = system->y[ind];
          z = system->z[ind];
        
          ytemp.push_back(y);
          ztemp.push_back(z);


        }
      }


    float y_com = accumulate(begin(ytemp), end(ytemp), 0.0);
    float z_com = accumulate(begin(ztemp), end(ztemp), 0.0);
    y_com = y_com / float(ytemp.size());
    z_com = z_com / float(ztemp.size());

    this->z_com_frames += z_com;


    for (int i=0; i<ytemp.size(); i++) {
        y = ytemp[i];
        z = ztemp[i];
        iybin = int((y - y_com + yshift)/this->dr);
        izbin = int((z + this->zshift)/this->dr);
        if (izbin >=0 && izbin < this->zbins && iybin >=0 && iybin < ybins) {
            density_yz[iybin][izbin] += 1.0;
            if (z > -this->zshift && (z+zshift)*(z+zshift) + (y - y_com)*(y - y_com)< R_dense*R_dense ) {
                this->density_bulk += 1.0;
            }
        }
    }

    //cout << "number of atoms in sel1: " << sel1->segments_ind[0].size() << endl;

    for (auto &segment:sel1->segments_ind) {
	for (int ind : segment) {
	    z = system->z[ind];
            izbin = int((z+z_shift)/dz);
            if (izbin >=0 and izbin < this->nbins) {
                if (which_density_profile == "mass") {
                    profile[izbin] += system->atoms[ind].mass; 
                } else if (which_density_profile == "number") {
                    profile[izbin] += 1.0; 
                }
             }
	}
    }



    if (abs(box[0])<1.0e-6 or abs(box[1])<1.0e-6 or abs(box[2])<1.0e-6) {
       cout << "Error: box dimensions are not given!" << endl;
       exit(1);
    }

    float scaling;

    if (which_density_profile == "mass") {
        scaling = 1.0/(box[0]*box[1]*dz*Avogadro*pow(Atocm,3.0));
    } else if (which_density_profile == "number") {
        scaling = 1.0/(box[0]*box[1]*dz*1.0e-3);
    }


    for (int i = 0; i < this->nbins; i++) {
        profile[i] *= scaling;
    }


    return profile;
}


ANALYSIS_DENSITY_PROFILE_DROPLET::~ANALYSIS_DENSITY_PROFILE_DROPLET()
{
    system = NULL;
    sel1 = NULL;
}

