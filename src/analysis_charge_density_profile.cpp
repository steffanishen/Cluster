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

#include <math.h>
#include "analysis_charge_density_profile.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_CHARGE_DENSITY_PROFILE::ANALYSIS_CHARGE_DENSITY_PROFILE(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, int nbins,string which_density_profile)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->nbins = nbins;
    this->which_density_profile = which_density_profile;
}

void ANALYSIS_CHARGE_DENSITY_PROFILE::init() {
}

vector<float> ANALYSIS_CHARGE_DENSITY_PROFILE::compute_vector() {
    sel1->anglezs.clear();
    float z;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;

    vector<float> profile(this->nbins,0.0);
    vector<double> box;
    box.resize(3);
    box[0] = system->box_first_frame[0];
    box[1] = system->box_first_frame[1]; 
    box[2] = system->box_first_frame[2];
    //box[0] = system->pbc[0];
    //box[1] = system->pbc[2];
    //box[2] = system->pbc[5];
    float dz = box[2]/float(this->nbins);
    float shift = box[2] * 0.5;

    //cout << "number of atoms in sel1: " << sel1->segments_ind[0].size() << endl;

    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
            if (system->z[ind] < -shift) {
	            z = system->z[ind] + box[2] + shift;
            } else if (system->z[ind] >= shift) {
	            z = system->z[ind] - box[2] + shift;
            } else {
	            z = system->z[ind] + shift;
            }
            izbin = int((z)/dz);

            if (izbin == this->nbins) izbin -= 1;
        // Meng comment: debug
            if (izbin < 0 or izbin >= this->nbins) {
                  cout << "z position is: " << system->z[ind] << "; lz is: " << box[2] << "; shift is: " << shift << "; izbin is: " << izbin << endl;
        //        exit(0);
            }
        // Meng comment: end of debug


            profile[izbin] += system->atoms[ind].charge; 
	    }
    }


    if (abs(box[0])<1.0e-6 or abs(box[1])<1.0e-6 or abs(box[2])<1.0e-6) {
       cout << "Error: box dimensions are not given!" << endl;
       exit(1);
    }

    float scaling;

    scaling = 1.0/(box[0]*box[1]*dz*1.0e-3);


    for (int i = 0; i < this->nbins; i++) {
        profile[i] *= scaling;
    }


    return profile;
}


ANALYSIS_CHARGE_DENSITY_PROFILE::~ANALYSIS_CHARGE_DENSITY_PROFILE()
{
    system = NULL;
    sel1 = NULL;
}

