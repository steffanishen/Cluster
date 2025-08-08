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
#include "analysis_energy.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_ENERGY::ANALYSIS_ENERGY(PSF *system, GROUP *sel1, GROUP *sel2, int vector1d, int vector2d, int voidf, string filename, float dist_crit)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->dist_crit = dist_crit;
}

void ANALYSIS_ENERGY::init() {
}

vector<float> ANALYSIS_ENERGY::compute_vector() {
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    float zave = 0.0;
    float energy_single = 0.0;

    vector<float> energy;
    vector<double> box;
    box.resize(3);
    box[0] = system->box_first_frame[0];
    box[1] = system->box_first_frame[1]; 
    box[2] = system->box_first_frame[2];
    //box[0] = system->pbc[0];
    //box[1] = system->pbc[2];
    //box[2] = system->pbc[5];
    float shift = box[2] * 0.5;

    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    if (sel2->NATOM == 0) error1.error_exit("ERROR: sel2 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose
//    cout << "M_PI" << M_PI << endl;// for debug purpose
 
    for (auto &segment:sel1->segments_ind) {
	for (int ind : segment) {
	    r[0] = system->x[ind];
	    r[1] = system->y[ind];
	    r[2] = system->z[ind];
            zave += r[2];
            float q0 = system->charge[ind];
            //cout <<"ind: " << ind << " charge: " << q0 << endl; //for debug purpose

    //cout << "sel1->atom_index: " << system->atom_index[ind] << endl; //for debug purpose
    //cout << "sel1->atomtype: " << system->atomtype[ind] << endl; //for debug purpose

            for (auto &segment:sel2->segments_ind) {
	        for (int ind1 : segment) {
	            r1[0] = system->x[ind1];
	            r1[1] = system->y[ind1];
	            r1[2] = system->z[ind1];
                    float q1 = system->charge[ind1];
                    
	            disp = getDistPoints(r, r1);
                    dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
                    dist = sqrt(dist2);
                    if (abs(dist) < 1e-6) error1.error_exit("ERROR: Atom overlapping! Check atom selection, sel1 should not overlap with sel2");
                    
                    energy_single += q0*q1/dist;
 
                }
            }
	}
    }

    zave = zave/float(sel1->NATOM); 

    float scaling = etoCoul*etoCoul/(4.0*M_PI*permittivity*1e-10*kBT);

    energy_single = energy_single * scaling;
  
    energy.push_back(zave);
    energy.push_back(energy_single);



    return energy;
}


ANALYSIS_ENERGY::~ANALYSIS_ENERGY()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
}

