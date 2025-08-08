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
#include "analysis_trajectory_extract.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_TRAJECTORY_EXTRACT::ANALYSIS_TRAJECTORY_EXTRACT(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
}

void ANALYSIS_TRAJECTORY_EXTRACT::init() {
}

vector<float> ANALYSIS_TRAJECTORY_EXTRACT::compute_vector() {
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r_ave(4,0.0);
    vector<float> disp(3,0.0);
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    this->iframe += 1;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;


    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose
//    cout << "M_PI" << M_PI << endl;// for debug purpose
 
    r_ave = {float(this->iframe),0.0,0.0,0.0};
    for (auto &segment:sel1->segments_ind) {
    	for (int ind : segment) {
    	    r[0] = system->x[ind];
    	    r[1] = system->y[ind];
    	    r[2] = system->z[ind];

            for (int dim = 0; dim < 3; dim++) {
                r_ave[dim+1] += r[dim];
            }
    //cout << "sel1->atom_index: " << system->atom_index[ind] << endl; //for debug purpose
    //cout << "sel1->atomtype: " << system->atomtype[ind] << endl; //for debug purpose

    

    	}
    }

    for (int dim = 1; dim <= 3; dim++) {
        r_ave[dim] /= float(sel1->NATOM);
    }

   cout << sel1->NATOM << endl;

    return r_ave;
}


ANALYSIS_TRAJECTORY_EXTRACT::~ANALYSIS_TRAJECTORY_EXTRACT()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
    this->rdf_count.clear();
}

