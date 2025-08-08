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
#include "analysis_shell.hpp"
#include <unistd.h>

#define PI 3.14159265

using namespace std;

ANALYSIS_SHELL::ANALYSIS_SHELL(PSF *system, GROUP *sel1, GROUP *sel2, int vector1d, int vector2d, int voidf, string filename, int nbins,float dist_crit)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->nbins = nbins;
    this->dist_crit = dist_crit;
}

void ANALYSIS_SHELL::init() {
}

vector<vector<float>> ANALYSIS_SHELL::compute_2D_vector() {
    sel1->anglezs.clear();
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;

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
    vector<int> atomids = sel1->atomids;
    vector<int> atomids1 = sel2->atomids;

   // cout << "atomids size:" << sel1->atomids.size() << endl;
 
    vector<vector<float>> nshell(atomids.size(),vector<float>(2,0.0));


        int THREADSIZE = atoi(std::getenv("OMP_NUM_THREADS"));
        #pragma omp parallel default(shared) num_threads(THREADSIZE)
        #pragma omp for schedule(dynamic) 
 
 
       
	for (int ind_sel=0; ind_sel < atomids.size(); ind_sel++) {
            int ind = atomids[ind_sel];
            vector<float> r(3,0.0);
	    r[0] = system->x[ind];
	    r[1] = system->y[ind];
	    r[2] = system->z[ind];
            float nshell_single = 0.0;
	    for (int ind_sel1=0; ind_sel1 < atomids1.size(); ind_sel1++) {
                    int ind1 = atomids1[ind_sel1]; 
                //    std::cout << "Current thread number: " << omp_get_thread_num() << std::endl;
                    vector<float> r1(3,0.0);

	            r1[0] = system->x[ind1];
	            r1[1] = system->y[ind1];
	            r1[2] = system->z[ind1];
                    
                    vector<float> dist(3,0.0);
                    float dist2;
	            dist = getDistPoints(r, r1);
                    dist2 = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];

                    if (dist2 < dist_crit2) {
                        nshell_single += 1.0; 
                    }
                }
                nshell[ind_sel][0] = r[2];
                nshell[ind_sel][1] = nshell_single;
            }
           // vector<float> nshell_temp;
            //nshell_temp.push_back(r[2]);
            //nshell_temp.push_back(nshell_single);



    if (abs(box[0])<1.0e-6 or abs(box[1])<1.0e-6 or abs(box[2])<1.0e-6) {
       cout << "Error: box dimensions are not given!" << endl;
       exit(1);
    }


    return nshell;
}


ANALYSIS_SHELL::~ANALYSIS_SHELL()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
}
