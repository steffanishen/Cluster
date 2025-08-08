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
#include "analysis_msd_region.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_MSD_REGION::ANALYSIS_MSD_REGION(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, int dtmax,int msd_xflag,int msd_yflag,int msd_zflag,float low,float high)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->dtmax = dtmax;
    this->msd_previous.resize(dtmax); //debug
    fill(this->msd_previous.begin(), this->msd_previous.end(),0.0);
    this->count.resize(dtmax); //debug
    fill(this->count.begin(), this->count.end(),0.0);
    this->msd_xflag = msd_xflag;
    this->msd_yflag = msd_yflag;
    this->msd_zflag = msd_zflag;
    this->low = low;
    this->high = high;
}

void ANALYSIS_MSD_REGION::init() {
}

vector<float> ANALYSIS_MSD_REGION::compute_vector() {
    vector<float> msd(this->dtmax,0.0);
    if (this->msd_xflag == 1) {
         msd = msd_region_x();
    } else if (this->msd_yflag == 1) {
         msd = msd_region_y();
    } else if (this->msd_zflag == 1) {
         msd = msd_region_z();
    } else {
         msd = msd_region();
    }
    return msd;
}

vector<float> ANALYSIS_MSD_REGION::msd_region() {
    sel1->anglezs.clear();
    float x;
    float y;
    float z;
    float x0;
    float y0;
    float z0;
    float dx;
    float dy;
    float dz;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;

    vector<float> msd(this->dtmax,0.0);
    vector<double> box;
    box.resize(3);

    int frame_stored = system->XS.size();
    int dtmax_current = min(this->dtmax,frame_stored);

 

    if (system->iframe > 0) {


        for (int i=0; i< dtmax_current; i++) {
            int dt = dtmax_current - i;
            for (auto &segment:sel1->segments_ind) {
	        for (int ind : segment) {
	            x = system->x[ind];
	            y = system->y[ind];
	            z = system->z[ind];
	            x0 = system->XS[i][ind];
	            y0 = system->YS[i][ind];
	            z0 = system->ZS[i][ind];
 
		    dx = x-x0;
                    dy = y-y0;
                    dz = z-z0;


                    this->msd_previous[dt] += dx*dx + dy*dy + dz*dz;
                    this->count[dt] += 1.0;
                }
	    }
        }

        for (int i=0; i< this->dtmax; i++) {
            if (this->count[i] > 0.0) {
                msd[i] = this->msd_previous[i]/this->count[i];
            }
        }



    }

    system->XS.push_back(new float[system->NATOM]);
    system->YS.push_back(new float[system->NATOM]);
    system->ZS.push_back(new float[system->NATOM]);

    memcpy(system->XS[frame_stored],system->x,system->NATOM*sizeof(float));
    memcpy(system->YS[frame_stored],system->y,system->NATOM*sizeof(float));
    memcpy(system->ZS[frame_stored],system->z,system->NATOM*sizeof(float));



    if (system->XS.size() >= this->dtmax) {
        system->XS.erase(system->XS.begin());
        system->YS.erase(system->YS.begin());
        system->ZS.erase(system->ZS.begin());
            
     }

    return msd;
}


vector<float> ANALYSIS_MSD_REGION::msd_region_x() {
    sel1->anglezs.clear();
    float x;
    float y;
    float z;
    float x0;
    float y0;
    float z0;
    float dx;
    float dy;
    float dz;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;

    vector<float> msd(this->dtmax,0.0);
    vector<double> box;
    box.resize(3);

    int frame_stored = system->XS.size();
    int dtmax_current = min(this->dtmax,frame_stored);

 
       for (auto &segment:sel1->segments_ind) {
 	   for (int ind : segment) {
	       x = system->x[ind];
	       y = system->y[ind];
	       z = system->z[ind];
               int tmin_msd_old = system->tmin_msd[ind];

           if (tmin_msd_old < 0) {
               if (x > low && x < high) {
                   system->tmin_msd[ind] = dtmax_current;
               }
           }

           if (tmin_msd_old >= 0) {
               
               if (x > low && x < high) {

                   for (int i=system->tmin_msd[ind]; i< dtmax_current; i++) {
                        int dt = dtmax_current - i;
	                x0 = system->XS[i][ind];
	                y0 = system->YS[i][ind];
                        z0 = system->ZS[i][ind];
 
		        dx = x-x0;
                        dy = y-y0;
                        dz = z-z0;

                        this->msd_previous[dt] += dx*dx + dy*dy + dz*dz;
                        this->count[dt] += 1.0;
                    }
                    if (dtmax_current == this->dtmax) system->tmin_msd[ind] -= 1;
                }
	   }

           if (x < low || x > high) {
               system->tmin_msd[ind] = -1;
           }
        }
    }

    for (int i=0; i< this->dtmax; i++) {
        if (this->count[i] > 0.0) {
            msd[i] = this->msd_previous[i]/this->count[i];
        }
    }


    system->XS.push_back(new float[system->NATOM]);
    system->YS.push_back(new float[system->NATOM]);
    system->ZS.push_back(new float[system->NATOM]);

    memcpy(system->XS[frame_stored],system->x,system->NATOM*sizeof(float));
    memcpy(system->YS[frame_stored],system->y,system->NATOM*sizeof(float));
    memcpy(system->ZS[frame_stored],system->z,system->NATOM*sizeof(float));


    if (system->XS.size() >= this->dtmax) {
        system->XS.erase(system->XS.begin());
        system->YS.erase(system->YS.begin());
        system->ZS.erase(system->ZS.begin());
            
     }

    return msd;
}

vector<float> ANALYSIS_MSD_REGION::msd_region_y() {
    sel1->anglezs.clear();
    float x;
    float y;
    float z;
    float x0;
    float y0;
    float z0;
    float dx;
    float dy;
    float dz;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;

    vector<float> msd(this->dtmax,0.0);
    vector<double> box;
    box.resize(3);

    int frame_stored = system->XS.size();
    int dtmax_current = min(this->dtmax,frame_stored);

 
       for (auto &segment:sel1->segments_ind) {
 	   for (int ind : segment) {
	       x = system->x[ind];
	       y = system->y[ind];
	       z = system->z[ind];
               int tmin_msd_old = system->tmin_msd[ind];

           if (tmin_msd_old < 0) {
               if (y > low && y < high) {
                   system->tmin_msd[ind] = dtmax_current;
               }
           }

           if (tmin_msd_old >= 0) {
               
               if (y > low && y < high) {

                   for (int i=system->tmin_msd[ind]; i< dtmax_current; i++) {
                        int dt = dtmax_current - i;
	                x0 = system->XS[i][ind];
	                y0 = system->YS[i][ind];
                        z0 = system->ZS[i][ind];
 
		        dx = x-x0;
                        dy = y-y0;
                        dz = z-z0;

                        this->msd_previous[dt] += dx*dx + dy*dy + dz*dz;
                        this->count[dt] += 1.0;
                    }
                    if (dtmax_current == this->dtmax) system->tmin_msd[ind] -= 1;
                }
	   }

           if (y < low || y > high) {
               system->tmin_msd[ind] = -1;
           }
        }
    }

    for (int i=0; i< this->dtmax; i++) {
        if (this->count[i] > 0.0) {
            msd[i] = this->msd_previous[i]/this->count[i];
        }
    }


    system->XS.push_back(new float[system->NATOM]);
    system->YS.push_back(new float[system->NATOM]);
    system->ZS.push_back(new float[system->NATOM]);

    memcpy(system->XS[frame_stored],system->x,system->NATOM*sizeof(float));
    memcpy(system->YS[frame_stored],system->y,system->NATOM*sizeof(float));
    memcpy(system->ZS[frame_stored],system->z,system->NATOM*sizeof(float));


    if (system->XS.size() >= this->dtmax) {
        system->XS.erase(system->XS.begin());
        system->YS.erase(system->YS.begin());
        system->ZS.erase(system->ZS.begin());
            
     }

    return msd;
}


vector<float> ANALYSIS_MSD_REGION::msd_region_z() {
    sel1->anglezs.clear();
    float x;
    float y;
    float z;
    float x0;
    float y0;
    float z0;
    float dx;
    float dy;
    float dz;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;

    vector<float> msd(this->dtmax,0.0);
    vector<double> box;
    box.resize(3);

    int frame_stored = system->XS.size();
    int dtmax_current = min(this->dtmax,frame_stored);

 
       for (auto &segment:sel1->segments_ind) {
 	   for (int ind : segment) {
	       x = system->x[ind];
	       y = system->y[ind];
	       z = system->z[ind];
               int tmin_msd_old = system->tmin_msd[ind];

           if (tmin_msd_old < 0) {
               if (z > low && z < high) {
                   system->tmin_msd[ind] = dtmax_current;
               }
           }

           if (tmin_msd_old >= 0) {
               
               if (z > low && z < high) {

                   for (int i=system->tmin_msd[ind]; i< dtmax_current; i++) {
                        int dt = dtmax_current - i;
	                x0 = system->XS[i][ind];
	                y0 = system->YS[i][ind];
                        z0 = system->ZS[i][ind];
 
		        dx = x-x0;
                        dy = y-y0;
                        dz = z-z0;

                        this->msd_previous[dt] += dx*dx + dy*dy + dz*dz;
                        this->count[dt] += 1.0;
                    }
                    if (dtmax_current == this->dtmax) system->tmin_msd[ind] -= 1;
                }
	   }

           if (z < low || z > high) {
               system->tmin_msd[ind] = -1;
           }
        }
    }

    for (int i=0; i< this->dtmax; i++) {
        if (this->count[i] > 0.0) {
            msd[i] = this->msd_previous[i]/this->count[i];
        }
    }


    system->XS.push_back(new float[system->NATOM]);
    system->YS.push_back(new float[system->NATOM]);
    system->ZS.push_back(new float[system->NATOM]);

    memcpy(system->XS[frame_stored],system->x,system->NATOM*sizeof(float));
    memcpy(system->YS[frame_stored],system->y,system->NATOM*sizeof(float));
    memcpy(system->ZS[frame_stored],system->z,system->NATOM*sizeof(float));


    if (system->XS.size() >= this->dtmax) {
        system->XS.erase(system->XS.begin());
        system->YS.erase(system->YS.begin());
        system->ZS.erase(system->ZS.begin());
            
     }

    return msd;
}


ANALYSIS_MSD_REGION::~ANALYSIS_MSD_REGION()
{
    system = NULL;
    sel1 = NULL;
}

