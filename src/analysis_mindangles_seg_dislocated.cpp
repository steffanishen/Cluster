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
#include "analysis_mindangles_seg_dislocated.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_MINDANGLES_SEG_DISLOCATED::ANALYSIS_MINDANGLES_SEG_DISLOCATED(PSF *system, GROUP *sel1, int whichN,int vector1d,int vector2d,int voidf,string filename)
{
    this->system = system;
    this->sel1 = sel1;
    this->whichN = whichN;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
}

void ANALYSIS_MINDANGLES_SEG_DISLOCATED::init() {
}

vector<float> ANALYSIS_MINDANGLES_SEG_DISLOCATED::compute_vector() {
    vector<float> mindA_seg;
    mindA_seg.clear();
    if (sel1->anglezs.size() == 0) {
	error1.error_exit("Error: the anglezs are not computed yet!");
    }
    int istart = this->whichN;
    for (int segi = 0; segi < sel1->anglezs.size(); segi++) {
   	float mindA = 360.0;
   	    //choose the first N atom
   	    //calculate the smallest angle
   	int Nangles = sel1->anglezs[segi].size();
	int segi_pre;
	if (segi == 0) {
	    segi_pre = Nangles - 1;
   	} else {
	    segi_pre = segi - 1;
	}

	for (int Ni = 0; Ni < Nangles; Ni++) {
	    float dA = dangle(sel1->anglezs[segi][Ni], sel1->anglezs[segi_pre][istart]);
	    if (istart > 11) cerr << "istart: "<< istart << " >11!!!!!!!!!!!";
	    if (fabs(dA) < fabs(mindA)) {
		mindA = dA;
		istart = Ni;
	    }
	}

   	     //accumulate the angle
//   	    cout << "Minimum dangle: " << mindA << endl;
	mindA_seg.push_back(mindA);
//	    angle_file << mindA << " ";
    }
    return mindA_seg;
}


ANALYSIS_MINDANGLES_SEG_DISLOCATED::~ANALYSIS_MINDANGLES_SEG_DISLOCATED()
{
    system = NULL;
    sel1 = NULL;
}

