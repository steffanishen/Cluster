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
#include "analysis_mindangles_seg.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_MINDANGLES_SEG::ANALYSIS_MINDANGLES_SEG(PSF *system, GROUP *sel1, int vector1d, int vector2d,int voidf,string filename)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
}

void ANALYSIS_MINDANGLES_SEG::init() {
}

vector<float> ANALYSIS_MINDANGLES_SEG::compute_vector() {
    vector<float> mindA_seg;
    mindA_seg.clear();
    vector<float> danglezsi;
    danglezs.clear();
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

	danglezsi.clear();
	for (int Ni = 0; Ni < Nangles; Ni++) {
	    float mindAi = 360.0;
	    for (int Nj = 0; Nj < Nangles; Nj++) {
	    	float dA = dangle(sel1->anglezs[segi][Nj], sel1->anglezs[segi_pre][Ni]);
	    	if (fabs(dA) < fabs(mindAi)) {
		    mindAi = dA;
	    	}
	    }
	    danglezsi.push_back(mindAi);
	    if (fabs(mindAi) < fabs(mindA)) {
		mindA = mindAi;
	    }
	}
	danglezs.push_back(danglezsi);

   	     //accumulate the angle
//   	    cout << "Minimum dangle: " << mindA << endl;
	mindA_seg.push_back(mindA);
//	    angle_file << mindA << " ";
    }
    return mindA_seg;
}


ANALYSIS_MINDANGLES_SEG::~ANALYSIS_MINDANGLES_SEG()
{
    system = NULL;
    sel1 = NULL;
}

