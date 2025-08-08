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
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "analysis_rotate.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_ROTATE::ANALYSIS_ROTATE(PSF *system, GROUP *sel1, int vector1d, int vector2d,int voidf, float bl, string name1, string name2, string name3, string name4, string name_ref1, string name_ref2, string filename,int nbinsangle,int axisid)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->name1 = name1;
    this->name2 = name2;
    this->name3 = name3;
    this->name4 = name4;
    this->name_ref1 = name_ref1;
    this->name_ref2 = name_ref2;
    cout << "name_ref1: " << this->name_ref1 << "; name_ref2: " << this->name_ref2 << endl; //debug
    this->bl = bl;
    this->nbinsangle = nbinsangle;
    this->axisid = axisid;
}

void ANALYSIS_ROTATE::init() {
}





void ANALYSIS_ROTATE::compute_void() {
    FILE * fpdb;
    int ind_a1;
    int ind_a2;
    int ind_a3;
    int ind_a4;
    int ind_H;
    int ind_a1_begin;
    int ind_a1_second;
    int ind_a2_begin;
    int ind_a3_begin;
    int ind_a4_begin;
    int ind_H_begin;
    int ind_a4_ref;
    int ind_a3_pre;
    int ind_a4_pre;
    int ind_a1_pre;
    int ind_H_pre;
    int ind_a1_end;
    int ind_ref_H1;
    int ind_ref_H2;
    int ind_ref_a1;
    int ind_head;
    vector<float> shift(3);
    vector<float> vec1(3);
    vector<float> vec2(3);
    vector<float> vec_CC(3);
    vector<float> pos_pivot(3);
    vector<float> dist_rot(3);
  //  vector<float> next_bond(3);

    cout << "Nresidues: " << sel1->residues.size() << endl;


/////////////////////////////end determining the C-C bond///////////////////////

    for (int ind : sel1->residues[1]) {
	if (system->atomname[ind] == name1) {
             ind_a1_begin = ind;
	 }
         if (system->atomname[ind] == name2) {
	     ind_a2_begin = ind;
	 }
         if (system->atomname[ind] == name3) {
	     ind_a3_begin = ind;
	 }
         if (system->atomname[ind] == name4) {
	     ind_a4_begin = ind;
	 }
      //       cout << "testing C1-4: " << ind_a1 << " " << ind_a4 << endl;
         if (system->atomname[ind] == name_ref2) {
	     ind_H_begin = ind;
	 }
    }

    for (int ind : sel1->residues[2]) {
	if (system->atomname[ind] == name1) {
             ind_a1_second = ind;
	 }
    }

    float dangle = 2.0 * PI / float(nbinsangle);

    for (int ibin = 0; ibin < nbinsangle; ibin++) {
        string pdbfilename = filename + to_string(ibin) + ".pdb";
        const char *filechar = pdbfilename.c_str();
        fpdb = fopen(filechar,"w");
    	int restempid = 2;
        if (axisid == 1) {
            residue_rotate_random(ind_a3_begin, ind_a4_begin, ind_a4_begin, dangle, restempid, ind_a4_begin, ind_a2); //restempid means the first residue to rotate
        } else if (axisid == 2) {
            residue_rotate_random_regular(ind_a4_begin, ind_a1_second, ind_a1_second, dangle, restempid, ind_a2); //restempid means the first residue to rotate
        }
        for (int restempid = 3; restempid < sel1->residues.size(); restempid++) {
            if (axisid == 1) {
                residue_rotate_random_regular(ind_a3_begin, ind_a4_begin, ind_a4_begin, dangle, restempid, ind_a2);
	    } else if (axisid == 2) {
                residue_rotate_random_regular(ind_a4_begin, ind_a1_second, ind_a1_second, dangle, restempid, ind_a2);
	    }
        }
  

        for (int restempid = 0; restempid < sel1->residues.size(); restempid++) {
            for (int ind : sel1->residues[restempid]) {
                fprintf(fpdb,"%4s%7d  %-4s%-4s%1d%4d%12.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM", ind + 1, system->atomname[ind].c_str(), system->resname[ind].c_str(), 1, system->resid[ind], system->x[ind], system->y[ind], system->z[ind], 0,0);
            }
        }
        fclose(fpdb);    
    }

}


ANALYSIS_ROTATE::~ANALYSIS_ROTATE()
{
    system = NULL;
    sel1 = NULL;
}

