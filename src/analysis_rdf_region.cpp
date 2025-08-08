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
#include "analysis_rdf_region.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_RDF_REGION::ANALYSIS_RDF_REGION(PSF *system, GROUP *sel1, GROUP *sel2, int vector1d, int vector2d, int voidf, string filename, float dist_crit, float dr,float zlow, float zhigh)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->dist_crit = dist_crit;
    this->dr = dr;
    nbins = int(this->dist_crit/dr);
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
    this->zlow = zlow;
    this->zhigh = zhigh;


}

void ANALYSIS_RDF_REGION::init() {
}

vector<float> ANALYSIS_RDF_REGION::compute_vector() {
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
    float dr = this->dr;
    int nbins = int(this->dist_crit/dr);
    this->iframe += 1;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;

    vector<float> rdf(nbins,0.0);

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

    //cout << "sel1->atom_index: " << system->atom_index[ind] << endl; //for debug purpose
    //cout << "sel1->atomtype: " << system->atomtype[ind] << endl; //for debug purpose
            if (r[2] > this->zlow && r[2] < this->zhigh) {
                for (auto &segment:sel2->segments_ind) {
	                for (int ind1 : segment) {
	                    r1[0] = system->x[ind1];
	                    r1[1] = system->y[ind1];
	                    r1[2] = system->z[ind1];
                    
                        if (r[2] > this->zlow && r[2] < this->zhigh) {

	                        disp = getDistPoints(r, r1);
                            dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
                            if (dist2 < dist_crit2) {
                                dist = sqrt(dist2);
                                int ibin = int(dist/dr);
                                this->rdf_count[ibin] += 1.0;
                            }
                        }
 
                    }
                }
            }


	}
    }


    for (int ibin = 0; ibin <= nbins; ibin++) {
        float drp = dr * ibin;
        float drn = dr * (ibin+1);
        rdf[ibin] = this->rdf_count[ibin]/(this->iframe * 4.0/3.0 * PI * (drn*drn*drn - drp*drp*drp));
 //       cout <<"rdf: " << rdf[ibin] << endl;
    }

    for (int ibin = 0; ibin <= nbins; ibin++) {
        float nideal = rdf[nbins-1];
        rdf[ibin] /= nideal; 
    }
/*
*/

    return rdf;
}


ANALYSIS_RDF_REGION::~ANALYSIS_RDF_REGION()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
    this->rdf_count.clear();
}

