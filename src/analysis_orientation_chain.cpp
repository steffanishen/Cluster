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
#include "analysis_orientation_chain.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_ORIENTATION_CHAIN::ANALYSIS_ORIENTATION_CHAIN(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, int nbins, int every_n_frame,float dtheta)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->name0 = name0;
    this->name1 = name1;
    this->name2 = name2;
    this->name3 = name3;
    this->dist_crit = dist_crit;
    this->nbins = nbins;
//    this->rdf_count.resize(nbins);
//    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
    this->dtheta = dtheta;
    this->every_n_frame = every_n_frame;


}

void ANALYSIS_ORIENTATION_CHAIN::init() {
}

void ANALYSIS_ORIENTATION_CHAIN::output_density(vector<vector<float>> density_yz) {
    for (int ithetabin=0; ithetabin < this->thetabins; ithetabin++) {
        for (int iphibin=0; iphibin < this->thetabins; iphibin++) {
            float theta_contour = this->dtheta * float(ithetabin);
            float phi_contour = this->dtheta * float(iphibin);
            *this->density_file << theta_contour << " " << phi_contour << " " << density_yz[ithetabin][iphibin] << endl;
        }
    }
}

vector<float> ANALYSIS_ORIENTATION_CHAIN::compute_vector() {
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> r2(3,0.0);
    vector<float> r3(3,0.0);
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
    this->iframe += 1;
    float costheta = 0.0;
    float cosphi = 0.0;
    vector<double> box;
    vector<float> order_parameters(2);

    box.resize(3);
    if (this->iframe == 1){
        box[2] = system->box_first_frame[2];
        this->zshift = box[2] * 0.5;
        this->dz = box[2]/float(this->nbins);
        this->thetabins = int (180.0 / this->dtheta);
        this->density_yz.resize(this->thetabins,vector<float>(this->thetabins*2));
        this->costheta2 = 0.0;
    }



    float order_parameter;

    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose

    int nchains = sel1->segments.size();



    for (auto &segment:sel1->segments) {
        r = residue_com(segment[0]);
        int nres = segment.size();
        r1 = residue_com(segment[nres-1]);
        disp = getDistPoints(r,r1);


        vector<float> plane_norm_xy_projection(3);
        plane_norm_xy_projection[0] = disp[0];
        plane_norm_xy_projection[1] = disp[1];
        plane_norm_xy_projection[2] = 0.0;

        vector<float> nz = {0.0, 0.0, -1.0};
        //costheta = dot_product(plane_norm,nz)/(norm(plane_norm)*norm(nz));

        vector<float> nx = {1.0, 0.0, 0.0};
        cosphi = dot_product(plane_norm_xy_projection,nx)/(norm(plane_norm_xy_projection)*norm(nx));

        float theta = acos (costheta) * 180.0 / PI;
        float phi = acos (cosphi) * 180.0 / PI;

        if (disp[1] < 0.0) {
            phi = 360.0 - phi;
        //    theta = PI - theta;
        }

        costheta = dot_product(disp,nz)/(norm(disp)*norm(nz));

	//cout << "costheta2: " << costheta2 << endl;

        this->costheta2 += costheta * costheta;
        int itheta = int(theta/this->dtheta);
        int iphi = int(phi/this->dtheta);

        if (itheta < this->thetabins && iphi < this->thetabins*2) density_yz[itheta][iphi] += 1.0;

    }


	this->costheta2 = this->costheta2 / float(nchains);
    order_parameter = (3.0 * this->costheta2 - 1.0) * 0.5;


    if (this->iframe == this->every_n_frame * (system->nframes_tot/this->every_n_frame)) {
        for (int ithetabin = 0; ithetabin < this->thetabins; ithetabin ++) {
            float theta_temp = this->dtheta /180.0 * PI * float(ithetabin);
            if (theta_temp < 1.0e-9) theta_temp = 1.0e-7;
	    float sin_theta = sin(theta_temp);
            cout << "sin_theta: " << sin_theta << endl;
            for (int iphibin = 0; iphibin < this->thetabins*2; iphibin ++) {
                if (density_yz[ithetabin][iphibin] > 1.0e-9) {
                density_yz[ithetabin][iphibin] = density_yz[ithetabin][iphibin]/(float(this->every_n_frame) * float(sel1->segments.size()) * this->dtheta * this->dtheta * sin_theta);
                }
           }
        }
        output_density(density_yz);
        this->density_yz.clear();
        this->density_yz.resize(this->nbins,vector<float>(this->thetabins,0.0));
    }

    order_parameters[0] = this->iframe;
    order_parameters[1] = order_parameter;

    return order_parameters;
}


ANALYSIS_ORIENTATION_CHAIN::~ANALYSIS_ORIENTATION_CHAIN()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
 //   this->rdf_count.clear();
}

