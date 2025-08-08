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
//***************** use of this class must follow anglezs_rings **************************

#include <iostream>
#include <vector>
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include "analysis.hpp"
#ifndef ANALYSIS_ORIENTATION_CHAIN_HPP
#define	ANALYSIS_ORIENTATION_CHAIN_HPP

using namespace std;

class  ANALYSIS_ORIENTATION_CHAIN : public ANALYSIS
{

private:
    //private attributes
    //no private methods
    
public:
   
    // no public attributes
    // public methods
//    GROUP *sel1;
 //   int whichN;
    //
    int iframe;
    string name0;
    float dz;
    float zshift;
    float dtheta;
    int thetabins;
    int every_n_frame;
    float costheta2;
    vector<float> rdf_count_single_frame;
    vector<vector<float>> density_yz;
    vector<vector<float>> density_yz_phi;
    ofstream *density_file = new ofstream("density_theta_phi.dat");

    ANALYSIS_ORIENTATION_CHAIN(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename,  int nbins, int every_n_frame,float dtheta); //constructor
    
    void init();
    void output_density(vector<vector<float>> density_yz);

    vector<float> compute_vector();
     
    ~ANALYSIS_ORIENTATION_CHAIN();

};

#endif	/* DCD_R_HPP */

