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

#include <iostream>
#include <vector>
#include "analysis.hpp"
#include "analysis_scale.hpp"
#include "analysis_anglezs_rings.hpp"
#include "analysis_mindangles_seg.hpp"
#include "analysis_avedangles_seg.hpp"
#include "analysis_mindangles_seg_dislocated.hpp"
#include "analysis_random_walk.hpp"
#include "analysis_straight_chain.hpp"
#include "analysis_rotate.hpp"
#include "analysis_packing.hpp"
#include "analysis_density_profile.hpp"
#include "analysis_charge_density_profile.hpp"
#include "analysis_density_profile_droplet.hpp"
#include "analysis_shell.hpp"
#include "analysis_shell_output_coord.hpp"
#include "analysis_energy.hpp"
#include "analysis_rdf.hpp"
#include "analysis_rdf_region.hpp"
#include "analysis_orientation.hpp"
#include "analysis_orientation_chain.hpp"
#include "analysis_functionalize.hpp"
#include "analysis_patch_no_order.hpp"
#include "analysis_trajectory_extract.hpp"
#include "analysis_passivation.hpp"
#include "analysis_msd.hpp"
#include "analysis_msd_region.hpp"
#include "analysis_contact_angle.hpp"
#include "analysis_contact_angle_density_profile.hpp"
#include "input.hpp"
#ifndef ANALYSIS_POINTERS_HPP
#define	ANALYSIS_POINTERS_HPP

using namespace std;

class ANALYSIS_POINTERS
{

private:
    //private attributes
    //no private methods
    
public:
   
    ERROR1 error1;
    // no public attributes
    // public methods
    PSF *system;
    vector<PSF*> monomers;
    vector<GROUP*> sels;
    vector<ANALYSIS*> analysis;
    INPUT* input;

    ANALYSIS_POINTERS(PSF *system, vector<GROUP*> sels, INPUT* input, vector<PSF*> monomers); //constructor
    
    vector<ANALYSIS*> init();
    //vector<vector<float>> anglezs_rings(GROUP *sel1);   
//    vector<float> mindangles_seg(GROUP *sel1, int whichN);
     
    virtual ~ANALYSIS_POINTERS();

};

#endif	/* DCD_R_HPP */

