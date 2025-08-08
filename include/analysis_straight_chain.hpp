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
//***************** Function: rotate residues ***********************************************

#include <iostream>
#include <vector>
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include <stdio.h>
#include "analysis.hpp"
#ifndef ANALYSIS_STRAIGHT_CHAIN_HPP
#define	ANALYSIS_STRAIGHT_CHAIN_HPP

using namespace std;

class  ANALYSIS_STRAIGHT_CHAIN: public ANALYSIS
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
    ANALYSIS_STRAIGHT_CHAIN(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, float bl, string filename, vector<PSF*> monomers, float fringe); //constructor
    
    void init();

    vector<vector<float>> bond_next(int residueid, PSF *monomer);
    vector<vector<float>> bond_pre(int residueid, PSF *monomer);

    void compute_void();
     
    ~ANALYSIS_STRAIGHT_CHAIN();

};

#endif	/* DCD_R_HPP */

