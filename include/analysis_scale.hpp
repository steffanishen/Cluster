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
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include "analysis.hpp"
#ifndef ANALYSIS_SCALE_HPP
#define	ANALYSIS_SCALE_HPP

using namespace std;

class  ANALYSIS_SCALE : public ANALYSIS
{

private:
    //private attributes
    //no private methods
    
public:
   
    // no public attributes
    // public methods

    //
    ANALYSIS_SCALE(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, float scaling); //constructor
    
    void init();

    void compute_void();   
     
    ~ANALYSIS_SCALE();

};

#endif	/* DCD_R_HPP */

