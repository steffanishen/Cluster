/*
 *  read_dcd : c++ class + main file example for reading a CHARMM dcd file
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//***************** Partially contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

#include "group.hpp"

#define PI 3.14159265

using namespace std;

GROUP::GROUP()
{
    anglezs.clear();
}


GROUP::~GROUP()
{
    atoms.clear();
    residues.clear();
    segments.clear();
    segments_ind.clear();
    atomids.clear();
    anglezs.clear();
}
