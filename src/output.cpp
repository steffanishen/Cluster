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

//***************** Fully contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************
//****************** Note the x, y, z and pbc need to be carefully considered, currently we assume they are deallocated simultaneously with thathat in DCD when being used, but we'll work on this later******************************************

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "output.hpp"

using namespace std;

OUTPUT::OUTPUT()
{
    
}




OUTPUT::~OUTPUT()
{
    
}
