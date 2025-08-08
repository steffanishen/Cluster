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


#ifndef OUTPUT_HPP_INCLUDED
#define OUTPUT_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "analysis.hpp"
#include "error1.hpp"

using namespace std;

class OUTPUT
{
public:
    vector<string> filenames; 
    vector<ofstream> output_files;  // for opening input file


    ERROR1 error1;

    int whichN;

    OUTPUT();
   

    ~OUTPUT();
};

#endif // DCD_HPP_INCLUDED
