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

#ifndef ERROR1_HPP_INCLUDED
#define ERROR1_HPP_INCLUDED

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <stack>
#include <stdio.h>
#include <ctype.h>

using namespace std;

class ERROR1
{
protected:
    //protected attributes
     
public:


//    keywords_sel.insert(make_pair("resname",resname));

    ERROR1(){
    }

    void error_exit(string str);
    //ATOM(int atom_index = 0, vector<float> coord={0.0, 0.0, 0.0}, int segid = 0, int resid = 0, string resname = "NULL", string atomname = "NULL", string atomtype = "NULL", float charge = 0.0, float mass = 0.0, int beta = 0);
    
    
    ~ERROR1();
};

#endif // DCD_HPP_INCLUDED
