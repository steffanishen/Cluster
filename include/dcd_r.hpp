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

#include "dcd.hpp"
#include <iostream>
#ifndef DCD_R_HPP
#define	DCD_R_HPP

using namespace std;

class DCD_R : public DCD
{

private:
    //no private attributes
    //private methods
    void alloc();
    
public:
    
    // no public attributes
    // public methods
    DCD_R(string filename, PSF* sys); //constructor
    
    void read_header();
    void read_oneFrame();
    void printHeader() const;


    vector<float> getDist(ATOM a1, ATOM a2);
    vector<float> getDistPoints(vector<float> p1, vector<float> p2);
    void dcdrdebug();
    float angle_z(vector<float> dist);
    float dangle(float a1, float a2);
        
    ~DCD_R();

};

#endif	/* DCD_R_HPP */

