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

#include "input.hpp"

using namespace std;

INPUT::INPUT(string filename)
{
    
    input_script.exceptions(std::ifstream::failbit);
    try
    {
        input_script.open(filename,ios::in);
    }
    catch(std::ifstream::failure e)
    {
        cerr << "Exception opening/reading input script '" << filename << "' : " << std::endl;
        cerr << "Please check the path of the input script and if it exists." << endl;
    } 
    
}



void INPUT::read_parameters()
{
    if (!input_script.good()) this->error1.error_exit( "No input script!.");
    string lineStr; 
    string keyword;
    while(getline(input_script,lineStr)) {
        istringstream iss(lineStr);
	iss >> keyword;
	if (keyword == "End_of_header") break;
	cout << "keyword: " << keyword << endl;
    }
    this->monomers.clear();
    this->dcdnames.clear();

    if (keyword != "End_of_header") cerr << "Please check if there is an End_of_header keyword." << endl;

    while (getline (input_script,lineStr)) {
         istringstream iss(lineStr);
	 string keyword;
         string temp;
	 string monomer;
	 string dcdname;
	 vector<string> arg1;
	 int int_param;
	 float int_float;
         iss >> keyword;
	 if (keyword == "DCDFILE") {
	     iss >> dcdname;
	     this->dcdnames.push_back(dcdname);
	     cout << "keyword: " << dcdname << endl;
	 } else if (keyword == "PSFFILE") {
	     iss >> this->psfname;
	     cout << "keyword: " << psfname << endl;
	 } else if (keyword == "MONOMER") {
	     iss >> monomer;
	     this->monomers.push_back(monomer);
	 } else if (keyword == "GROUP") {
	     arg1.clear();
 	     while (iss >> temp) {
	         arg1.push_back(temp);
  	     }
	     this->group_selection.push_back(arg1);
	     arg1.clear();
	 } else if (keyword == "ANALYSIS") {
	     arg1.clear();
 	     while (iss >> temp) {
	         arg1.push_back(temp);
  	     }
	     this->analysis.push_back(arg1);
	     arg1.clear();
	 } else if (keyword == "whichN") {
	     iss >> this->whichN;
	     cout << "keyword: " << psfname << endl;
	 }
	 if (keyword == "End_of_file") break;
    }

}



INPUT::~INPUT()
{
    input_script.close();
    
}
