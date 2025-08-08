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

#ifndef ATOM_HPP_INCLUDED
#define ATOM_HPP_INCLUDED

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <stack>
#include <stdio.h>
#include <ctype.h>

using namespace std;

class ATOM
{
protected:
    //protected attributes
     
public:
    int atom_index; // Array storing indexes of moving atoms. Note the indices start from 0.
    vector<float> coord{0.0, 0.0, 0.0};

    int segid;
    string segname;
    int resid;
    string resname;
    string atomname;
    string atomtype;
    float charge;
    float mass;
    int beta;
    int atom_flag = 0;
    int debug_flag = 0;
    map<string,string> keywords_sel;

//    keywords_sel.insert(make_pair("resname",resname));

    ATOM(){
    }

    int string_comp(vector<string> str);
    int string_comp_atom_index(vector<string> str);
    int string_comp_segid(vector<string> str);
    int string_comp_resid(vector<string> str);
    int string_comp_resname(vector<string> str);
    int string_comp_segname(vector<string> str);
    int string_comp_atomname(vector<string> str);
    int string_comp_atomtype(vector<string> str);
    int string_comp_charge(vector<string> str);
    int string_comp_mass(vector<string> str);
    int string_comp_beta(vector<string> str);

    int select_atoms(vector<string> str);

    bool is_selection(string str);
    bool is_operator(string str);

    bool has_precedence(string op1, string op2);
    int  op_apply(string op1, int val1, int val2);

    int ind_comp(vector<string> str);

    int arithmetic_calc(vector<string> str);
    void error_exit(string str);
    //ATOM(int atom_index = 0, vector<float> coord={0.0, 0.0, 0.0}, int segid = 0, int resid = 0, string resname = "NULL", string atomname = "NULL", string atomtype = "NULL", float charge = 0.0, float mass = 0.0, int beta = 0);
    
    
    ~ATOM();
};

#endif // DCD_HPP_INCLUDED
