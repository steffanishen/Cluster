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
//***************** use of this class must follow anglezs_rings **************************

#include <iostream>
#include <vector>
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include "analysis.hpp"
#include "cluster.hpp"
#ifndef ANALYSIS_PASSIVATION_HPP
#define	ANALYSIS_PASSIVATION_HPP

using namespace std;

class  ANALYSIS_PASSIVATION : public ANALYSIS
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
    int iframe;

    int nclusters = 0;
    int nresidues = 0;
    vector<int> linkedlist;
    ofstream *file_temp;
    ofstream *seg_res_ids; 
    vector<GROUP*> sels;
    vector<CLUSTER*> clusters;
    vector<vector<int>> residue_cluster_ind;
    int crosslink_in_cluster = 0;

    ifstream *input_cluster;
    ofstream *output_cluster;

    ANALYSIS_PASSIVATION(PSF *system, vector<GROUP*> sels, int vector1d, int vector2d, int voidf, string filename); //constructor
    
    void init();

    string patchtype(string name1, string resname1, string name2, string resname2  );
    string patchtype(string name1, string resname1 );
    bool is_empty(std::ifstream *pFile);
    void select_atoms(GROUP *atoms_select);

    void compute_void();
     
    ~ANALYSIS_PASSIVATION();

};

#endif	/* DCD_R_HPP */

