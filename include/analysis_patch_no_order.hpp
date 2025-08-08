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
#ifndef ANALYSIS_PATCH_NO_ORDER_HPP
#define	ANALYSIS_PATCH_NO_ORDER_HPP

using namespace std;

class  ANALYSIS_PATCH_NO_ORDER : public ANALYSIS
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
    int nchains = 0;
    int nresidues = 0;
    vector<int> linkedlist;
    ofstream *file_temp;
    ofstream *seg_res_ids; 
    vector<GROUP*> sels;
    vector<CLUSTER*> clusters;
    vector<CLUSTER*> chains;
    vector<vector<int>> residue_cluster_ind;
    vector<vector<int>> residue_chain_ind;
    int crosslink_in_cluster = 0;

    ifstream *input_cluster;
    ofstream *output_cluster;

    ANALYSIS_PATCH_NO_ORDER(PSF *system, GROUP *sel1, GROUP *sel2, vector<GROUP*> sels, int vector1d, int vector2d, int voidf, string input_cluster_name, string output_cluster_name, string filename, float dist_crit, float dr,float cellsizex, float cellsizey, float cellsizez, int crosslink_in_cluster, int pre_crosslinking); //constructor
    
    void init();

    vector<vector<vector<int>>> head_cell(vector<vector<int>> segments_ind);
    int neighbor_cell_ind(int i, int i_incr, int n);
    string patchtype(string name1, string resname1, string name2, string resname2  );
    void flagallinres(int segid, int resid);
    void flagallinresiffunctionalized(int segid, int resid);
    void flagfunctionalizationiffunctionalized(int segid, int resid);
    void flagfunctionalization(int segid, int resid);
    void read_cluster_file();
    bool is_empty(std::ifstream *pFile);
    void initialize_clusters();
    void initialize_chains();
    void find_initial_clusters();
    void find_initial_chains();
    int form_ring_styr(int ind1,int ind2);
    int form_ring_dvb(int ind1,int ind2);
    void organize_clusters(int cluster_id, int cluster_main_id);
    void reduce_clusters();
    void reduce_clusters_corr();
    void reduce_chains_corr();
    void merge_clusters(int cluster1, int cluster2);
    void merge_chains(int cluster1, int cluster2);
    void select_atoms(GROUP *atoms_select);

    void compute_void();
     
    ~ANALYSIS_PATCH_NO_ORDER();

};

#endif	/* DCD_R_HPP */

