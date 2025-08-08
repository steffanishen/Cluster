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


#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <algorithm>

#include <math.h>
#include "analysis_patch_no_order.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_PATCH_NO_ORDER::ANALYSIS_PATCH_NO_ORDER(PSF *system, GROUP *sel1, GROUP *sel2, vector<GROUP*> sels, int vector1d, int vector2d, int voidf, string input_cluster_name, string output_cluster_name, string filename, float dist_crit, float dr, float cellsizex, float cellsizey, float cellsizez, int crosslink_in_cluster, int pre_crosslinking)
{
    this->system = system;
    this->sels = sels;

//    this->sel1 = sel1;
//    this->sel2 = sel2;

    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->dist_crit = dist_crit;
    this->dr = dr;
    nbins = int(this->dist_crit/dr);
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
    this->system->crosslinking_flag.resize(system->NATOM,0);
    this->system->functionalizing_flag.resize(system->NATOM,0);
    this->cellsizex = cellsizex;
    this->cellsizey = cellsizey;
    this->cellsizez = cellsizez;
    //char fileSpec[filename.length()+1];
    //snprintf(fileSpec, sizeof(fileSpec),"%s",filename.c_str());
    this->file_temp = new ofstream(filename);
    this->seg_res_ids = new ofstream("seg_res_ids.dat");
   // fstream *input_cluster = new fstream(input_cluster_name);
    this->crosslink_in_cluster = crosslink_in_cluster;
    this->pre_crosslinking = pre_crosslinking;

    input_cluster = new ifstream(input_cluster_name);
    linkedlist.clear(); //Meng: clear the linked list
/*
    try
   {
      input_cluster.open(input_cluster_name.c_str());
      if (input_cluster.fail()) throw input_cluster_name;  // the exception being checked
      /*
      while (input_cluster.good())  // check next character
      {
         cout << descrip << ' ' << price << endl;
         inFile >> descrip >> price;
      }
      inFile.close();
      return 0;
   }
   catch(string e)
   {
      cout << e << " was not successfully opened.\n Please check that the file currently exists." << endl;
      exit(1);
   }
*/

/*
    input_cluster->exceptions(std::ifstream::failbit);
    try
    {
        input_cluster->open(input_cluster_name,ios::in);
    }
    catch(ifstream::failure e)
    {
        cerr << "Exception opening/reading file '" << input_cluster_name << "' : " << std::endl;
        cerr << "Please check the path of the file and if it exists." << endl;
    }

*/
    
    output_cluster = new ofstream(output_cluster_name);
   // clusters.resize(system->segments.size(), vector<CLUSTER*>(system->segments[0].size(), new CLUSTER(-1)));
    residue_cluster_ind.resize(system->segments.size(), vector<int>(system->segments[0].size() ));
    residue_chain_ind.resize(system->segments.size(), vector<int>(system->segments[0].size(),-1 ));

    system->segid_ind.clear();
    system->resid_ind.clear();


// make sure all segid and resid start from 0 but still keep the order from the psf file.

    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            for (int i2 = 0; i2 < system->segments[i][i1].size(); i2++) {
                int atom_id_temp = system->segments[i][i1][i2];
                system->segid_ind.push_back(i);
                system->resid_ind.push_back(i1);
            }
        }
    }


}

void ANALYSIS_PATCH_NO_ORDER::init() {
}

vector<vector<vector<int>>> ANALYSIS_PATCH_NO_ORDER::head_cell(vector<vector<int>> segments_ind) {
    vector<vector<vector<int>>> head;
    head.resize(xcount,vector<vector<int>>(ycount,vector<int>(zcount,-1)));//,vector<int>(zcount)));

    int ixcell,iycell,izcell;
     for (auto &segment:segments_ind) {
	    for (int ind : segment) {
            ixcell = int ((system->x[ind] + system->pbc[0]*0.5)/cellsize);
            if (ixcell < 0) ixcell = 0;
            else if (ixcell > xcount-1) ixcell = xcount - 1;

            iycell = int ((system->y[ind] + system->pbc[2]*0.5)/cellsize);
            if (iycell < 0) iycell = 0;
            else if (iycell > ycount-1) iycell = ycount - 1;

            izcell = int ((system->z[ind] + system->pbc[5]*0.5)/cellsize);
            if (izcell < 0) izcell = 0;
            else if (izcell > zcount-1) izcell = zcount - 1;

            linkedlist[ind] = head[ixcell][iycell][izcell];
            head[ixcell][iycell][izcell] = ind;
        }
    }
    return head;
}


//Find the neighbor cell index considering the periodic boundary condition
int ANALYSIS_PATCH_NO_ORDER::neighbor_cell_ind(int i, int i_incr, int n) {
    int i_adj;
    i_adj = i + i_incr;
    if (i_adj < 0) i_adj += n;
    if (i_adj >= n) i_adj -= n;
    return i_adj;
}


string ANALYSIS_PATCH_NO_ORDER::patchtype(string name1, string resname1, string name2, string resname2) {
    string patchtype;
    if (name1 == "C2" && resname1 == "STYR" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST1";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST2";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST4";
    else if (name1 == "C2" && resname1 == "DVB" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST2";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C1" && resname2 == "DVB" ) patchtype = "PST2";
    else if (name1 == "C1" && resname1 == "DVB" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C2" && resname2 == "DVB" ) patchtype = "PST3";
    else if (name1 == "C3" && resname1 == "DVB" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PSTVB7";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C3" && resname2 == "DVB" ) patchtype = "PSTVB8";
    else if (name1 == "C4" && resname1 == "DVB" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PSTVB3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C4" && resname2 == "DVB" ) patchtype = "PSTVB4";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO31";
    else if (name1 == "CE1"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO11";
    else if (name1 == "CD1"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO21";
    else if (name1 == "CD2"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO41";
    else if (name1 == "CE2"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO51";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO32";
    else if (name1 == "CE1"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO12";
    else if (name1 == "CD1"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO22";
    else if (name1 == "CD2"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO42";
    else if (name1 == "CE2"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO52";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO33";
    else if (name1 == "CE1"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO13";
    else if (name1 == "CD1"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO23";
    else if (name1 == "CD2"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO43";
    else if (name1 == "CE2"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO53";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO34";
    else if (name1 == "CE1"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO14";
    else if (name1 == "CD1"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO24";
    else if (name1 == "CD2"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO44";
    else if (name1 == "CE2"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO54";
    else error1.error_exit("Cannot find the patch type!!"); 
    return patchtype;
}

// Flag all atoms of a residue for crosslinking
void ANALYSIS_PATCH_NO_ORDER::flagallinres(int segid_ind, int resid_ind) {
    int natoms = system->segments[segid_ind][resid_ind].size();
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid_ind][resid_ind][i];
        system->crosslinking_flag[flag_index] = 1;
    }
}

// Flag all atoms of a residue of a functional group for functionalization 
void ANALYSIS_PATCH_NO_ORDER::flagallinresiffunctionalized(int segid_ind, int resid_ind) {
    int natoms = system->segments[segid_ind][resid_ind].size();
    int flag_temp = 0;
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid_ind][resid_ind][i];
        if (system->atomname[flag_index] == "C1" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C2" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C3" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C4" && system->charge[flag_index] > -0.15) flag_temp = 1;
    }

    if (flag_temp == 1) {
        for (int i = 0; i < natoms; i++) {
            int flag_index = system->segments[segid_ind][resid_ind][i];
            system->functionalizing_flag[flag_index] = 1;
        }
    }
}

// Flag all atoms of a residue of a chain residue (STYR) or crosslinker (DVB) for functionalization 
 
void ANALYSIS_PATCH_NO_ORDER::flagfunctionalizationiffunctionalized(int segid_ind, int resid_ind) {
    int natoms = system->segments[segid_ind][resid_ind].size();
    int flag_temp = 0;
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid_ind][resid_ind][i];
        if ((system->atomname[flag_index] == "CD1" || system->atomname[flag_index] == "CD2" || system->atomname[flag_index] == "CE1" || system->atomname[flag_index] == "CE2") && (system->resname[flag_index] == "STYR" || system->resname[flag_index] == "DVB") && system->charge[flag_index] > -0.05) flag_temp = 1;
        else if (system->atomname[flag_index] == "CG" && system->resname[flag_index] == "STYR" && system->charge[flag_index] > -0.05 ) flag_temp = 1;
    }

    if (flag_temp == 1) {
        for (int i = 0; i < natoms; i++) {
            int flag_index = system->segments[segid_ind][resid_ind][i];
            if (((system->atomname[flag_index] == "CD1" || system->atomname[flag_index] == "CD2" || system->atomname[flag_index] == "CE1" || system->atomname[flag_index] == "CE2") && (system->resname[flag_index] == "STYR" || system->resname[flag_index] == "DVB")) || (system->atomname[flag_index] == "CG" && system->resname[flag_index] == "STYR") ) {
                system->functionalizing_flag[flag_index] = 1;
            }
        }
    }
}

void ANALYSIS_PATCH_NO_ORDER::flagfunctionalization(int segid_ind, int resid_ind) {
    int natoms = system->segments[segid_ind][resid_ind].size();

    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid_ind][resid_ind][i];
        system->functionalizing_flag[flag_index] = 1;
    }
}


bool ANALYSIS_PATCH_NO_ORDER::is_empty(std::ifstream *pFile)
{
    return pFile->peek() == std::ifstream::traits_type::eof();
}


//Initialize each residue as a cluster
void ANALYSIS_PATCH_NO_ORDER::initialize_clusters() {
    int icluster = 0;
    nresidues = 0;
    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            residue_cluster_ind[i][i1] = icluster;
            CLUSTER *cluster_temp = new CLUSTER(icluster);
            cluster_temp->residue_members.push_back({i,i1});
            clusters.push_back(cluster_temp);


            icluster++;
            nresidues++;
        }
    }
    nclusters = icluster;

}

//Initialize each residue of STYR as a chain, and identify their DVB head and tail caps

void ANALYSIS_PATCH_NO_ORDER::initialize_chains() {
    int ichain = 0;
    nresidues = 0;
    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            int ind1 = system->segments[i][i1][0];
            if (system->resname[ind1] == "STYR") {
                residue_chain_ind[i][i1] = ichain;
                CLUSTER *chain_temp = new CLUSTER(ichain);
                chain_temp->residue_members.push_back({i,i1});
                
// Stores the segid and resid of the DVB crosslinker in the chain. But can a chain be connected to multiple DVB?
// A: A "chain" is either a linear STYR section between 2 DVBs or a dangling STYR chain attached 1 DVB, or an isolated STYR chain
// Note: Without DVB, STYR can only be linear since each STYR monomer has two polymerization sites only.

                for (int i2 = 0; i2 < system->segments[i][i1].size(); i2++) {
                    int atom_index = system->segments[i][i1][i2];
                    if (system->atomname[atom_index] == "C1") {
                        for (auto &ibonded: system->ibond_atom[atom_index]) {
                            if (system->resname[ibonded] == "DVB") {
                                chain_temp->head_cap[0] = system->segid_ind[ibonded];
                                chain_temp->head_cap[1] = system->resid_ind[ibonded];
                                chain_temp->head_C1 = atom_index;
                            }
                        }

                    } else if (system->atomname[atom_index] == "C2") {
                        for (auto &ibonded: system->ibond_atom[atom_index]) {
                            if (system->resname[ibonded] == "DVB") {
                                chain_temp->tail_cap[0] = system->segid_ind[ibonded];
                                chain_temp->tail_cap[1] = system->resid_ind[ibonded];
                                chain_temp->tail_C2 = atom_index;
                            }
                        }

                    } 
                        
                }

                chains.push_back(chain_temp);

                ichain++;

            }

            nresidues++;
        }
    }
    nchains = ichain;

}



// Add the segid and resid (cluster_id) of the current cluster to its ancestor cluster (cluster_main_id)
// // Add the segid and resid (cluster_id) of the kid cluster of the current cluster to its ancestor cluster (cluster_main_id)
// // Deprecated, seems it's not being used given some issues, for example, the members are not informed of their updated clusters
void ANALYSIS_PATCH_NO_ORDER::organize_clusters(int cluster_id, int cluster_main_id) {
    if (clusters[cluster_id]->kid_clusters.size() > 0) {
        clusters[cluster_main_id]->residue_members.insert(clusters[cluster_main_id]->residue_members.end(),clusters[cluster_id]->residue_members.begin(),clusters[cluster_id]->residue_members.end());
        for (int i = 0; i < clusters[cluster_id]->kid_clusters.size(); i++) {
            organize_clusters(clusters[cluster_id]->kid_clusters[i]->cluster_ind, cluster_main_id);
        }
    }
}

void ANALYSIS_PATCH_NO_ORDER::reduce_clusters() {
    int icluster = 0;
    int segid_temp;
    int resid_temp;
    for (int i = 0; i < clusters.size(); i++) {
        if (clusters[i]->parent_cluster_ind == -1) {
            clusters[i]->cluster_ind = icluster;
            for (int i1 = 0; i1 < clusters[i]->residue_members.size(); i1++) {
                segid_temp = clusters[i]->residue_members[i1][0];
                resid_temp = clusters[i]->residue_members[i1][1];
                residue_cluster_ind[segid_temp][resid_temp] = icluster;
            }
            clusters[i]->kid_clusters.clear();
            icluster++;
        }
    }

    auto it = clusters.begin();
    while (it != clusters.end()) {
        if ( (*it)->parent_cluster_ind != -1) {
//            cout << "parent_cluster_ind: " << (*it)->parent_cluster_ind << endl;
            it = clusters.erase(it);
        } else {
            ++it;
        }
    }
    nclusters = clusters.size();

}


// Remove empty clusters
void ANALYSIS_PATCH_NO_ORDER::reduce_clusters_corr() {
    int icluster = 0;
    int segid_temp;
    int resid_temp;
    for (int i = 0; i < clusters.size(); i++) {
        if (clusters[i]->cluster_ind != -1) {
            clusters[i]->cluster_ind = icluster;
            for (int i1 = 0; i1 < clusters[i]->residue_members.size(); i1++) {
                segid_temp = clusters[i]->residue_members[i1][0];
                resid_temp = clusters[i]->residue_members[i1][1];
                residue_cluster_ind[segid_temp][resid_temp] = icluster;
            }
            icluster++;
        }
    }

    auto it = clusters.begin();
    while (it != clusters.end()) {
        if ( (*it)->cluster_ind == -1) {
//            cout << "parent_cluster_ind: " << (*it)->parent_cluster_ind << endl;
            it = clusters.erase(it);
        } else {
            ++it;
        }
    }
    nclusters = clusters.size();

}

//deprecated.
void ANALYSIS_PATCH_NO_ORDER::reduce_chains_corr() {
    int ichain = 0;
    int segid_temp;
    int resid_temp;
    for (int i = 0; i < chains.size(); i++) {
        if (chains[i]->cluster_ind != -1) {
            chains[i]->cluster_ind = ichain;
            for (int i1 = 0; i1 < chains[i]->residue_members.size(); i1++) {
                segid_temp = chains[i]->residue_members[i1][0];
                resid_temp = chains[i]->residue_members[i1][1];
                residue_chain_ind[segid_temp][resid_temp] = ichain;
            }
            ichain++;
        }
    }

    auto it = chains.begin();
    while (it != chains.end()) {
        if ( (*it)->cluster_ind == -1) {
//            cout << "parent_cluster_ind: " << (*it)->parent_cluster_ind << endl;
            it = chains.erase(it);
        } else {
            ++it;
        }
    }
    nchains = chains.size();
        
}





// merge clusters, add members of the cluster of larger id to the cluster of the smaller id
void ANALYSIS_PATCH_NO_ORDER::merge_clusters(int cluster1, int cluster2) {
    int minid = min(cluster1,cluster2);
    int maxid = max(cluster1,cluster2);
    //clusters[minid]->residue_members.insert(clusters[minid]->residue_members.end(),clusters[maxid]->residue_members.begin(),clusters[maxid]->residue_members.end());
    for (int i = 0; i < clusters[maxid]->residue_members.size(); i++) {
        int segid_temp = clusters[maxid]->residue_members[i][0];
        int resid_temp = clusters[maxid]->residue_members[i][1];
        residue_cluster_ind[segid_temp][resid_temp] = minid;
        clusters[minid]->residue_members.push_back({segid_temp,resid_temp});
    }
    
    clusters[maxid]->residue_members.clear();
    clusters[maxid]->cluster_ind = -1;
}

// merge the chains
// update the head cap of chain 1 and tail cap of chain 2
// Clear the chain of larger chain id
void ANALYSIS_PATCH_NO_ORDER::merge_chains(int chain1, int chain2) {
    int minid = min(chain1,chain2);
    int maxid = max(chain1,chain2);
    //clusters[minid]->residue_members.insert(clusters[minid]->residue_members.end(),clusters[maxid]->residue_members.begin(),clusters[maxid]->residue_members.end());
 // Meng: Seems to be opposite, debug!!
 //The following debugs!
 /*
    chains[chain1]->head_C1 = chains[chain2]->head_C1;
    chains[chain1]->head_cap[0] = chains[chain2]->head_cap[0];
    chains[chain1]->head_cap[1] = chains[chain2]->head_cap[1]; 
    chains[chain2]->tail_C2 = chains[chain1]->tail_C2;
    chains[chain2]->tail_cap[0] = chains[chain1]->tail_cap[0];
    chains[chain2]->tail_cap[1] = chains[chain1]->tail_cap[1];
*/
    chains[chain2]->head_C1 = chains[chain1]->head_C1;
    chains[chain2]->head_cap[0] = chains[chain1]->head_cap[0];
    chains[chain2]->head_cap[1] = chains[chain1]->head_cap[1]; 
    chains[chain1]->tail_C2 = chains[chain2]->tail_C2;
    chains[chain1]->tail_cap[0] = chains[chain2]->tail_cap[0];
    chains[chain1]->tail_cap[1] = chains[chain2]->tail_cap[1];



    for (int i = 0; i < chains[maxid]->residue_members.size(); i++) {
        int segid_temp = chains[maxid]->residue_members[i][0];
        int resid_temp = chains[maxid]->residue_members[i][1];
        residue_chain_ind[segid_temp][resid_temp] = minid;
        chains[minid]->residue_members.push_back({segid_temp,resid_temp});
    }

    
    chains[maxid]->residue_members.clear();
    chains[maxid]->cluster_ind = -1;
    chains[maxid]->head_cap.resize(2,-1);
    chains[maxid]->tail_cap.resize(2,-1);
    chains[maxid]->head_C1 = -1;
    chains[maxid]->tail_C2 = -1;
}

//Find the clusters of initial structure before crosslinking
void ANALYSIS_PATCH_NO_ORDER::find_initial_clusters() {
    for (int i = 0; i < system->NBONDS; i++) {
        int iatom = system->ibond[i][0];
            int jatom = system->ibond[i][1];
            int segid_temp0 = system->segid_ind[iatom];
            int resid_temp0 = system->resid_ind[iatom];
            int segid_temp1 = system->segid_ind[jatom];
            int resid_temp1 = system->resid_ind[jatom];
            int cluster1;
            int cluster2;
            if (segid_temp0 != segid_temp1 || resid_temp0 != resid_temp1) {
                cluster1 = residue_cluster_ind[segid_temp0][resid_temp0];
                cluster2 = residue_cluster_ind[segid_temp1][resid_temp1]; 
                if (cluster1 != cluster2) {
                    merge_clusters(cluster1,cluster2);
                }
            }
    }

    reduce_clusters_corr();

}

//Find the chains of initial structure before crosslinking
void ANALYSIS_PATCH_NO_ORDER::find_initial_chains() {
    for (int i = 0; i < system->NBONDS; i++) {
        int iatom = system->ibond[i][0];
        if (system->resname[iatom] == "STYR") {
            int jatom = system->ibond[i][1];
            if (system->resname[jatom] == "STYR") {
                int segid_temp0 = system->segid_ind[iatom];
                int resid_temp0 = system->resid_ind[iatom];
                int segid_temp1 = system->segid_ind[jatom];
                int resid_temp1 = system->resid_ind[jatom];
                int chain1;
                int chain2;
                if (segid_temp0 != segid_temp1 || resid_temp0 != resid_temp1) {
                    chain1 = residue_chain_ind[segid_temp0][resid_temp0];
                    chain2 = residue_chain_ind[segid_temp1][resid_temp1]; 
                    if (chain1 != chain2) {
                        if (system->atomname[iatom] == "C1" && system->atomname[jatom] == "C2") {
                            merge_chains(chain1,chain2);
                        } else if (system->atomname[iatom] == "C2" && system->atomname[jatom] == "C1") {
                            merge_chains(chain2,chain1);
                        }
                    }
                }
            }
        }
    }

    reduce_chains_corr();

}


//See if the crosslinking forms a ring
int ANALYSIS_PATCH_NO_ORDER::form_ring_styr(int ind1,int ind2) {
    int form_ring = 0;
    int segid1_ind = system->segid_ind[ind1];
    int segid2_ind = system->segid_ind[ind2];
    int resid1_ind = system->resid_ind[ind1];
    int resid2_ind = system->resid_ind[ind2];

    int chain1 = residue_chain_ind[segid1_ind][resid1_ind];
    int chain2 = residue_chain_ind[segid2_ind][resid2_ind];

    if (chain1 == chain2) form_ring = 1;

    if (system->atomname[ind1] == "C1" && system->atomname[ind2] == "C2") {
        if (!(chains[chain1]->tail_cap[0] == -1 && chains[chain2]->head_cap[0] == -1)) {
            if (chains[chain1]->tail_cap[0] == chains[chain2]->head_cap[0] && chains[chain1]->tail_cap[1] == chains[chain2]->head_cap[1]) {
                form_ring = 1;
            }
        }
        
    } else if (system->atomname[ind1] == "C2" && system->atomname[ind2] == "C1") {
        if (!(chains[chain1]->head_cap[0] == -1 && chains[chain2]->tail_cap[0] == -1)) {
            if (chains[chain1]->head_cap[0] == chains[chain2]->tail_cap[0] && chains[chain1]->head_cap[1] == chains[chain2]->tail_cap[1]) {
                form_ring = 1;
            }
        }
    } else {
        cout << "STYR1: " << system->atomname[ind1] << "; STYR2: " << system->atomname[ind2] << endl;
        error1.error_exit("The atomnames for crosslinking type is wrong!!");
    }

    return form_ring;
}


int ANALYSIS_PATCH_NO_ORDER::form_ring_dvb(int ind1,int ind2) {
    int form_ring = 0;
    int segid1_ind = system->segid_ind[ind1];
    int segid2_ind = system->segid_ind[ind2];
    int resid1_ind = system->resid_ind[ind1];
    int resid2_ind = system->resid_ind[ind2];

    int ichain = residue_chain_ind[segid1_ind][resid1_ind];

    if (system->atomname[ind1] == "C1") {
        if (chains[ichain]->tail_cap[0] == segid2_ind && chains[ichain]->tail_cap[1] == resid2_ind) {
            form_ring = 1;
        }
        
    } else if (system->atomname[ind1] == "C2") {
        if (chains[ichain]->head_cap[0] == segid2_ind && chains[ichain]->head_cap[1] == resid2_ind) {
            form_ring = 1;
        }
    } else {
        error1.error_exit("The atomnames for crosslinking type is wrong!!");
    }

    return form_ring;
}



void ANALYSIS_PATCH_NO_ORDER:: select_atoms(GROUP *atoms_select) {
//    vector<int> flagarray(NATOM,0);
    int NATOM; // Number of atoms
    vector<ATOM> atoms;

    string temp;
    vector<string> arg1;
    int segid_run = -1;
    int resid_run = -1;
    vector<int> residues_item;
    vector<vector<int>> segments_item;  

    for (auto &atom : atoms_select->atoms) {
        int atom_index = atom.atom_index;
        int flag = 0;
        for (auto &ibonded: system->ibond_atom[atom_index]) {
            int segid1_ind = system->segid_ind[atom_index];        
            int resid1_ind = system->resid_ind[atom_index];       
            if (! (system->segid[ibonded] == system->segid[atom_index] && system->resid[ibonded] == system->resid[atom_index])) {
                system->crosslinking_flag[atom_index] = 1;
                flag = 1;
                if (this->pre_crosslinking == 1 && system->resname[atom_index] == "STYR") {
                    flagallinres(segid1_ind,resid1_ind);
                }
            }
        }
        if (flag == 0 && system->functionalizing_flag[atom_index] == 0) atoms.push_back(atom);
    }



    atoms_select->atoms.clear();
    copy(atoms.begin(),atoms.end(),back_inserter(atoms_select->atoms));

    atoms_select->residues.clear();
    atoms_select->segments.clear();
    atoms_select->segments_ind.clear();
    atoms_select->atomids.clear();

    atoms_select->NATOM = atoms_select->atoms.size();

    for (size_t i = 0; i < atoms_select->NATOM; i++ ) {
        ATOM atom = atoms_select->atoms[i];
	if (i == 0) {
            residues_item.push_back(atom.atom_index);
	        resid_run = system->resid[atom.atom_index];
	        segid_run = system->segid[atom.atom_index];
            if (i == atoms_select->NATOM-1) {
	            atoms_select->residues.push_back(residues_item);
	            segments_item.push_back(residues_item);
	            atoms_select->segments.push_back(segments_item);                
            }
        } else if (i == size_t(atoms_select->NATOM -1)) {
  	    if (system->resid[atom.atom_index] != resid_run || system->segid[atom.atom_index] != segid_run) {
	        atoms_select->residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = system->resid[atom.atom_index];
		residues_item.clear();
	    }
	    if (system->segid[atom.atom_index] != segid_run) {
	        atoms_select->segments.push_back(segments_item);
		segid_run = system->segid[atom.atom_index];
		segments_item.clear();
	    }
	    residues_item.push_back(atom.atom_index);
	    atoms_select->residues.push_back(residues_item);
	    segments_item.push_back(residues_item);
	    atoms_select->segments.push_back(segments_item);
	    residues_item.clear();
	    segments_item.clear();
        } else {
  	    if (system->resid[atom.atom_index] != resid_run || system->segid[atom.atom_index] != segid_run) {
	        atoms_select->residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = system->resid[atom.atom_index];
		residues_item.clear();
	    }
	    if (system->segid[atom.atom_index] != segid_run) {
	        atoms_select->segments.push_back(segments_item);
		segid_run = system->segid[atom.atom_index];
		segments_item.clear();
	    }
            residues_item.push_back(atom.atom_index);
        }
//        cout << "debug: we monitor atom " << atom.atom_index <<" the size of segments_item " <<  segments_item.size() << endl; 
    }

    for (auto &segment : atoms_select->segments) {
	vector<int> segment_temp;
	for (auto &residue : segment) {
	    for (auto &atomid : residue) {
	        segment_temp.push_back(atomid);
	        atoms_select->atomids.push_back(atomid); //MSe comment: group vector with collapsed atomid
	    }
	}
	atoms_select->segments_ind.push_back(segment_temp);
    }


}






void ANALYSIS_PATCH_NO_ORDER::compute_void() {
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    vector<float> distances;
    vector<vector<int>> candidate_pairs;
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    float dr = this->dr;
    int nbins = int(this->dist_crit/dr);
    this->iframe += 1;
    vector<vector<vector<vector<int>>>> heads;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;

    

    vector<float> rdf(nbins,0.0);
    int nsels = sels.size();
    if ( (nsels % 2) != 0) error1.error_exit("ERROR: Odd number of groups. Please select pairs of groups for patches!!");
    int npairs = nsels/2;

// Consider if residues are NC4 functionalized or not
    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            int ind1 = system->segments[i][i1][0];
            if (system->resname[ind1] == "NC4") flagallinresiffunctionalized(i, i1);
            if (system->resname[ind1] == "STYR" || system->resname[ind1] == "DVB"  ) flagfunctionalizationiffunctionalized(i,i1);
        }
    }

// Remove the crosslinked or functionalized atoms
    for (auto &crosslinking_sel: sels) {
        select_atoms(crosslinking_sel);
    }


//Check if the updated selection groups are empty
    vector<int> is_not_empty_group;
    for (int iselpair = 0; iselpair < nsels/2; iselpair++) {
        int is_not_empty_group_temp = 1;
        int isel1 = iselpair*2;
        int isel2 = isel1 + 1;
        if (sels[isel1]->NATOM == 0) {
            is_not_empty_group_temp = 0;
            cout << "Select group " << isel1 << " is empty!!" << endl;
        }
        if (sels[isel2]->NATOM == 0) {
            is_not_empty_group_temp = 0;
            cout << "Select group " << isel2 << " is empty!!" << endl;
        }
        is_not_empty_group.push_back(is_not_empty_group_temp);

        //cout << "segments_ind size: " << sels[i]->segments_ind.size()<< endl;
    }
    
    int sum_is_not_empty_group = 0;
    for (int i = 0; i < nsels/2; i++) {
        sum_is_not_empty_group += is_not_empty_group[i];
    }

    if (sum_is_not_empty_group == 0) error1.error_exit("ERROR: sels don't contain any atoms!");



    if (system->pbc[0] < 0.01 || system->pbc[2] < 0.01 || system->pbc[5] < 0.01 ) error1.error_exit("ERROR: Box size not specified!");
//    cout << "M_PI" << M_PI << endl;// for debug purpose

    xcount = int(system->pbc[0]/cellsizex);
    ycount = int(system->pbc[2]/cellsizey);
    zcount = int(system->pbc[5]/cellsizez);

    linkedlist.resize(system->NATOM,-1);

//    vector<vector<vector<int>>> head1 = head_cell(sel1->segments_ind);
//    vector<vector<vector<int>>> head2 = head_cell(sel2->segments_ind);

    for (int i = 0; i < nsels; i++) {
        vector<vector<vector<int>>> headtemp = head_cell(sels[i]->segments_ind);
        heads.push_back(headtemp);
    }

    int xcount = heads[0].size();
    int ycount = heads[0][0].size();
    int zcount = heads[0][0][0].size();


    for (unsigned i = 0; i < sels.size(); i++) {
        for (auto &segment:sels[i]->segments_ind) {
	        for (int ind : segment) {
                if (system->functionalizing_flag[ind] != 1 && system->crosslinking_flag[ind] != 1) {
                    cout << "Atomname: " << system->atomname[ind] << "; Resname: " << system->resname[ind] << "; index: " << system->atom_index[ind] << " ; segid: " << system->segid[ind] << "; resid: " << system->resid[ind] << " ; Charge: " << system->charge[ind] << endl;
                }
            }
        }
    }



// Initialize the residue clusters
    initialize_clusters();

// Initialize the residue chains
    initialize_chains();

int max_n_chains = chains.size();
int max_n_clusters = clusters.size();

// Merge the initial clusters and merge the initial chains
    find_initial_clusters(); 
    find_initial_chains(); 


    cout << "nresidues: " << nresidues << "; nclusters: " << nclusters << endl;

// Find the distance between possible croslinking pairs
    for (int iselpair = 0; iselpair < nsels/2; iselpair++) {
        int isel1 = iselpair*2;
        int isel2 = isel1 + 1;
        for (int i = 0; i < xcount; i++) {
            for (int j = 0; j < ycount; j++) {
                for (int k = 0; k < zcount; k++) {
                    int ind1 = heads[isel1][i][j][k];
                    while (1) {
                        if (ind1 < 0) break;
                            r[0] = system->x[ind1];
	                    r[1] = system->y[ind1];
	                    r[2] = system->z[ind1];
                        for (int iprime = max(-1,-xcount+1); iprime <= min(1,xcount-1); iprime++) {
                            int i2 = neighbor_cell_ind(i,iprime,xcount);
                            for (int jprime = max(-1,-ycount+1); jprime <= min(1,ycount-1); jprime++) {
                                int j2 = neighbor_cell_ind(j,jprime,ycount);
                                for (int kprime = max(-1,-zcount+1); kprime <= min(1,zcount-1); kprime++) {
                                    int k2 = neighbor_cell_ind(k,kprime,zcount);
                                    int ind2 = heads[isel2][i2][j2][k2];
                                    while (1) {
                                        if (ind2 < 0) break;
                                        if (!(system->segid[ind1] == system->segid[ind2] && system->resid[ind1] == system->resid[ind2])) {


                                            r1[0] = system->x[ind2];
    	                                    r1[1] = system->y[ind2];
    	                                    r1[2] = system->z[ind2];

                                            disp = getDistPoints(r, r1);
                                            dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
                                            dist = sqrt(dist2);
                                    
                                            distances.push_back(dist);
                                            candidate_pairs.push_back({ind1,ind2});

                                        }
                                        ind2 = linkedlist[ind2];
                                    }
                                }
                            }
                        }
                        ind1 = linkedlist[ind1];
                    }
                }
            }
        }
    }

    //vector<float> vector_test = distances;

    vector<int> a_ind = heapSort(distances,distances.size());




    for (int i = 0; i < a_ind.size(); i++ ) {
        int i_sorted = a_ind[i];
        float local_dist = distances[i];
        int ind1 = candidate_pairs[i_sorted][0];
        int ind2 = candidate_pairs[i_sorted][1];
        string segid1 = system->segname[ind1];        
        string segid2 = system->segname[ind2]; 
        int resid1 = system->resid[ind1];       
        int resid2 = system->resid[ind2];       


        int segid1_ind = system->segid_ind[ind1];        
        int segid2_ind = system->segid_ind[ind2]; 
        int resid1_ind = system->resid_ind[ind1];       
        int resid2_ind = system->resid_ind[ind2];       

        int segid_temp;
        int resid_temp;
        
        int functionalizing = 0;
        int functionalized = 0;

        int cluster1;
        int cluster2;

        int chain1 = residue_chain_ind[segid1_ind][resid1_ind];
        int chain2 = residue_chain_ind[segid2_ind][resid2_ind];

        int form_ring = 0;
        int ichain;


// begin debugging
//            if (system->resid[ind1] == 649 || system->resid[ind2] == 649) cout << "Debugging print resid 649: " << system->resid[ind1] << "; ind2: " << ind2 << "; system->crosslinking_flag[ind2]: " << system->crosslinking_flag[ind2] << "; system->resid[ind2]: " << system->resid[ind2] << endl;
//end debugging


        if (local_dist < dist_crit && system->crosslinking_flag[ind1] == 0 && system->crosslinking_flag[ind2] == 0) {
        /*    if (system->resname[ind1] == "NC4") {
                functionalizing = 1;
                segid_temp = system->segid[ind1];
                resid_temp = system->resid[ind1];
                flagallinresifcrosslinked(segid_temp, resid_temp);
            }
            if (system->resname[ind2] == "NC4") {
                functionalizing = 1;
                segid_temp = system->segid[ind2];
                resid_temp = system->resid[ind2];
                flagallinresifcrosslinked(segid_temp, resid_temp);
            }
        */



            if (system->resname[ind1] == "NC4" || system->resname[ind2] == "NC4") {
                functionalizing = 1;
                if (system->functionalizing_flag[ind1] == 1 || system->functionalizing_flag[ind2] == 1) functionalized = 1;
            }


            //if not (functionalizing and functionzlized) 
                if (! (functionalizing == 1 && functionalized == 1) ) {
            // if the two residues are not in the same cluster
                    cluster1 = residue_cluster_ind[segid1_ind][resid1_ind];
                    cluster2 = residue_cluster_ind[segid2_ind][resid2_ind]; 
                    if (cluster1 != cluster2) {
                        cout << "merge cluster1: " << cluster1 << " and cluster2: " << cluster2 << endl;
                        merge_clusters(cluster1,cluster2);
                        if (system->resname[ind1] == "STYR" && system->resname[ind2] == "STYR") {
                            chain1 = residue_chain_ind[segid1_ind][resid1_ind];
                            chain2 = residue_chain_ind[segid2_ind][resid2_ind];

                            if (system->atomname[ind1] == "C1" && system->atomname[ind2] == "C2") {
                                merge_chains(chain1,chain2);
                            } else if (system->atomname[ind1] == "C2" && system->atomname[ind2] == "C1") {
                                merge_chains(chain2,chain1);
                            }
                        } else if (system->resname[ind1] == "STYR" && system->resname[ind2] == "DVB") {
                            ichain = residue_chain_ind[segid1_ind][resid1_ind];
                            if (system->atomname[ind1] == "C1") {
                                chains[ichain]->head_cap[0] = segid2_ind;
                                chains[ichain]->head_cap[1] = resid2_ind;
                            } else if (system->atomname[ind1] == "C2") { 
                                chains[ichain]->tail_cap[0] = segid2_ind;
                                chains[ichain]->tail_cap[1] = resid2_ind;
                            }
                        } else if (system->resname[ind2] == "STYR" && system->resname[ind1] == "DVB") {
                            ichain = residue_chain_ind[segid2_ind][resid2_ind];
                            if (system->atomname[ind2] == "C1") {
                                chains[ichain]->head_cap[0] = segid1_ind;
                                chains[ichain]->head_cap[1] = resid1_ind;
                            } else if (system->atomname[ind2] == "C2") { 
                                chains[ichain]->tail_cap[0] = segid1_ind;
                                chains[ichain]->tail_cap[1] = resid1_ind;
                            }
                        } 

                        string patchtype = this->patchtype(system->atomname[ind1],system->resname[ind1],system->atomname[ind2],system->resname[ind2]);
                        *this->file_temp << "patch " << patchtype << " " << segid1 << ":" << resid1 << " " << segid2 << ":" << resid2 << endl;
                        cout << "patch " << patchtype << " " << segid1 << ":" << resid1 << " " << segid2 << ":" << resid2 << "; distance: " << local_dist << endl;
                        //*this->seg_res_ids << segid1 << " " << resid1 << " " << segid2 << " " << resid2 << endl;
                        *this->seg_res_ids << ind1 << " " << ind2 << endl;
                        system->crosslinking_flag[ind1] = 1;
                        system->crosslinking_flag[ind2] = 1;

                        if (system->resname[ind1] == "NC4") {
                            flagallinres(segid1_ind, resid1_ind);
                            flagfunctionalization(segid2_ind, resid2_ind);
                        } else if (system->resname[ind2] == "NC4") {
                            flagallinres(segid2_ind, resid2_ind);
                            flagfunctionalization(segid1_ind, resid1_ind);
                        }

                        if (this->pre_crosslinking == 1) {
                            if (system->resname[ind1] == "STYR") {
                                flagallinres(segid1_ind, resid1_ind);
                            } else if (system->resname[ind2] == "STYR") {
                                flagallinres(segid2_ind, resid2_ind);
                            }
                        }
                        
                    } else if (crosslink_in_cluster == 1) {

                        if (system->resname[ind1] == "STYR" && system->resname[ind2] == "STYR") {
                            chain1 = residue_chain_ind[segid1_ind][resid1_ind];
                            chain2 = residue_chain_ind[segid2_ind][resid2_ind];

                            form_ring = form_ring_styr(ind1,ind2);
                            if (form_ring == 1) continue;

                            if (system->atomname[ind1] == "C1" && system->atomname[ind2] == "C2") {
                                    merge_chains(chain1,chain2);
                                } else if (system->atomname[ind1] == "C2" && system->atomname[ind2] == "C1") {
                                    merge_chains(chain2,chain1);
                                }
                            } else if (system->resname[ind1] == "STYR" && system->resname[ind2] == "DVB") {

                                ichain = residue_chain_ind[segid1_ind][resid1_ind];

                                form_ring = form_ring_dvb(ind1,ind2);
                                if (form_ring == 1) continue;

                                if (system->atomname[ind1] == "C1") {
                                    chains[ichain]->head_cap[0] = segid2_ind;
                                    chains[ichain]->head_cap[1] = resid2_ind;
                                } else if (system->atomname[ind1] == "C2") { 
                                    chains[ichain]->tail_cap[0] = segid2_ind;
                                    chains[ichain]->tail_cap[1] = resid2_ind;
                                }
                            } else if (system->resname[ind2] == "STYR" && system->resname[ind1] == "DVB") {

                                ichain = residue_chain_ind[segid2_ind][resid2_ind];

                                form_ring = form_ring_dvb(ind2,ind1);
                                if (form_ring == 1) continue;

                                if (system->atomname[ind2] == "C1") {
                                    chains[ichain]->head_cap[0] = segid1_ind;
                                    chains[ichain]->head_cap[1] = resid1_ind;
                                } else if (system->atomname[ind2] == "C2") { 
                                    chains[ichain]->tail_cap[0] = segid1_ind;
                                    chains[ichain]->tail_cap[1] = resid1_ind;
                                }
                            }

                        string patchtype = this->patchtype(system->atomname[ind1],system->resname[ind1],system->atomname[ind2],system->resname[ind2]);
                        *this->file_temp << "patch " << patchtype << " " << segid1 << ":" << resid1 << " " << segid2 << ":" << resid2 << endl;
                        cout << "patch " << patchtype << " " << segid1 << ":" << resid1 << " " << segid2 << ":" << resid2 << "; distance: " << local_dist << endl;
                        //*this->seg_res_ids << segid1 << " " << resid1 << " " << segid2 << " " << resid2 << endl;
                        *this->seg_res_ids << ind1 << " " << ind2 << endl;
                        system->crosslinking_flag[ind1] = 1;
                        system->crosslinking_flag[ind2] = 1;


                        if (system->resname[ind1] == "NC4") {
                            flagallinres(segid1_ind, resid1_ind);
                            flagfunctionalization(segid2_ind, resid2_ind);
                        } else if (system->resname[ind2] == "NC4") {
                            flagallinres(segid2_ind, resid2_ind);
                            flagfunctionalization(segid1_ind, resid1_ind);
                        }

                        if (this->pre_crosslinking == 1) {
                            if (system->resname[ind1] == "STYR") {
                                flagallinres(segid1_ind, resid1_ind);
                            } else if (system->resname[ind2] == "STYR") {
                                flagallinres(segid2_ind, resid2_ind);
                            }
                        }


                    }
                }
            //}

        }
    }

    cout << "N_pairs: " << a_ind.size() << endl; 

vector<int> member_in_chains(max_n_chains + 1,0);
vector<int> member_in_clusters(max_n_clusters + 1,0);

for (auto &chain_temp: chains) {
    int n_member_in_chains = chain_temp->residue_members.size();
    member_in_chains[n_member_in_chains] += 1;
}

for (auto &cluster_temp: clusters) {
    int n_member_in_clusters = cluster_temp->residue_members.size();
    member_in_clusters[n_member_in_clusters] += 1;
}

ofstream *f_chain_distribution = new ofstream("chain_distribution.dat");
ofstream *f_cluster_distribution = new ofstream("cluster_distribution.dat");

for (int i = 0; i < max_n_chains + 1; i++ ) {
    *f_chain_distribution << member_in_chains[i] << endl;
}

for (int i = 0; i < max_n_clusters + 1; i++ ) {
    *f_cluster_distribution << member_in_clusters[i] << endl;
}


}




ANALYSIS_PATCH_NO_ORDER::~ANALYSIS_PATCH_NO_ORDER()
{
    system = NULL;
    //sel1 = NULL;
    //sel2 = NULL;
    this->rdf_count.clear();
//    fclose(this->outfile_box);
    this->seg_res_ids->close();
    this->file_temp->close();
}

