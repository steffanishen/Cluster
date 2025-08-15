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
#include <deque>
#include <numeric>
#include <algorithm>

#include <math.h>
#include "analysis_persistence.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_PERSISTENCE::ANALYSIS_PERSISTENCE(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, string xyz_filename, string networks_xyz_filename, string chains_xyz_filename, string rings_xyz_filename, string bondlength_filename,  float dist_crit, float cellsize)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->xyz_filename = xyz_filename;
    xyz_file = new ofstream (xyz_filename.c_str());
    networks_xyz_file = new ofstream (networks_xyz_filename.c_str());
    chains_xyz_file = new ofstream (chains_xyz_filename.c_str());
    rings_xyz_file = new ofstream (rings_xyz_filename.c_str());
    bondlength_file = new ofstream (bondlength_filename.c_str());
    this->dist_crit = dist_crit;
    this->cellsize = cellsize;
    //this->rdf_count.resize(nbins);
    //fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    cout << "Natoms: " << sel1->NATOM << endl;
    this->iframe = 0;
    this->prot_id.resize(sel1->NATOM);
    fill(this->prot_id.begin(),this->prot_id.end(),0);


}

void ANALYSIS_PERSISTENCE::init() {
}

void ANALYSIS_PERSISTENCE::wrap_positions() {

    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
	        system->x[ind] = remainder(system->x[ind],system->pbc[0]);
	        system->y[ind] = remainder(system->y[ind],system->pbc[2]);
	        system->z[ind] = remainder(system->z[ind],system->pbc[5]);

	    }

    }
}

//vector<vector<int>> ANALYSIS_PERSISTENCE::adjacency_list() {
void ANALYSIS_PERSISTENCE::adjacency_list() {

    vector<int> linkedlist(sel1->NATOM,-1);
//    vector<vector<int>> adj_list(sel1->NATOM);
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    float dist_crit2 = this->dist_crit*this->dist_crit;

// First, assign the protein id to the representative atom
    int i = 0;
    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {

            prot_id[i] = ind;

    //cout << "sel1->atom_index: " << system->atom_index[ind] << endl; //for debug purpose
    //cout << "sel1->atomtype: " << system->atomtype[ind] << endl; //for debug purpose
            i++;

	    }
    }

// Second, build the linked cell list
    int xcount = int(system->pbc[0]/cellsize);
    int ycount = int(system->pbc[2]/cellsize);
    int zcount = int(system->pbc[5]/cellsize);

    int ncells = xcount * ycount * zcount;

    vector<int> head(ncells,-1);

    i = 0;
    int ixcell,iycell,izcell,icell;

     for (auto &segment:sel1->segments_ind) {
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

            icell = ixcell * ycount * zcount + iycell * zcount + izcell;

            linkedlist[i] = head[icell];
            head[icell] = i;

            i++;
        }
    }

    // debug
    //for (i = 0; i < ncells; ++i) {
    //    cout << "i: " << i << " head: " << head[i] << endl;
    //}
    // end debug


// Third, build the adjacency list using the linked cell list

    i = 0; 

    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
	        r[0] = system->x[ind];
	        r[1] = system->y[ind];
	        r[2] = system->z[ind];

            //cout << "ind: " << ind << endl; //debug

            int cidx = static_cast<int>((system->x[ind] + system->pbc[0]*0.5) / cellsize);
            int cidy = static_cast<int>((system->y[ind] + system->pbc[2]*0.5) / cellsize);
            int cidz = static_cast<int>((system->z[ind] + system->pbc[5]*0.5) / cellsize);
            for (int ix : {(cidx - 1 + xcount) % xcount, cidx, (cidx + 1) % xcount}) {
                for (int iy : {(cidy - 1 + ycount) % ycount, cidy, (cidy + 1) % ycount}) {
                    for (int iz : {(cidz - 1 + zcount) % zcount, cidz, (cidz + 1) % zcount}) {
                        icell = ix * ycount * zcount + iy * zcount + iz;
                        int j = head[icell];
            
                        while (j >= 0) {
                            if (j > i) {
                               // cout << "icell: " << icell << " j: " << j << endl;

                                int ind1 = prot_id[j];
	                            r1[0] = system->x[ind1];
	                            r1[1] = system->y[ind1];
	                            r1[2] = system->z[ind1];

                               // cout << "ind: " << ind << " ind1: " << ind1 << endl; //debug

	                            vector<float> disp = getDistPoints(r, r1);
                                float dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];

                                if (dist2 < dist_crit2) {
                                    this->adj_list[i].push_back(j);
                                    this->adj_list[j].push_back(i);
                                }
                            /*
                            */

                            }
                            

                            j = linkedlist[j];

                        }
                    }
                }
            }
        i++;
        }
    }


    i = 0;
    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
            if (this->adj_list[i].size() > 2) hubs.push_back(i);
            else if (this->adj_list[i].size() == 1) ends.push_back(i);
            else if (this->adj_list[i].size() == 2) mids.push_back(i);
            i++;
        }
    }


//    return adj_list;
}

vector<float> ANALYSIS_PERSISTENCE::compute_vector() {
    this->hubs.clear();
    this->ends.clear();
    this->mids.clear();
    this->edgesVisited.clear();
    this->adj_list.clear();
    this->adj_list.resize(sel1->NATOM);
    vector<vector<int>> networks;
    vector<vector<int>> chains;
    vector<vector<int>> rings;
    vector<float> polymerbondlength;
    vector<int> chainlength;
    vector<float> correlationij;
    vector<int> correlation_count;
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    vector<vector<int>> clusters;
    vector<int> cluster_id(sel1->NATOM, 0); 
    this->iframe += 1;

    vector<float> cluster_size;

    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose

    // Before everything, wrap the atom coordinates over PBC
    wrap_positions();

    // Next, build the adjacency list
        // First, create the linked cell list
            // Assign atoms to cells
            // Obtain headers of each cell, and linked list of each atom
        // Second, build the adjancency list for each atom
    //vector<vector<int>> adj_list = adjacency_list();
    adjacency_list(); // make adjacency list a class property

    //creat lists of hubs, ends and mids based on the size of adjacency lists
    // Use Depth first search (DFS) to identify chains
    // Step 1: Handle hubs
        while (!hubs.empty()) {
            int hub = hubs.front(); hubs.pop_front();
            for (int nb : adj_list[hub]) {
                if (edge_visited(hub, nb)) continue;
                mark_edge(hub, nb);
                std::vector<int> chain = traverse_chain(hub, nb);
        // We remove dimers or trimers from the network type clusters
                if (chain.size() > 4) networks.push_back(chain);
            }
        }

        // Step 2: Handle remaining ends
        while (!ends.empty()) {
            int end = ends.front(); ends.pop_front();
            for (int nb : adj_list[end]) {
                if (edge_visited(end, nb)) continue;
                mark_edge(end, nb);
                std::vector<int> chain = traverse_chain(end, nb);
                chains.push_back(chain);
            }
        }

    // Step 3: Analyez the rings
        while (!mids.empty()) {
            int mid = mids.front(); mids.pop_front();
            for (int nb : adj_list[mid]) {
                if (edge_visited(mid, nb)) continue;
                mark_edge(mid, nb);
                std::vector<int> ring = traverse_chain(mid, nb);
                rings.push_back(ring);
            }
        }




    *xyz_file << sel1->NATOM << " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    //*xyz_file << sel1->NATOM << endl;//<< " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    *xyz_file << endl;

    int nhubs = 0;
    int i = 0;
    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
            if (adj_list[i].size() > 2) {
                *xyz_file << "hubs " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 
            //*xyz_file << cluster_id[i] << " " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << " " << system->pbc[1] << endl; 
                nhubs += 1;
            } else if (adj_list[i].size() == 1) {
                *xyz_file << "ends " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 
            } else if (adj_list[i].size() == 2) {
                *xyz_file << "mids " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 
            } else {
                *xyz_file << "monomers " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 
            }
            i++;
        }
    }

    int n_networks = 0;
    for (auto &chain: networks) {
        n_networks += chain.size();
    }

    *networks_xyz_file << n_networks << " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    //*xyz_file << sel1->NATOM << endl;//<< " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    *networks_xyz_file << endl;

    for (auto &chain: networks) {
        for (int id_chain: chain){
            int ind = prot_id[id_chain];
                *networks_xyz_file << "networks " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 

        }
    }


    int n_chains = 0;
    for (auto &chain: chains) {
        n_chains += chain.size();
    }

    *chains_xyz_file << n_chains << " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    //*xyz_file << sel1->NATOM << endl;//<< " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    *chains_xyz_file << endl;

    for (auto &chain: chains) {
        for (int id_chain: chain){
            int ind = prot_id[id_chain];
                *chains_xyz_file << "chains " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 

        }
    }


    int n_rings = 0;
    for (auto &chain: rings) {
        n_rings += chain.size();
    }

    *rings_xyz_file << n_rings << " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    //*xyz_file << sel1->NATOM << endl;//<< " " << system->pbc[0] << " " << system->pbc[2] << " " << system->pbc[5]  << endl;
    *rings_xyz_file << endl;

    for (auto &chain: rings) {
        for (int id_chain: chain){
            int ind = prot_id[id_chain];
                *rings_xyz_file << "rings " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 

        }
    }

    for (auto &chain: networks) {
        chainlength.push_back(chain.size());
        for (int id_chain = 0; id_chain < chain.size()-1; id_chain++){
            int ind = prot_id[id_chain] - 1;
            int ind1 = ind + 1;
	        r[0] = system->x[ind];
	        r[1] = system->y[ind];
	        r[2] = system->z[ind];
	        r1[0] = system->x[ind1];
	        r1[1] = system->y[ind1];
	        r1[2] = system->z[ind1];
            polymerbondlength.push_back(norm(getDistPoints(r, r1)));
        }
    }

    for (auto &chain: chains) {
        chainlength.push_back(chain.size());
        for (int id_chain = 0; id_chain < chain.size()-1; id_chain++){
            int ind = prot_id[id_chain] - 1;
            int ind1 = ind + 1;
	        r[0] = system->x[ind];
	        r[1] = system->y[ind];
	        r[2] = system->z[ind];
	        r1[0] = system->x[ind1];
	        r1[1] = system->y[ind1];
	        r1[2] = system->z[ind1];
            polymerbondlength.push_back(norm(getDistPoints(r, r1)));
        }
    }

    for (auto &chain: rings) {
        chainlength.push_back(chain.size());
        for (int id_chain = 0; id_chain < chain.size()-1; id_chain++){
            int ind = prot_id[id_chain] - 1;
            int ind1 = ind + 1;
	        r[0] = system->x[ind];
	        r[1] = system->y[ind];
	        r[2] = system->z[ind];
	        r1[0] = system->x[ind1];
	        r1[1] = system->y[ind1];
	        r1[2] = system->z[ind1];
            polymerbondlength.push_back(norm(getDistPoints(r, r1)));
        }
    }

    float average_bondlength = 0.0;

    auto max_chainlength = max_element(chainlength.begin(),chainlength.end());
    if (max_chainlength != chainlength.end()) {
        cout << "max_chainlength: " << *max_chainlength << endl;
        correlationij.resize(*max_chainlength);
        correlation_count.resize(*max_chainlength);
        fill(correlation_count.begin(),correlation_count.end(),0);

        for (auto bl : polymerbondlength) {
            *bondlength_file << bl << " ";
            average_bondlength += bl;
        }
        *bondlength_file << endl;

        average_bondlength /= static_cast<float>(polymerbondlength.size());
        cout << "average_bondlength: " << average_bondlength << endl;

    }


// We can finally analyze the persistence length now!
    if (max_chainlength != chainlength.end()) {
        vector<vector<int>> all_chains;
        //copy(networks.begin(), networks.end(), back_inserter(all_chains));
        all_chains.insert(all_chains.end(),networks.begin(),networks.end());
        all_chains.insert(all_chains.end(),chains.begin(),chains.end());
        all_chains.insert(all_chains.end(),rings.begin(),rings.end());

        for (auto &chain: chains) {
            for (int id_chain = 0; id_chain < chain.size()-2; id_chain++){
                int ind = prot_id[chain[id_chain]];
                int ind1 = prot_id[chain[id_chain + 1]];
                r[0] = system->x[ind];
                r[1] = system->y[ind];
                r[2] = system->z[ind];
                r1[0] = system->x[ind1];
                r1[1] = system->y[ind1];
                r1[2] = system->z[ind1];
                vector<float> disp = getDistPoints(r, r1);

                for (int id_chain1 = id_chain; id_chain1 < chain.size()-1; id_chain1++){
                    ind = prot_id[chain[id_chain1]];
                    ind1 = prot_id[chain[id_chain1 + 1]];
                    r[0] = system->x[ind];
                    r[1] = system->y[ind];
                    r[2] = system->z[ind];
                    r1[0] = system->x[ind1];
                    r1[1] = system->y[ind1];
                    r1[2] = system->z[ind1];
                    vector<float> disp1 = getDistPoints(r, r1);
                    float costhetaij = dot_product(disp,disp1)/(norm(disp)*norm(disp1));
                    correlationij[id_chain1 - id_chain] += costhetaij;
                    correlation_count[id_chain1 - id_chain] += 1;


                }
                
                

            }
        }

        for (auto id_along_chain=0; id_along_chain < *max_chainlength; id_along_chain++) {
            if (correlation_count[id_along_chain] > 0) {
                correlationij[id_along_chain] /= static_cast<float>(correlation_count[id_along_chain]);
            }
        }


    }



    return correlationij;
}


ANALYSIS_PERSISTENCE::~ANALYSIS_PERSISTENCE()
{
    system = NULL;
    sel1 = NULL;
    this->prot_id.clear();
    this->rdf_count.clear();
}

