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
#include "analysis_cluster.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_CLUSTER::ANALYSIS_CLUSTER(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, string xyz_filename, float dist_crit, float cellsize)
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->xyz_filename = xyz_filename;
    xyz_file = new ofstream (xyz_filename.c_str());
    this->dist_crit = dist_crit;
    this->cellsize = cellsize;
    //this->rdf_count.resize(nbins);
    //fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    cout << "Natoms: " << sel1->NATOM << endl;
    this->iframe = 0;
    this->prot_id.resize(sel1->NATOM);
    fill(this->prot_id.begin(),this->prot_id.end(),0);


}

void ANALYSIS_CLUSTER::init() {
}

void ANALYSIS_CLUSTER::wrap_positions() {

    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
	        system->x[ind] = remainder(system->x[ind],system->pbc[0]);
	        system->y[ind] = remainder(system->y[ind],system->pbc[2]);
	        system->z[ind] = remainder(system->z[ind],system->pbc[5]);

	    }

    }
}

vector<vector<int>> ANALYSIS_CLUSTER::adjacency_list() {

    vector<int> linkedlist(sel1->NATOM,-1);
    vector<vector<int>> adj_list(sel1->NATOM);
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
                                    adj_list[i].push_back(j);
                                    adj_list[j].push_back(i);
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





    return adj_list;
}

vector<float> ANALYSIS_CLUSTER::compute_vector() {
    sel1->anglezs.clear();
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
    vector<vector<int>> adj_list = adjacency_list();

    //debug
    //cout << "adjacency list: " << endl;
    //for (vector<int> adj_l : adj_list) {
    //    cout << adj_l.size() << endl;
    //}
    //end debug

    // Third, do breadth first search (BFS) to identify clusters
        // Initialize the remaining proteins
    std::deque<int> remaining_proteins(sel1->NATOM);
    std::iota(remaining_proteins.begin(), remaining_proteins.end(), 0); // Fill with 0, 1, ..., natoms-1
        // Starting from current cluster_id, identify the clusters. Keep isolated proteins' cluster_id=0
    int current_cluster_id = 0;
                //cout << "cluster_id: " << current_cluster_id << " remaining proteins: " << remaining_proteins.size() << endl;//debug

    while (!remaining_proteins.empty()) { // Equivalent to `while remaining_proteins:`
                //cout << "cluster_id: " << current_cluster_id << endl;//debug
        int i = remaining_proteins.front(); // Get the first element
        remaining_proteins.pop_front();    // Remove the first element (popleft)

        std::vector<int> current_cluster; // Keep track of atoms assigned to the current cluster
        current_cluster.push_back(i); // Add atom 'i' to the current cluster
        
        std::deque<int> queue; // BFS queue
        queue.push_back(i);

        // Check if `i` has neighbors.
        // In the original Python, `len(adj_list[i]) > 0` is used. 
        // In C++, the BFS will effectively only run if neighbors exist.
        // So this explicit check is not strictly necessary but can be kept for clarity.
        if (!adj_list[i].empty()) { 
            current_cluster_id++;
            cluster_id[i] = current_cluster_id;


            while (!queue.empty()) {
                int current_node = queue.front();
                queue.pop_front();

                for (int neighbor : adj_list[current_node]) { // Iterate through neighbors
                    // Check if neighbor not in current_cluster
                    // Python's `in` operator on lists can be slow.
                    // If `current_cluster` becomes very large, consider using `std::set` or `std::unordered_set`
                    // for faster lookup. However, for a single cluster's atoms, `std::find` on a `std::vector` is often acceptable.
                    if (std::find(current_cluster.begin(), current_cluster.end(), neighbor) == current_cluster.end()) {
                        current_cluster.push_back(neighbor);
                        queue.push_back(neighbor);
                        cluster_id[neighbor] = current_cluster_id;

                        // This is the most problematic part for direct translation: `remaining_proteins.remove(neighbor)`
                        // `std::deque::erase` can be relatively slow if the element is not at the ends.
                        // A more efficient approach for the `remaining_proteins` set is often a boolean `visited` array, 
                        // or to handle removal differently. For now, a direct translation is shown.
                        remaining_proteins.erase(std::remove(remaining_proteins.begin(), remaining_proteins.end(), neighbor), remaining_proteins.end());
                    }
                }
            }
        }
        
        clusters.push_back(current_cluster); // Add the completed cluster to the list of all clusters
    }

    for (vector<int> cluster : clusters) {
        cluster_size.push_back(cluster.size());
    }






    *xyz_file << sel1->NATOM << endl;
    *xyz_file << endl;

    int i = 0;
    for (auto &segment:sel1->segments_ind) {
	    for (int ind : segment) {
            *xyz_file << clusters[cluster_id[i]].size() << " " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; 
            //*xyz_file << cluster_id[i] << " " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << " " << system->pbc[1] << endl; 
            i++;
        }
    }


    return cluster_size;
}


ANALYSIS_CLUSTER::~ANALYSIS_CLUSTER()
{
    system = NULL;
    sel1 = NULL;
    this->prot_id.clear();
    this->rdf_count.clear();
}

