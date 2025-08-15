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

//***************** Contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************
//***************** use of this class must follow anglezs_rings **************************

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <utility>  // for std::pair
#include <algorithm>
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include "analysis.hpp"
#ifndef ANALYSIS_PERSISTENCE_HPP
#define	ANALYSIS_PERSISTENCE_HPP

using namespace std;

class  ANALYSIS_PERSISTENCE : public ANALYSIS
{

private:
    //private attributes
    //and private methods
    long long make_edge_key(int a, int b) {
        return ((long long)std::min(a,b) << 32) | (unsigned int)std::max(a,b);
    }

    void mark_edge(int a, int b) {
        edgesVisited[make_edge_key(a,b)] = true;
    }

    bool edge_visited(int a, int b) {
        return edgesVisited[make_edge_key(a,b)];
    }


    //std::vector<int> traverse_chain(int start, int next, vector<vector<int>> adjacency) {
    std::vector<int> traverse_chain(int start, int next) {
        // Identify a line section that either start from a hub/end and end at a hub/end or
        std::vector<int> chain = {start, next};
        int prev = start;
        int curr = next;
        while (adj_list[curr].size() == 2) {
            int next_node = (adj_list[curr][0] == prev) ? adj_list[curr][1] : adj_list[curr][0];
            if (edge_visited(curr, next_node)) break;
            mark_edge(curr, next_node);
            chain.push_back(next_node);
            prev = curr;
            curr = next_node;
        }
        return chain;
    }
    
public:
   
    // no public attributes
    // public methods
//    GROUP *sel1;
 //   int whichN;
    //
    vector<vector<int>> adj_list;
    ofstream *xyz_file;
    ofstream *networks_xyz_file;
    ofstream *chains_xyz_file;

    std::deque<int> hubs, ends, mids;
    std::unordered_map<long long, bool> edgesVisited;
    //std::unordered_set<std::pair<int, int>, hash<std::pair<int, int>>> visited_edges; // track edges
    //std::unordered_set<std::pair<int,int>, pair_hash> visited_edges;

    int iframe;
    string xyz_filename;
    vector<int> prot_id;

    ANALYSIS_PERSISTENCE(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, string xyz_filename, string networks_xyz_filename, string chains_xyz_filename,  float dist_crit, float cellsize); //constructor
    
    //vector<vector<int>> adjacency_list();
    void adjacency_list();

    void init();

    vector<float> compute_vector();
    void wrap_positions();
     
    ~ANALYSIS_PERSISTENCE();

};

#endif	/* DCD_R_HPP */

