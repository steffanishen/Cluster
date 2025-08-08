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

#include <iostream>
#include <vector>
#include "psf.hpp"
#include "atom.hpp"
#include "error1.hpp"
#include "quarternion.hpp"
#include <cmath>
#include <algorithm>
#include <functional>
#ifndef ANALYSIS_HPP
#define	ANALYSIS_HPP

using namespace std;

class ANALYSIS 
{

private:
    //private attributes
    //no private methods
    
public:
   
    PSF *system; // the system object including atom/particle information, must be defined before running the simulation and the box information
    GROUP *sel1;
    GROUP *sel2;
    vector<PSF*> monomers;
    ERROR1 error1;
    int dtmax;
    int whichN;
    int vector1d;
    int vector2d;
    int voidf;
    int nbinsangle;
    int nbins;
    int axisid;
    string name1;
    string name2;
    string name3;
    string name4;
    float bl; // the length of the bond between adjacent monomers
    string filename;
    vector<vector<float>> danglezs; // the minimum dangle for every N of each of all the segments
    string name_ref1;
    string name_ref2;
    vector<float> count; 
    int wall_rw = 0;
    float scaling;
    float lx;
    float ly;
    float lz;
    float sigma = 1.0;
    float epsilon = 1.0;
    vector<int> hitwall = {0,0,0,0,0,0}; //determine if the residue gos outside the wall
//    int hitwall = 0;
    vector<float> x_temp;
    vector<float> y_temp;
    vector<float> z_temp;
    vector<float> msd_previous;
    vector<float> rdf_count;
    int msd_xflag = 0;
    int msd_yflag = 0;
    int msd_zflag = 0;
    float low;
    float high; 
    float E_res = 0.0;
    float fringe = 0.0;
    float r2_a2 = 0.0;
    string which_density_profile = "mass";
    float dist_crit;
    float dr;
    float cellsize = 10.0;
    float cellsizex = 10.0;
    float cellsizey = 10.0;
    float cellsizez = 10.0;
    int pre_crosslinking = 0; //If pre_crosslinking = 1, then form precursor chains before crosslinking. In practice, precrosslink dvb and styr to guarantee formation of networks later on
    int xcount; 
    int ycount; 
    int zcount; 

    // no public attributes
    // public methods
    ANALYSIS(); //constructor
    
    virtual void init()=0;
    vector<float> getDist(ATOM a1, ATOM a2);  //calculate the distance between two atoms/particles considering periodic boundary condition
    vector<float> getDistPoints(vector<float> p1, vector<float> p2); //calculate the distance between two points considering periodic boundary condition
    float dot_product(vector<float> v1,vector<float> v2);
    vector<float> cross_product(vector<float> v1,vector<float> v2);
    vector<float> scalar_product(vector<float> v1, float scalar1);
    vector<float> vector_sum(vector<float> v1, vector<float> v2);
    vector<double> vector_sum(vector<double> v1, vector<double> v2);
    vector<float> vector_subtract(vector<float> v1, vector<float> v2);
    vector<float> vector_quarternion_prod(vector<float> v1, QUARTERNION q);
    vector<float> rotate(vector<float> v1, float theta, vector<float> s);
    vector<float> rotate(vector<float> v1, float sintheta, float costheta, vector<float> s);
    vector<vector<float>> matrix_mult(vector<vector<float>> m1,vector<vector<float>> m2); //Matrix multiplication
    float norm(vector<float> v1);
    float angle_z(vector<float> dist); //calculate the angle of a vector projected on the x-y plane with respect to vector (0,1)
    float dangle(float a1, float a2); // calculate the difference between two angles
    vector<float> bond_rotate_vec(int ind1, int ind2, int ind3, vector<float> vec_ref, vector<float> vec_ref1, vector<float> vec_bf_rot); //rotate vec_bf_rot based on the rotation from vec_ref to vec_ind1_ind2. More in detail: rotate vec_bf_rot along two orthogonal axes: 1. the cross product of vec_ref and vec_CH, and 2. vec_CH, with the angle two planes sharing vec_CH. 
    vector<int> cis_trans_id(vector<PSF*> monomers, string resname); // determine the psf ID for cis-trans choice of different monomers
    float wall_pos(float dist);
   
    template <typename mytype1>
    vector<mytype1> vector_unique(vector<mytype1> vec1) {
        vector<mytype1> vec2;
        int i = 0;
        for (mytype1 vec1in : vec1) {
            if (i == 0) {
	        vec2.push_back(vec1in);
    	    } else {
	        if (vec1in != vec1[i-1]) {
	            vec2.push_back(vec1in);
	        }
	    }
	    i++;
        }
        return vec2;
    }

    template <typename T>
    std::vector<T> vector_add(const std::vector<T>& a, const std::vector<T>& b)
    {
        assert(a.size() == b.size());

        std::vector<T> result;
        result.reserve(a.size());

        std::transform(a.begin(), a.end(), b.begin(), 
        std::back_inserter(result), std::plus<T>());
        return result;
    }


    void residue_rotate_reset(vector<float> vec_CC, int ind_pivot, vector<float> vec_ref, int residueid);
    void residue_rotate_random(int ind1, int ind2, int ind_pivot, float theta, int residueid, int ind_a4_pre, int ind_a2);
    void residue_rotate_random_regular(int ind1, int ind2, int ind_pivot, float theta, int residueid, int ind_a2);
    void update_res(int residueid, int ind_a4_pre);
    void update_res_regular(int residueid);
    vector<float> linspace(float rlow,float rhigh,float dr);

    virtual vector<vector<float>> bond_next(int residueid, PSF *monomer);
    virtual vector<vector<float>> bond_pre(int residueid, PSF *monomer);
    virtual vector<vector<float>> compute_2D_vector();   
    virtual vector<float> compute_vector();

    void heapify(vector<float>& arr, vector<int>& a_ind, int n, int i);
    vector<int> heapSort(vector<float>& arr, int n);
    vector<float> residue_com(vector<int> residue);

    virtual void compute_void();
 
    virtual ~ANALYSIS();

};

#endif	/* DCD_R_HPP */

