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
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include "analysis.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS::ANALYSIS()
{
    this->system = NULL;
    this->sel1 = NULL;
    this->sel2 = NULL;
    this->whichN = 0;
}

//get distance between atoms
vector<float> ANALYSIS::getDist(ATOM a1, ATOM a2) {
    vector<float> dist;    
    vector<double> box;
    box.resize(3);
    box[0] = system->pbc[0];
    box[1] = system->pbc[2];
    box[2] = system->pbc[5];
    dist.resize(3);
    for (int i = 0; i <= 2; i++) {
        dist[i] = a1.coord[i] - a2.coord[i];
        dist[i] = dist[i] - box[i] * round( dist[i]/box[i] );
    }
    return dist;
}

//get distance between vectors
vector<float> ANALYSIS::getDistPoints(vector<float> p1, vector<float> p2) {
    vector<float> dist;    
    vector<double> box;
    box.resize(3);
    box[0] = system->pbc[0];
    box[1] = system->pbc[2];
    box[2] = system->pbc[5];
    dist.resize(3);
    for (int i = 0; i <= 2; i++) {
        dist[i] = p1[i] - p2[i];
        dist[i] = dist[i] - box[i] * round( dist[i]/box[i] );
    }
    return dist;
}

vector<vector<float>> ANALYSIS::matrix_mult(vector<vector<float>> m1,vector<vector<float>> m2) {
    int nrows = m1.size();
    int ncols = m2[0].size();
    vector<vector<float>> prod(nrows,vector<float>(ncols));

    for (int i=0; i < m1.size();i++) {
        for (int j=0; j < m2[i].size(); j++) {
	    prod[i][j] = 0.0;
        }
    }

    for (int i=0; i < m1.size();i++) {
        for (int j=0; j < m2[i].size(); j++) {
	    for (int k=0; k < m2.size(); k++) {
	        prod[i][j] = prod[i][j] + m1[i][k] * m2[k][j];
	        //cout << "prod: " << prod[i][j] << endl;
	    }
	}
    }
  
 
    return prod;
}

float ANALYSIS::dot_product(vector<float> v1,vector<float> v2) {
    float dot;
    dot = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        dot = dot + v1[i]*v2[i];
    }
    return dot;
}

vector<float> ANALYSIS::cross_product(vector<float> v1,vector<float> v2) {
    vector<float> cross(3);
    cross[0] =  v1[1]*v2[2] - v1[2]*v2[1];
    cross[1] = -v1[0]*v2[2] + v1[2]*v2[0];
    cross[2] =  v1[0]*v2[1] - v1[1]*v2[0];
    return cross;
}

vector<float> ANALYSIS::scalar_product(vector<float> v1, float scalar1) {
    vector<float> v2(v1.size());
    for (int i = 0; i < v1.size(); i++) {
        v2[i] = v1[i] * scalar1;
    }
    return v2;
}

vector<float> ANALYSIS::vector_sum(vector<float> v1, vector<float> v2) {
    vector<float> v3(v1.size());
    for (int i = 0; i < v1.size(); i++) {
	v3[i] = v1[i] + v2[i];
    }
    return v3;
}

vector<double> ANALYSIS::vector_sum(vector<double> v1, vector<double> v2) {
    vector<double> v3(v1.size());
    for (int i = 0; i < v1.size(); i++) {
	v3[i] = v1[i] + v2[i];
    }
    return v3;
}


vector<float> ANALYSIS::vector_subtract(vector<float> v1, vector<float> v2) {
    vector<float> v3(v1.size());
    for (int i = 0; i < v1.size(); i++) {
	v3[i] = v1[i] - v2[i];
    }
    return v3;
}

float ANALYSIS::norm(vector<float> v1) {
    float norm1 = dot_product(v1,v1);
    norm1 = sqrt(norm1);
    return norm1;
}

vector<float> ANALYSIS::vector_quarternion_prod(vector<float> v1, QUARTERNION q) { // v1 should be a 3d vector

    if (v1.size() != 3) {
        cout << "vector size is wrong!! Should be 3D for quarternion operations!" << endl;
	exit(1);
    }

    vector<float> t1(v1.size());
    vector<float> t2(v1.size());
    vector<float> t3(v1.size());
    vector<float> t4(v1.size());
    vector<float> r2(v1.size());

    t1 = scalar_product(v1,(q.q0)*(q.q0));
    r2 = t1;

    float f2 = dot_product(q.qv,q.qv);
    t2 = scalar_product(v1,-f2);
    r2 = vector_sum(r2,t2);

    vector<float> vt3 = cross_product(q.qv,v1);
    t3 = scalar_product(vt3,2*q.q0);
    r2 = vector_sum(r2,t3);

    float f4 = dot_product(q.qv,v1);
    t4 = scalar_product(q.qv,2*f4);
    r2 = vector_sum(r2,t4);

    return r2;
}

vector<float> ANALYSIS::rotate(vector<float> v1,float theta, vector<float> s) {
    QUARTERNION q;
    q.q0 = cos(0.5*theta);
    q.qv = scalar_product(s,sin(0.5*theta));
    vector<float> r2;
    r2 = vector_quarternion_prod(v1,q);
    return r2;
}


vector<float> ANALYSIS::rotate(vector<float> v1,float sintheta, float costheta, vector<float> s) {
    vector<float> t1(v1.size());
    vector<float> t2(v1.size());
    vector<float> t3(v1.size());
    vector<float> t4(v1.size());
    vector<float> r2(v1.size());

    float f1 = (1 + costheta) * 0.5;
    t1 = scalar_product(v1,f1);
    r2 = t1;
    float f2 = -(1-costheta) * 0.5 * dot_product(s,s);
    t2 = scalar_product(v1,f2);
    r2 = vector_sum(r2,t2);
    float f3 = sintheta;
    t3 = scalar_product(cross_product(s,v1),f3);
    r2 = vector_sum(r2,t3);
    float f4 = (1 - costheta) * dot_product(s,v1);
    t4 = scalar_product(s,f4);
    r2 = vector_sum(r2,t4);

    return r2;
    
}


vector<vector<float>> ANALYSIS::bond_next(int residueid, PSF *monomer) {

    vector<vector<float>> bond_vectors;
    return bond_vectors;

}

vector<vector<float>> ANALYSIS::bond_pre(int residueid, PSF *monomer) {

    vector<vector<float>> bond_vectors;
    return bond_vectors;

}

float ANALYSIS::dangle(float a1, float a2) {

    float dA;

    dA = a1 - a2;
    if ( dA > 180.0 ) {
	dA = dA - 360.0;
    } else if (dA < -180.0) {
	dA = dA + 360.0;
    }

    return dA;
}


float ANALYSIS::angle_z(vector<float> dist) {
    float an_z;
    dist.resize(3);
    if (dist[0] > -1.0e-9 && dist[0] < 1.0e-9) {
        if (dist[1] >= 0.0) {
            an_z = 90.0;
        } else {
            an_z = -90.0;
        }
    } else {
        an_z = atan(dist[1]/dist[0]) * 180.0/PI;
        if (dist[0] < 0.0) {
            an_z = an_z + 180.0;
        }
    }

    return an_z;
}


vector<float> ANALYSIS::bond_rotate_vec(int ind1, int ind2, int ind3, vector<float> vec_ref, vector<float> vec_ref1, vector<float> vec_bf_rot) { //rotate vec_bf_rot along two orthogonal axes: 1. the cross product of vec_ref and vec_CH, and 2. vec_CH, with the angle two planes sharing vec_CH. 

    vector<float> vec_rot(vec_ref.size());
    vector<float> vec1 = {system->x[ind1],system->y[ind1],system->z[ind1]};
    vector<float> vec2 = {system->x[ind2],system->y[ind2],system->z[ind2]};
    vector<float> vec3 = {system->x[ind3],system->y[ind3],system->z[ind3]};
    vector<float> vec_CH = vector_subtract(vec2,vec1);
    vector<float> vec_CC_pre = vector_subtract(vec3,vec1);
    vector<float> vec_CC_pre_rot(3);
   	
    vector<float> s1 = cross_product(vec_ref,vec_CH);
    float len_s1 = norm(s1);
    float len_ref = norm(vec_ref);	    
    float len_CH = norm(vec_CH);	    

    float sintheta = len_s1/(len_ref*len_CH);
    float costheta = dot_product(vec_ref,vec_CH)/(len_ref*len_CH);


    if (len_s1 > 1e-6) {
	s1 = scalar_product(s1,1.0/len_s1);   
	vec_rot = rotate(vec_bf_rot,sintheta,costheta,s1);
	vec_CC_pre_rot = rotate(vec_ref1,sintheta,costheta,s1);
    } else {
	vec_rot = scalar_product(vec_bf_rot, costheta);
	vec_CC_pre_rot = scalar_product(vec_ref1, costheta);
    }



    vector<float> n1 = cross_product(vec_CC_pre,vec_CH); // put the axis the first, actually, doesn't matter since the crossproduct is invariant if both of them are reversed
    vector<float> n2 = cross_product(vec_CC_pre_rot,vec_CH);
    vector<float> s2 = cross_product(n2,n1);// the second rotation axis is vec_CH
    float len_s2 = norm(s2);
    float len_n1 = norm(n1);
    float len_n2 = norm(n2);

    sintheta = len_s2/(len_n2*len_n1);
    costheta = dot_product(n2,n1)/(len_n2*len_n1); 
   
    
 
    if (len_s2 > 1e-6) {
	s2 = scalar_product(s2,1.0/len_s2);
        vec_rot = rotate(vec_rot,sintheta,costheta,s2); 
    } else {
        vec_rot = scalar_product(vec_rot,costheta);
    }
   

    return vec_rot;

}

vector<int> ANALYSIS::cis_trans_id(vector<PSF*> monomers, string resname) {
    vector<int> indices;
    indices.clear();
    int i = 0;
    for (PSF *monomer : monomers) {
        if (monomer->resname[1] == resname) {
	    indices.push_back(i);
	}
        i++;
    }
    return indices;
}

vector<float> ANALYSIS::linspace(float rlow,float rhigh,float dr) {
    int npoints = int((rhigh - rlow)/dr);
    vector<float> array1(npoints,0.0);
    for (int i = 0; i < array1.size(); i++) {
        array1[i] = rlow + dr * float(i);
    }
    return array1;
}


void ANALYSIS::residue_rotate_reset(vector<float> vec_CC, int ind_pivot, vector<float> vec_ref, int residueid) { // rotate the atoms in residue[residueid] along the plane normal of the plane containing vec_ref and vec_CC from vec_ref to vec_CC
    
   	
    //cout << "entering rotation reset!!!!!!" << endl; //debug
    //cout << "nsize vec_ref: " << vec_ref.size() << " nsize vec_CC: " << vec_CC.size() << endl; //debug
    vector<float> s1 = cross_product(vec_ref,vec_CC); // similarly, put the direction rotating to on the 2nd
    //cout << "determining the rotation axis!!!!!!" << endl; //debug


    float len_s1 = norm(s1);
    float len_ref = norm(vec_ref);	    
    float len_CC = norm(vec_CC);	    

    float sintheta = len_s1/(len_ref*len_CC);
    float costheta = dot_product(vec_ref,vec_CC)/(len_ref*len_CC);

    vector<float> pos_pivot = {system->x[ind_pivot],system->y[ind_pivot],system->z[ind_pivot]};
    vector<float> dist_rot;

    if (len_s1 > 1e-6) {
        s1 = scalar_product(s1,1.0/len_s1);   
    }

//    cout << "before rotation reset!!!!!!" << endl; //debug

    for (int ind : sel1->residues[residueid]) {
	vector<float> dist = {system->x[ind]-pos_pivot[0],system->y[ind]-pos_pivot[1],system->z[ind]-pos_pivot[2]};
        if (len_s1 > 1e-6) {
	    dist_rot = rotate(dist,sintheta,costheta,s1);
	} else {
	    dist_rot = scalar_product(dist,costheta);
  	}
	system->x[ind] = pos_pivot[0] + dist_rot[0];
	system->y[ind] = pos_pivot[1] + dist_rot[1];
	system->z[ind] = pos_pivot[2] + dist_rot[2];
//	if (residueid == 92 || residueid == 93) {
//	    cout << "after resetting: " << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; //debug
//	}
    }

//    cout << "after rotation reset!!!!!!" << endl; //debug

}

// On the inside of the box, the random walk has no constrictions rather than the rigid angle; close to the walls, the random walk inclines to the inside. We take the algorithm more radical than a Monte Carlo for the constriction.
void ANALYSIS::residue_rotate_random(int ind1, int ind2, int ind_pivot, float theta, int residueid, int ind_a4_pre, int ind_a2) {

    x_temp.clear();
    y_temp.clear();
    z_temp.clear();

    fill(hitwall.begin(), hitwall.end(),0);
    E_res = 0.0;
    int hitwall_indicator = 0;

    vector<float> vec1 = {system->x[ind1],system->y[ind1],system->z[ind1]}; //ind1 = ind_a3_pre
    vector<float> vec2 = {system->x[ind2],system->y[ind2],system->z[ind2]};  //ind2 = ind_a4_pre
    vector<float> vec_CC = vector_subtract(vec2,vec1);
    vector<float> sp = scalar_product(vec_CC,1.0/norm(vec_CC));

    vector<float> pos_pivot(3);
    pos_pivot[0] = system->x[ind_pivot];
    pos_pivot[1] = system->y[ind_pivot];
    pos_pivot[2] = system->z[ind_pivot];

    int ind_temp = 0;
    int ind_temp_a2 = 0;
    r2_a2 = 0.0;

    for (int ind: sel1->residues[residueid - 1]) {
   	if (ind > ind_a4_pre) {
 	//  cout << "testing: " << ind << endl; //debug
	    vector<float> dist = {system->x[ind]-pos_pivot[0],system->y[ind]-pos_pivot[1],system->z[ind]-pos_pivot[2]};
	    vector<float> dist_rot = rotate(dist,sin(theta),cos(theta),sp);
	    x_temp.push_back(pos_pivot[0] + dist_rot[0]);
	    y_temp.push_back(pos_pivot[1] + dist_rot[1]);
	    z_temp.push_back(pos_pivot[2] + dist_rot[2]);

	  //	    cout << system->x[ind] << " " << system->y[ind] << " " << system->z[ind] << endl; //debug
	    ind_temp++;
	 }
    }


    for (int ind : sel1->residues[residueid]) {
        vector<float> dist = {system->x[ind]-pos_pivot[0],system->y[ind]-pos_pivot[1],system->z[ind]-pos_pivot[2]};
	vector<float> dist_rot = rotate(dist,sin(theta),cos(theta),sp);
	x_temp.push_back(pos_pivot[0] + dist_rot[0]);
	y_temp.push_back(pos_pivot[1] + dist_rot[1]);
	z_temp.push_back(pos_pivot[2] + dist_rot[2]);


        if (ind == ind_a2) {
	    ind_temp_a2 = ind_temp;
//	    cout << "1st seg: " << endl; //debug
	}
	ind_temp++;
    }

    if (x_temp[ind_temp_a2] > system->pbc[0] * 0.5 - fringe) {
	hitwall[0] = 1;
   	float r = x_temp[ind_temp_a2] - system->pbc[0]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (x_temp[ind_temp_a2] < -system->pbc[0] * 0.5 + fringe) {
	hitwall[1] = 1;
	float r = -x_temp[ind_temp_a2] - system->pbc[0]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (y_temp[ind_temp_a2] > system->pbc[2] * 0.5 - fringe) {
	hitwall[2] = 1;
	float r = y_temp[ind_temp_a2] - system->pbc[2]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (y_temp[ind_temp_a2] < -system->pbc[2] * 0.5 + fringe) {
	hitwall[3] = 1;
	float r = -y_temp[ind_temp_a2] - system->pbc[2]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (z_temp[ind_temp_a2] > system->pbc[5] * 0.5 - fringe) {
	hitwall[4] = 1;
	float r = z_temp[ind_temp_a2] - system->pbc[5]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (z_temp[ind_temp_a2] < -system->pbc[5] * 0.5 + fringe) {
	hitwall[5] = 1;
	float r = -z_temp[ind_temp_a2] - system->pbc[5]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 

    if (hitwall_indicator > 3) {cout << "Hitwall messed!!!!!!!!!!!!!!1" << endl;}
    if (hitwall_indicator > 0) {
	r2_a2 = x_temp[ind_temp_a2] * x_temp[ind_temp_a2] + y_temp[ind_temp_a2] * y_temp[ind_temp_a2] + z_temp[ind_temp_a2] * z_temp[ind_temp_a2];
    }
    //cout << "Hitwall for " << hitwall_indicator << " times." << endl;

}

void ANALYSIS::update_res(int residueid, int ind_a4_pre) {
    int ind_temp = 0;
    for (int ind: sel1->residues[residueid - 1]) {
   	if (ind > ind_a4_pre) {
            system->x[ind] = x_temp[ind_temp];
            system->y[ind] = y_temp[ind_temp];
            system->z[ind] = z_temp[ind_temp];
	    ind_temp++;
	}
    }

    for (int ind : sel1->residues[residueid]) {
        system->x[ind] = x_temp[ind_temp];
        system->y[ind] = y_temp[ind_temp];
        system->z[ind] = z_temp[ind_temp];
	ind_temp++;
    }

    x_temp.clear();
    y_temp.clear();
    z_temp.clear();

}



void ANALYSIS::residue_rotate_random_regular(int ind1, int ind2, int ind_pivot, float theta, int residueid, int ind_a2) {

// Note: 1. ind1, ind2 determin the rotation axis; 2. ind_pivot determines the rotation pivot; 3. ind_a2 determines whether the residue will hit the wall
    x_temp.clear();
    y_temp.clear();
    z_temp.clear();

    fill(hitwall.begin(), hitwall.end(),0);
    E_res = 0.0;

    vector<float> vec1 = {system->x[ind1],system->y[ind1],system->z[ind1]}; //ind1 = ind_a3_pre
    vector<float> vec2 = {system->x[ind2],system->y[ind2],system->z[ind2]};  //ind2 = ind_a4_pre
    vector<float> vec_CC = vector_subtract(vec2,vec1);
    vector<float> sp = scalar_product(vec_CC,1.0/norm(vec_CC));

    vector<float> pos_pivot(3);
    pos_pivot[0] = system->x[ind_pivot];
    pos_pivot[1] = system->y[ind_pivot];
    pos_pivot[2] = system->z[ind_pivot];

    int ind_temp = 0;
    int ind_temp_a2 = 0;
    r2_a2 = 0.0;

    int hitwall_indicator = 0;

    for (int ind : sel1->residues[residueid]) {
        vector<float> dist = {system->x[ind]-pos_pivot[0],system->y[ind]-pos_pivot[1],system->z[ind]-pos_pivot[2]};
	vector<float> dist_rot = rotate(dist,sin(theta),cos(theta),sp);
	x_temp.push_back(pos_pivot[0] + dist_rot[0]);
	y_temp.push_back(pos_pivot[1] + dist_rot[1]);
	z_temp.push_back(pos_pivot[2] + dist_rot[2]);
//        cout << "pbc lx: " << system->pbc[0] << endl;
//        cout << "pbc ly: " << system->pbc[2] << endl;
//        cout << "pbc lz: " << system->pbc[5] << endl;

        if (ind == ind_a2) {
	    ind_temp_a2 = ind_temp;
//	    cout << "2nd seg: " << endl; //debug
	} 
	ind_temp++;
    }

    if (x_temp[ind_temp_a2] > system->pbc[0] * 0.5 - fringe) {
	hitwall[0] = 1;
   	float r = x_temp[ind_temp_a2] - system->pbc[0]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (x_temp[ind_temp_a2] < -system->pbc[0] * 0.5 + fringe) {
	hitwall[1] = 1;
	float r = -x_temp[ind_temp_a2] - system->pbc[0]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (y_temp[ind_temp_a2] > system->pbc[2] * 0.5 - fringe) {
	hitwall[2] = 1;
	float r = y_temp[ind_temp_a2] - system->pbc[2]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (y_temp[ind_temp_a2] < -system->pbc[2] * 0.5 + fringe) {
	hitwall[3] = 1;
	float r = -y_temp[ind_temp_a2] - system->pbc[2]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (z_temp[ind_temp_a2] > system->pbc[5] * 0.5 - fringe) {
	hitwall[4] = 1;
	float r = z_temp[ind_temp_a2] - system->pbc[5]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 
    if (z_temp[ind_temp_a2] < -system->pbc[5] * 0.5 + fringe) {
	hitwall[5] = 1;
	float r = -z_temp[ind_temp_a2] - system->pbc[5]*0.5;
	E_res += wall_pos(r);
	hitwall_indicator++;
    } 

    if (hitwall_indicator > 3) {cout << "Hitwall messed!!!!!!!!!!!!!!1" << endl;}
 //   cout << "Hitwall for " << hitwall_indicator << " times." << endl;
    if (hitwall_indicator > 0) {
	r2_a2 = x_temp[ind_temp_a2] * x_temp[ind_temp_a2] + y_temp[ind_temp_a2] * y_temp[ind_temp_a2] + z_temp[ind_temp_a2] * z_temp[ind_temp_a2];
    }

}

void ANALYSIS::update_res_regular(int residueid) {
    int ind_temp = 0;
    for (int ind : sel1->residues[residueid]) {
        system->x[ind] = x_temp[ind_temp];
        system->y[ind] = y_temp[ind_temp];
        system->z[ind] = z_temp[ind_temp];
	ind_temp++;
    }

    x_temp.clear();
    y_temp.clear();
    z_temp.clear();

}




float ANALYSIS::wall_pos(float r) {
    float E = 0.0;
//    float r6 = r2 * r2 * r2;
//    float r12 = r6 * r6;
//    float sigma2 = sigma * sigma;
//    float sigma6 = sigma2 * sigma2 * sigma2;
//    float sigma12 = sigma6 * sigma6;
//    E = 4.0 * epsilon * ( r2 / sigma2 );
    E = -4.0 * epsilon * sigma * sigma * sigma * sigma * sigma / (r * r * r * r * r) ;
    return E;
}

vector<float> ANALYSIS::compute_vector() {
    vector<float> placeholder;
    return placeholder;
}

vector<vector<float>> ANALYSIS::compute_2D_vector() {
    vector<vector<float>> placeholder;
    return placeholder;
}

void ANALYSIS::heapify(vector<float>& arr, vector<int>& a_ind, int n, int i) {
        // Find largest among root, left child and right child
        int largest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        float placeholder;
        int placeholder_int;

        
        if (left < n && arr[left] > arr[largest]) largest = left;
        
        if (right < n && arr[right] > arr[largest]) largest = right;

    // Swap and continue heapifying if root is not largest
        if (largest != i) {
            placeholder = arr[i];
            arr[i] = arr[largest]; 
            arr[largest] = placeholder;

            placeholder_int = a_ind[i];
            a_ind[i] = a_ind[largest]; 
            a_ind[largest] = placeholder_int;

            heapify(arr, a_ind, n, largest);
            }
    }


vector<int> ANALYSIS::heapSort(vector<float>& arr, int n)
{
    // Build heap (rearrange array)
    float placeholder;
    int placeholder_int;

    vector<int> a_ind;
    a_ind.resize(n);

    for (int i = 0; i < n; i++) a_ind[i] = i;

    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, a_ind, n, i);
  
    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        placeholder = arr[0];
        arr[0] = arr[i];
        arr[i] = placeholder;

        placeholder_int = a_ind[0];
        a_ind[0] = a_ind[i];
        a_ind[i] = placeholder_int;
        
  
        // call max heapify on the reduced heap
        heapify(arr, a_ind, i, 0);
    }
    return a_ind;
}


vector<float> ANALYSIS::residue_com(vector<int> residue)
{
    vector<float> com(3,0.0);
    float mass = 0.0;
    for (int ind : residue) {
        mass += system->atoms[ind].mass;
        com[0] += system->x[ind] * system->atoms[ind].mass;
        com[1] += system->y[ind] * system->atoms[ind].mass;
        com[2] += system->z[ind] * system->atoms[ind].mass;
    }
    
    for (int i = 0; i < 3; i++) {
        com[i] = com[i]/mass;
    }

    return com;
}

void ANALYSIS::compute_void() {
}


ANALYSIS::~ANALYSIS()
{
    this->system = NULL;
    this->sel1 = NULL;
    this->sel2 = NULL;
    danglezs.clear();
}

