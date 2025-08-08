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


#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unordered_map>
#include "analysis_straight_chain.hpp"


#define PI 3.14159265

using namespace std;

//ANALYSIS_STRAIGHT_CHAIN::ANALYSIS_STRAIGHT_CHAIN(PSF *system, GROUP *sel1, int vector1d, int vector2d,int voidf, float bl, string name1, string name2, string name3, string name4, string name_ref1, string name_ref2, string filename, vector<PSF*> monomers)
ANALYSIS_STRAIGHT_CHAIN::ANALYSIS_STRAIGHT_CHAIN(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, float bl, string filename, vector<PSF*> monomers, float fringe) //constructor
{
    this->system = system;
    this->sel1 = sel1;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
//    this->name1 = name1;
//    this->name2 = name2;
//    this->name3 = name3;
//    this->name4 = name4;
//    this->name_ref1 = name_ref1;
//    this->name_ref2 = name_ref2;
//    cout << "name_ref1: " << this->name_ref1 << "; name_ref2: " << this->name_ref2 << endl; //debug
    this->bl = bl;
    for (unsigned i = 0; i < monomers.size(); i++) {
        this->monomers.push_back(monomers[i]);
    }
    this->fringe = fringe;
}

void ANALYSIS_STRAIGHT_CHAIN::init() {
}


vector<vector<float>> ANALYSIS_STRAIGHT_CHAIN::bond_next(int residueid, PSF *monomer) {

    vector<vector<float>> bond_vectors;
    int ind_ref_H;
    int ind_a1;
    int ind_a2;
    int ind_a1_next;
    vector<float> Vch(3);
    vector<float> Vcc(3);
    vector<float> Vcc_pre(3);

    bond_vectors.clear();

    for (int ind : monomer->residues[residueid]) {
         if (monomer->atomname[ind] == monomer->name1) {
	     ind_a1 = ind;
	 }

         if (monomer->atomname[ind] == monomer->name2) {
	     ind_a2 = ind;
	 }
      //       cout << "testing C1-4: " << ind_a1 << " " << ind_a4 << endl;
         if (monomer->atomname[ind] == monomer->name_ref1) {
             ind_ref_H = ind; 
	 }
    }

    residueid++;


    for (int ind : monomer->residues[residueid]) {
         if (monomer->atomname[ind] == monomer->name1) {
             ind_a1_next = ind;
	 }
    }


    vector<float> vec0 = {monomer->x[ind_a1],monomer->y[ind_a1],monomer->z[ind_a1]};
    vector<float> vec1 = {monomer->x[ind_a2],monomer->y[ind_a2],monomer->z[ind_a2]};
    vector<float> vec2 = {monomer->x[ind_ref_H],monomer->y[ind_ref_H],monomer->z[ind_ref_H]};
    vector<float> vec3 = {monomer->x[ind_a1_next],monomer->y[ind_a1_next],monomer->z[ind_a1_next]};


    Vch = vector_subtract(vec2,vec1); // the reference C2-H3 bond
    Vcc = vector_subtract(vec3,vec1); // the reference C2-C1 bond
    Vcc_pre = vector_subtract(vec0,vec1); // the reference C1-C2 bond
    
    bond_vectors.push_back(Vch);
    bond_vectors.push_back(Vcc);
    bond_vectors.push_back(Vcc_pre);

    return bond_vectors;

}


vector<vector<float>> ANALYSIS_STRAIGHT_CHAIN::bond_pre(int residueid, PSF *monomer) {

    vector<vector<float>> bond_vectors;
    int ind_ref_H;
    int ind_a4_ref;
    int ind_a1;
    vector<float> Vch(3);
    vector<float> Vcc(3);

    bond_vectors.clear();

    for (int ind : monomer->residues[residueid]) {
         if (monomer->atomname[ind] == name1) {
	     ind_a1= ind;
	 }
      //       cout << "testing C1-4: " << ind_a1 << " " << ind_a4 << endl;
  //       if (monomer->atomname[ind] == name_ref1) {
             ind_ref_H = ind; 
//	 }
    }


    residueid--;

    for (int ind : monomer->residues[residueid]) {
         if (monomer->atomname[ind] == name4) {
             ind_a4_ref = ind;
	 }
    }

    vector<float> vec1 = {monomer->x[ind_a4_ref],monomer->y[ind_a4_ref],monomer->z[ind_a4_ref]};
//    vector<float> vec2 = {monomer->x[ind_ref_H],monomer->y[ind_ref_H],monomer->z[ind_ref_H]};
    vector<float> vec3 = {monomer->x[ind_a1],monomer->y[ind_a1],monomer->z[ind_a1]};


    //Vch = vector_subtract(vec2,vec3); // the reference CH bond
    Vcc = vector_subtract(vec3,vec1); // the reference C2-C1 bond
    
   // bond_vectors.push_back(Vch);
    bond_vectors.push_back(Vcc);

    return bond_vectors;

}



void ANALYSIS_STRAIGHT_CHAIN::compute_void() {
    FILE * fpdb;
    const char *filechar = filename.c_str();
    fpdb = fopen(filechar,"w");
    int ind_a1 = -1;
    int ind_a2 = -1;
    int ind_a3 = -1;
    int ind_a4 = -1;
    int ind_H = -1;
    int ind_a4_ref = -1;
    int ind_a3_pre = -1;
    int ind_a4_pre = -1;
    int ind_a1_pre = -1;
    int ind_H_pre = -1;
    int ind_a1_end = -1;
    int ind_ref_H1 = -1;
    int ind_ref_H2 = -1;
    int ind_ref_a1 = -1;
    int ind_head = -1;
    vector<float> shift(3);
    vector<float> vec1(3);
    vector<float> vec2(3);
    vector<float> vec_CC(3);
    vector<float> pos_pivot(3);
    vector<float> dist_rot(3);
    int moid;
    int ctsel;
    int moid_pre;
  //  vector<float> next_bond(3);

    cout << "Nresidues: " << sel1->residues.size() << endl;


///////////determine the initial orientation of the C-C bond////////////////////////////////////////

    int resid_ref = 1;
    vector<string> residue_names_redundant;

    //cout << "before determining special atoms: " << endl; //debug

    for (PSF *monomer : monomers) {
	monomer->determine_names();
	residue_names_redundant.push_back(monomer->monomer_resname);
    }
  
 
    sort(residue_names_redundant.begin(),residue_names_redundant.end());
    vector<string> residue_names = vector_unique(residue_names_redundant);

    unordered_map<string,vector<int>> ctids;
    for (string residuename : residue_names) {
	vector<int> ctid = cis_trans_id(monomers,residuename);
	ctids[residuename] = ctid;
    }


// begin commented out
// end commented out

/////////////////////////////end determining the C-C bond///////////////////////

    srand(time(NULL));

//    cout << "before loop: " << endl; //debug

    for (int restempid = 0; restempid < sel1->residues.size(); restempid++) {
        
        if (restempid > 0) {
            ind_a3_pre = ind_a3;
            ind_a4_pre = ind_a4;
            ind_H_pre = ind_H;
	    moid_pre = moid;
	}

// determine the reference monomer
        vector<int> monomerids = ctids[system->resname[sel1->residues[restempid][0]]];
        if (monomerids.size() == 1) {
	    ctsel = 0;
	} else {
	    float ctsel_random = (float) rand()/RAND_MAX;
	    ctsel = (int) round(ctsel_random);            
   	}

	moid = monomerids[ctsel];

// determine the reference atoms
        for (int ind : sel1->residues[restempid]) {
	     if (system->atomname[ind] == monomers[moid]->name1) {
                 ind_a1 = ind;
	     }
	    if (monomers[moid]->monomerflag < 3) {
             if (system->atomname[ind] == monomers[moid]->name2) {
	         ind_a2 = ind;
		 for (int ib = 0; ib < system->nbonds_atom[ind]; ib++) {
		     int ind_temp = system->ibond_atom[ind][ib];
                     if (system->atomname[ind_temp] == monomers[moid]->name_pre && ind_temp < ind) {
   	                 ind_a3 = ind_temp;
	             }
		 }
		 if (monomers[moid]->ind_a1 > monomers[moid]->ind_a2) {
		     ind_a3 = monomers[moid]->ind_a1;
		 }
	     }
// begin comment out
             if (system->atomname[ind] == monomers[moid]->name2) {
	         ind_a4 = ind;
	     }
// end comment out
             if (system->atomname[ind] == monomers[moid]->name_ref1) {
	         ind_H = ind;
	     }
	    }
        }


	int atid = monomers[moid]->ind_a1;
	if (monomers[moid]->ind_a1 > monomers[moid]->ind_a2 && restempid != sel1->residues.size()-1) {
	    atid = 0;
	}
	for (int ind: sel1->residues[restempid]) {
	    system->x[ind] = monomers[moid]->x[atid];
	    system->y[ind] = monomers[moid]->y[atid];
	    system->z[ind] = monomers[moid]->z[atid];
	    atid++;
 	}

//	cout << " monomerflag: " << monomers[moid]->monomerflag << " monomername " << monomers[moid]->monomer_resname << endl; 

        if (restempid > 0 ) {
   
//step1: attach the bond connecting this and the next monomer and the relaxed orientation based on the other atoms in the current monomer connected to C4 
 	    if (restempid == 1) {cout << "a4_pre: " << ind_a4_pre << " ind_H_pre: " << ind_H_pre << " a3_pre: " << ind_a3_pre << endl;}	  
 
	    vec_CC = bond_rotate_vec(ind_a4_pre, ind_H_pre, ind_a3_pre, monomers[moid_pre]->Vch, monomers[moid_pre]->Vcc, monomers[moid_pre]->Vcc_next);



            shift[0] = (system->x[ind_a4_pre] + vec_CC[0]) - system->x[ind_a1];
            shift[1] = (system->y[ind_a4_pre] + vec_CC[1]) - system->y[ind_a1];
            shift[2] = (system->z[ind_a4_pre] + vec_CC[2]) - system->z[ind_a1];


            for (int ind : sel1->residues[restempid]) {
 	        system->x[ind] = system->x[ind] + shift[0];
 	        system->y[ind] = system->y[ind] + shift[1];
 	        system->z[ind] = system->z[ind] + shift[2];
            }



///////////////rotate C2-C1 to default orientation/////////////////////////////
//	    cout <<restempid << " before rotation reset!! " << " resname " << system->resname[sel1->residues[restempid][0]]  << endl; //debug

// 	    cout << "Vcc_pre size: " << monomers[moid]->Vcc_pre.size() << " monomerflag: " << monomers[moid]->monomerflag << " monomername " << monomers[moid]->monomer_resname << endl;

	    residue_rotate_reset(vec_CC,ind_a1,monomers[moid]->Vcc_pre,restempid);
//	    cout <<restempid << " after rotation reset!! " << " resname " << system->resname[sel1->residues[restempid][0]]  << endl; //debug
/////////////////////////end default orientation///////////////////////////////

////////random rotation part///////////////////////////////////////////////
//step 2: rotate the next monomer (and the bonds connecting C4 in the current monomer) along the C3-C4 axis

// Note we don't differentiate between single or double bonds here. We will work on this later.
/*
	    float theta_random = 2 * PI * (float) rand()/RAND_MAX;

//            cout << "output theta: " << theta_random << " " << cos(theta_random) << endl;// debug
            fill(hitwall.begin(), hitwall.end(),0);
            residue_rotate_random(ind_a3_pre, ind_a4_pre, ind_a4_pre, theta_random, restempid, ind_a4_pre, ind_a2);


	    if (hitwall[0] == 0 && hitwall[1] == 0 && hitwall[2] == 0 && hitwall[3] == 0 && hitwall[4] ==0 && hitwall[5] == 0) {
   	        update_res(restempid, ind_a4_pre);
  	    } else {

		vector<float> x_temp_pre = x_temp;
		vector<float> y_temp_pre = y_temp;
		vector<float> z_temp_pre = z_temp;
		float E_res_pre = E_res;
  	        float r2_a2_pre = r2_a2;
		for (int mcstep = 0; mcstep < 50; mcstep++) {
	            theta_random = 2 * PI * float(mcstep) * 0.1;
                    residue_rotate_random(ind_a3_pre, ind_a4_pre, ind_a4_pre, theta_random, restempid, ind_a4_pre, ind_a2);
		//    if (E_res < E_res_pre) {
		    if (r2_a2 < r2_a2_pre) {
			r2_a2_pre = r2_a2;
			x_temp_pre = x_temp;
			y_temp_pre = y_temp;
			z_temp_pre = z_temp;
//			cout << "Find smaller r2 pre: " << endl; //debug
		    }
		}
		x_temp = x_temp_pre;
		y_temp = y_temp_pre;
		z_temp = z_temp_pre;

   	        update_res(restempid, ind_a4_pre);
	    }
*/
/////////step 3: rotate the monomer along the C4_pre-C1 bond
/*
	    if ( ind_a1 != ind_a2 ) {
                theta_random = 2 * PI * (float) rand()/RAND_MAX;          
                fill(hitwall.begin(), hitwall.end(),0);
                residue_rotate_random_regular(ind_a4_pre, ind_a1, ind_a1, theta_random, restempid, ind_a2);


	        if (hitwall[0] == 0 && hitwall[1] == 0 && hitwall[2] == 0 && hitwall[3] == 0 && hitwall[4] ==0 && hitwall[5] == 0) {
	    	    update_res_regular(restempid);
		} else {

		    vector<float> x_temp_pre = x_temp;
		    vector<float> y_temp_pre = y_temp;
		    vector<float> z_temp_pre = z_temp;
		    float E_res_pre = E_res;
  	            float r2_a2_pre = r2_a2;
		    for (int mcstep = 0; mcstep < 50; mcstep++) {
	                theta_random = 2 * PI * float(mcstep) * 0.1;
                	residue_rotate_random_regular(ind_a4_pre, ind_a1, ind_a1, theta_random, restempid, ind_a2);
		   //     if (E_res < E_res_pre) {
		        if (r2_a2 < r2_a2_pre) {
		//	    E_res_pre = E_res;
			    r2_a2_pre = r2_a2;
			    x_temp_pre = x_temp;
			    y_temp_pre = y_temp;
			    z_temp_pre = z_temp;
			    //cout << "Find smaller r2 : " << endl; //debug
		        }
		    }
		    x_temp = x_temp_pre;
		    y_temp = y_temp_pre;
		    z_temp = z_temp_pre;

   	            update_res_regular(restempid);
		}
	    }
*/
//	    cout << "after random rotation!!" << endl; //debug
///////////////end rotation part////////////////////
        }


    }


    for (int restempid = 0; restempid < sel1->residues.size(); restempid++) {
        for (int ind : sel1->residues[restempid]) {
            fprintf(fpdb,"%4s%7d  %-4s%-4s%1d%4d%12.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM", ind + 1, system->atomname[ind].c_str(), system->resname[ind].c_str(), 1, system->resid[ind], system->x[ind], system->y[ind], system->z[ind], 0,0);
        }
    }

    fclose(fpdb);    
}


ANALYSIS_STRAIGHT_CHAIN::~ANALYSIS_STRAIGHT_CHAIN()
{
    system = NULL;
    sel1 = NULL;
}

