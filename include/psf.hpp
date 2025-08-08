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


#ifndef PSF_HPP_INCLUDED
#define PSF_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include "atom.hpp"
#include "group.hpp"

using namespace std;

class PSF
{
public:
    ifstream psff;  // for opening dcd file, MS comment: Note this dcdf is an fstream object, not a DCD or DCD_R object.
    
    bool psf_first_read; // at first read of the coordinates if there are some frozen atoms the number of coordinates to read is different than for other frames
    
    int  NTITLE; //how many "title lines" in psf 
    char *TITLE; //each "title line"
    
    /*content of ICNTRL : non detailed ones are 0 */
    int nframes_tot;
    
    int NATOM; // Number of atoms
    int NBONDS; // Number of bonds
    int NANGLES; // Number of angles 
    int NDIHEDRALS; // Number of proper dihedrals
    int NIMPROPERS; // Number of improper dihedrals
    
    int *atom_index; // Array storing indexes of moving atoms. Note the index start from 0.
    int *segid;
    string *segname;
    int *resid;
    vector<int> segid_ind;
    vector<int> resid_ind;
    int *tmin_msd;
    string *resname;
    string *atomname;
    string *atomtype;
    float *charge;
    float *mass;
    int *beta;
    vector<ATOM> atoms;
    vector<vector<int>> residues; 
    vector<vector<vector<int>>> segments; 
    int **ibond;
    int **iangle;
    int **idihedral;
    int **improper;
    int *nbonds_atom;
    vector<vector<int>> ibond_atom; // The indices of atoms bonded to each of the atoms
    string name1 = "NULL";
    string name2 = "NULL";
    string name_ref1 = "NULL";
    string name_ref2 = "NULL";
    string monomer_resname = "NULL";
    int ind_a1 = -1;
    int ind_a2 = -1;  //the atom connected to the next monomer
    int ind_ref1 = -1;
    int ind_ref2 = -1;
    string name_pre = "NULL";
    string name_next = "NULL";
    vector<float> Vch = {0.0, 0.0, 0.0};
    vector<float> Vcc = {0.0, 0.0, 0.0};
    vector<float> Vcc_pre = {0.0, 0.0, 0.0};
    vector<float> Vcc_next = {0.0, 0.0, 0.0};
    int monomerflag = 0;

    float *x;
    float *y;
    float *z;

    vector<int> crosslinking_flag;
    vector<int> functionalizing_flag;

    vector<float*> XS;
    vector<float*> YS;
    vector<float*> ZS;

    vector<vector<vector<int>>> head;

    int iframe;

    double pbc[6];  //a matrix of 6 real defining the periodic boundary conditions : only useful if QCRYS is not 0
    double box_first_frame[3];  //a matrix of 6 real defining the periodic boundary conditions : only useful if QCRYS is not 0
    
    //coordinates stored in simple precision (IMPORTANT)
    
    //protected methods
    void alloc();
    void alloc_bonds();
    void alloc_angles();
    void alloc_dihedrals();
    void alloc_impropers();
    void determine_names();
    vector<float> vector_subtract(vector<float> v1, vector<float> v2);

    PSF(string filename);
    
    
    void read_header();
    void read_atoms();
    void read_bonds();
    void read_angles();
    void read_dihedrals();
    void read_impropers();

    void read_psf();


    GROUP select_atoms(vector<string> str);


    ~PSF();
};

#endif // DCD_HPP_INCLUDED
