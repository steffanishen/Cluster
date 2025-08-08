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
#include "analysis_passivation.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_PASSIVATION::ANALYSIS_PASSIVATION(PSF *system, vector<GROUP*> sels, int vector1d, int vector2d, int voidf, string filename)
{
    this->system = system;
    this->sels = sels;

//    this->sel1 = sel1;
//    this->sel2 = sel2;

    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    //this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
    this->file_temp = new ofstream(filename);
    this->seg_res_ids = new ofstream("seg_res_ids.dat");
   // fstream *input_cluster = new fstream(input_cluster_name);


    system->segid_ind.clear();
    system->resid_ind.clear();

    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            for (int i2 = 0; i2 < system->segments[i][i1].size(); i2++) {
                int atom_id_temp = system->segments[i][i1][i2];
                system->segid_ind.push_back(i);
                system->resid_ind.push_back(i1);
            }
        }
    }

    cout << "Initialized passivation!" << endl;


}

void ANALYSIS_PASSIVATION::init() {
}



string ANALYSIS_PASSIVATION::patchtype(string name1, string resname1, string name2, string resname2) {
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

string ANALYSIS_PASSIVATION::patchtype(string name1, string resname1) {
    string patchtype;
    if (name1 == "C1" && ( resname1 == "STYR" || resname1 == "DVB")) patchtype = "PSS1";
    else if (name1 == "C2" && ( resname1 == "STYR" || resname1 == "DVB")) patchtype = "PSS2";
    else if (name1 == "C3" && resname1 == "DVB" ) patchtype = "PSS3";
    else if (name1 == "C4" && resname1 == "DVB" ) patchtype = "PSS4";
    else error1.error_exit("Cannot find the patch type!!"); 
    return patchtype;
}




bool ANALYSIS_PASSIVATION::is_empty(std::ifstream *pFile)
{
    return pFile->peek() == std::ifstream::traits_type::eof();
}




void ANALYSIS_PASSIVATION:: select_atoms(GROUP *atoms_select) {
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
            if (! (system->segid[ibonded] == system->segid[atom_index] && system->resid[ibonded] == system->resid[atom_index])) flag = 1;
        }
        if (flag == 0) atoms.push_back(atom);
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






void ANALYSIS_PASSIVATION::compute_void() {
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    vector<float> distances;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    float dr = this->dr;
    this->iframe += 1;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;
    

    int nsels = sels.size();


// Remove the crosslinked or functionalized atoms
    for (auto &crosslinking_sel: sels) {
        select_atoms(crosslinking_sel);
    }



    int sum_is_not_empty_group = 0;
    for (int isel = 0; isel < nsels; isel++) {
        if (sels[isel]->NATOM == 0) {
            cout << "Select group " << isel << " is empty!!" << endl;
        }
        sum_is_not_empty_group += sels[isel]->NATOM;
    }
    

    if (sum_is_not_empty_group == 0) error1.error_exit("ERROR: sels don't contain any atoms!");



    for (unsigned i = 0; i < sels.size(); i++) {
        for (auto &segment:sels[i]->segments_ind) {
	        for (int ind : segment) {
                cout << "Atomname: " << system->atomname[ind] << "; Resname: " << system->resname[ind] << "; index: " << system->atom_index[ind] << " ; segid: " << system->segid[ind] << "; resid: " << system->resid[ind] << " ; Charge: " << system->charge[ind] << endl;
                string patchtype = this->patchtype(system->atomname[ind],system->resname[ind]);
                string segid1 = system->segname[ind];
                int resid1 = system->resid[ind];
                *this->file_temp << "patch " << patchtype << " " << segid1 << ":" << resid1 << endl;
            }
        }
    }




}


ANALYSIS_PASSIVATION::~ANALYSIS_PASSIVATION()
{
    system = NULL;
    //sel1 = NULL;
    //sel2 = NULL;
    this->rdf_count.clear();
//    fclose(this->outfile_box);
    this->seg_res_ids->close();
    this->file_temp->close();
}

