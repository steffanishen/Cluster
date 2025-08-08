/*
 *  read_dcd : c++ class + main file example for reading a CHARMM dcd file
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//***************** Mostly contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************
//Note in general, arrays are faster than vectors.



#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include "array_tools.hpp"

#include "atom.hpp"
#include "dcd_r.hpp"
#include "dcd.hpp"
#include "psf.hpp"
#include "analysis.hpp"
#include "analysis_pointers.hpp"
#include "input.hpp"
#include "analysis_anglezs_rings.hpp"
#include "analysis_mindangles_seg.hpp"

using namespace std;

int main(int argc, char* argv[])
{                
    // instance of a new object DCD_R attached to a dcd file 
    clock_t t;

    string input_file = argv[1];

//    ANALYSIS analysis; //& means the address of sys, exactly what pointer means. The rule of thumb for pointers is use * to access/change the VALUE of the POINTEr, and use & to access the ADDRESS of the VARIABLE.
    INPUT input(input_file);
    cout << "input: " << input_file << endl;


    input.read_parameters();


    string psfname = input.psfname;
    PSF sys(psfname);


    sys.read_psf();

    vector<PSF*> monomers;
    vector<DCD_R*> monomersdcd;

    int monomercount = 0;
    for (string monomer : input.monomers) {
	string monomerpsf = monomer + ".psf";
	string monomerdcd = monomer + ".dcd";
	monomers.push_back(new PSF(monomerpsf));
	monomers[monomercount]->read_psf();
	monomersdcd.push_back(new DCD_R(monomerdcd,monomers[monomercount]));
        monomersdcd[monomercount]->read_header();
        monomersdcd[monomercount]->read_oneFrame();
	monomercount++;
    }

    //for (PSF* monomer: monomers) {
    //    monomer->read_psf();
    //}
    
//******************************* test atom selection ******************************************************************************
// ****************************** selection guideline: need a space between '(' and selection option, a space between seletion argument and ')', and a space between the selection option and any arguments. 
//


    t = clock();
    
    vector<GROUP> sels;
    vector<GROUP*> psels;
    

    for (auto &group_sel : input.group_selection) {
        sels.push_back(sys.select_atoms(group_sel));
    }
 
    for (unsigned i = 0; i < sels.size(); i++) {
 	psels.push_back(&sels[i]);
    }
    
    //begin debugging: output the selection
    for (unsigned i = 0; i < sels.size(); i++) {
        for (auto &segment:psels[i]->segments_ind) {
	        for (int ind : segment) {
     //           cout << "Atomname: " << sys.atomname[ind] << "; Resname: " << sys.resname[ind] << "; Charge: " << sys.charge[ind] << endl; 
            }
        }
    }
    //end debugging

//    sel1 = sys.select_atoms("resname IMP and atomtype NE");

    t = clock() - t;

    cout << "Natoms selected: " << sels[0].atoms.size() << endl;

   
    cout << "before reading dcd" << endl;


    // read the header and print it
    // Note: DCD_R group has to be defined after psf_read, since PSF.x ..., PSF.pbc has to be allocated in psf_read first before assignment in DCD_R.

    vector<DCD_R*> dcdfs;

    int dcdcount = 0;
    int framecount = 0;
    int nframes_temp;
    for (string dcdname : input.dcdnames) {
	dcdfs.push_back(new DCD_R(dcdname,&sys));//MS changed: the first argument specifies the dcd file name. Note this dcdf is a DCD_R object, while dcdf in dcdr.hpp is a fstream object.
        dcdfs[dcdcount]->read_header();
        cout << "finished reading dcd header of " << dcdname << endl;
        dcdfs[dcdcount]->printHeader();
        nframes_temp = dcdfs[dcdcount]->getNFILE(); 
        cout << "Number of frames: " << nframes_temp << endl;
        dcdcount++;
        framecount += nframes_temp; 
    }

    sys.nframes_tot = framecount;

    FILE * outfile_box; 
    outfile_box = fopen("box.dat","w");
    fprintf(outfile_box,"lx   ly   lz\n");
    // in this loop the coordinates are read frame by frame
//    ofstream angle_file ("angles_cpp.dat");
//    ofstream angle_file_all ("angles_cpp_all.dat");

    t = clock();

/////////////////////////////////////////////////Initialization of analysis objects goes here!//////////////////////

    cout << "before assigning pointers" << endl; // Meng debug
    ANALYSIS_POINTERS ana_pointers(&sys,psels,&input,monomers);
    vector<ANALYSIS*> analysis = ana_pointers.init(); 


    cout << "before assigning files" << endl; // Meng debug
    vector<ofstream*> files;    

    for (int compID = 0; compID < analysis.size() ; compID++) {
        ofstream *file_temp = new ofstream (analysis[compID]->filename);
	files.push_back(file_temp);
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

//    cout << "Current index: " << analysis[0]->system->NATOM  << endl;

/*
    int iframe = 0;
    ERROR1 error1;
    for (DCD_R* dcdf:dcdfs) {
        for(int i=0;i<dcdf->getNFILE();i++)
        {

            dcdf->read_oneFrame();
            if (dcdf->NATOM != sys.NATOM) {
                cout << "dcd NATOM: " << dcdf->NATOM << "; psf NATOM: " << sys.NATOM << endl; 
                error1.error_exit("The number of atoms in dcd doesn't match that of psf!!!");
            }

            sys.iframe = iframe;
        // analysis goes here 
            fprintf(outfile_box, "%10.6f %10.6f %10.6f\n", sys.pbc[0], sys.pbc[2], sys.pbc[5]);
            if (iframe == 0) {
               sys.box_first_frame[0] = sys.pbc[0];
               sys.box_first_frame[1] = sys.pbc[2];
               sys.box_first_frame[2] = sys.pbc[5];
            }

	           cout << "frame: " << iframe << " " << sys.x[0] << endl;


///////////////// Analysis, will work on mapping this later///////////////////////////////


//////////////////// We save output for 2D vectors for later ///////////////////////////
            for (int compID = 0; compID < analysis.size() ; compID++) {
	        if (analysis[compID]->vector1d == 1) {
	        vector<float> oned_vector_output = analysis[compID]->compute_vector();
//	        *files[compID] << "frame: " << i << " : ";
                for (int segi = 0; segi < oned_vector_output.size(); segi++) {
	            *files[compID] << oned_vector_output[segi] << " ";
	        }
	        *files[compID] << endl;

	        } else if (analysis[compID]->vector2d == 1) {
	        vector<vector<float>> angles = analysis[compID]->compute_2D_vector();	 	
                    for (int segi = 0; segi < angles.size(); segi++) {
                        for (int segj = 0; segj < angles[segi].size(); segj++)
                        {
	                    *files[compID] << angles[segi][segj] << " ";
                        }
	                *files[compID] << endl;
	            }
	        } else if (analysis[compID]->voidf== 1) {
	            analysis[compID]->compute_void();	 	
	        } 
 	    }
            iframe++;
        }
    }

/////////////////////////////////////////////////////////////////////////////////////
    for (unsigned i = 0; i < sels.size(); i++) {
 	psels[i] = NULL;
    }

*/

//////////////testing matrix multiplication////////////////////////////////
/*
    vector<vector<float>> m1 = {{1.0,2.0,3.0,4.0},{5.0,6.0,7.0,8.0}};
    vector<vector<float>> m2 = {{1.0,2.0},{3.0,4.0},{5.0,6.0},{7.0,8.0}};
    vector<vector<float>> prod = analysis[0]->matrix_mult(m1,m2);
    for (int i=0; i<prod.size(); i++) {
        for (int j=0; j<prod[i].size(); j++) {
            cout << prod[i][j] << " ";
        }
        cout << endl;
    }
*/
//////////////////end testing////////////////////////////////////////////

////////////////testing dot and cross productus/////////////////


//...............end testing////////////////////////////////////
    t = clock() - t;

    cout << "reading dcd cost time: " << ((double) t) / CLOCKS_PER_SEC << " seconds" << endl;
//    fclose(outfile_box);   
 
    return EXIT_SUCCESS;
}
