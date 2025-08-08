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
#include "analysis_pointers.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_POINTERS::ANALYSIS_POINTERS(PSF *system, vector<GROUP*> sels, INPUT* input, vector<PSF*> monomers)
{
    this->system = system;

    for (unsigned i = 0; i < sels.size(); i++) {
        this->sels.push_back(sels[i]);
    }
    this->input = input;

    for (unsigned i = 0; i < monomers.size(); i++) {
        this->monomers.push_back(monomers[i]);
    }
}

vector<ANALYSIS*> ANALYSIS_POINTERS::init() {
    for (auto &analysis_opt : input->analysis) {
        int vector1d;
 	int vector2d;
 	int voidf;
	string filename;
	if (analysis_opt[0] == "anglezs_rings") {
		int groupid;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  
		}
		analysis.push_back(new ANALYSIS_ANGLEZS_RINGS(system,sels[groupid],vector1d,vector2d,voidf,filename));

        } else if (analysis_opt[0] == "mindangles_seg_dislocated") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    } else if (analysis_opt[argid] == "whichN") {
			whichN = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_MINDANGLES_SEG_DISLOCATED(system,sels[groupid],whichN,vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "scale") {
		int groupid;
		float scaling;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "scaling") {
			scaling = stoi(analysis_opt[argid+1]);
		    }  
		}
	        analysis.push_back(new ANALYSIS_SCALE(system,sels[groupid],vector1d,vector2d,voidf,filename,scaling));
        } else if (analysis_opt[0] == "mindangles_seg") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_MINDANGLES_SEG(system,sels[groupid],vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "avedangles_seg") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_AVEDANGLES_SEG(system,sels[groupid],vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "density_profile") {
		int groupid;
	 	int whichN;
                int nbins;
                string which_density_profile;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "nbins") {
			nbins = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "which_density_profile") {
			which_density_profile = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_DENSITY_PROFILE(system,sels[groupid],vector1d,vector2d,voidf,filename,nbins,which_density_profile));

        } else if (analysis_opt[0] == "charge_density_profile") {
		int groupid;
	 	int whichN;
                int nbins;
                string which_density_profile;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "nbins") {
			nbins = stoi(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_CHARGE_DENSITY_PROFILE(system,sels[groupid],vector1d,vector2d,voidf,filename,nbins,which_density_profile));









        } else if (analysis_opt[0] == "density_profile_droplet") {
		    int groupid;
	 	    int whichN;
            float zshift;
            float dr;
            int nbins;
            string which_density_profile;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "nbins") {
			nbins = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "which_density_profile") {
			which_density_profile = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "zshift") {
			zshift = stof(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "dr") {
			dr = stof(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_DENSITY_PROFILE_DROPLET(system,sels[groupid],vector1d,vector2d,voidf,filename,nbins,which_density_profile,zshift,dr));


        } else if (analysis_opt[0] == "nshell") {
		int groupid;
		int groupid1;
	 	int whichN;
                int nbins;
                float dist_crit;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    } else if (analysis_opt[argid] == "group1") {
			groupid1 = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "nbins") {
			nbins = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "dist_crit") {
			dist_crit = stof(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_SHELL(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,nbins,dist_crit));

        } else if (analysis_opt[0] == "nshell_coord") {
		int groupid;
		int groupid1;
	 	int whichN;
                int nbins;
                float dist_crit;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    } else if (analysis_opt[argid] == "group1") {
			groupid1 = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "nbins") {
			nbins = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "dist_crit") {
			dist_crit = stof(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_SHELL_OUTPUT_COORD(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,nbins,dist_crit));




        } else if (analysis_opt[0] == "inter_energy") {
		int groupid;
		int groupid1;
	 	int whichN;
                float dist_crit;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    } else if (analysis_opt[argid] == "group1") {
			groupid1 = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "dist_crit") {
			dist_crit = stof(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_ENERGY(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,dist_crit));
        } else if (analysis_opt[0] == "rdf") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                float dr;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "dist_crit") {
                        dist_crit = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_RDF(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,dist_crit,dr));

        
        } else if (analysis_opt[0] == "rdf_region") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                float dr;
                float zlow;
                float zhigh;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "dist_crit") {
                        dist_crit = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "zlow") {
                        zlow = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "zhigh") {
                        zhigh = stof(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_RDF_REGION(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,dist_crit,dr,zlow,zhigh));

        
        } else if (analysis_opt[0] == "orientation") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                int nbins;
                float dtheta;
                string name0;
                string name1;
                string name2;
                string name3;
                int every_n_frame;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name0") {
                        name0 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name1") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name2") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name3") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "nbins") {
                        nbins = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "every_n_frame") {
                        every_n_frame = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dtheta") {
                        dtheta = stof(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_ORIENTATION(system,sels[groupid],vector1d,vector2d,voidf,filename,name0,name1,name2,name3,nbins,every_n_frame,dtheta));


        } else if (analysis_opt[0] == "orientation_chain") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                int nbins;
                float dtheta;
                int every_n_frame;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "nbins") {
                        nbins = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "every_n_frame") {
                        every_n_frame = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dtheta") {
                        dtheta = stof(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_ORIENTATION_CHAIN(system,sels[groupid],vector1d,vector2d,voidf,filename,nbins,every_n_frame,dtheta));







        } else if (analysis_opt[0] == "contact_angle") {
                int groupid;
                int groupid1;
                int whichN;
                float zshift;
                float zlower;
                float dr;
                int every_n_frame;
                string contact_angle_filename;
                string fitting_function;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "contact_angle_filename") {
                        contact_angle_filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "zshift") {
                        zshift = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "zlower") {
                        zlower = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "every_n_frame") {
                        every_n_frame = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "fitting_function") {
                        fitting_function = analysis_opt[argid+1];
                    }
                }
                analysis.push_back(new ANALYSIS_CONTACT_ANGLE(system,sels[groupid],vector1d,vector2d,voidf,filename,contact_angle_filename,zshift,dr,zlower,fitting_function,every_n_frame));


        } else if (analysis_opt[0] == "contact_angle_density_profile") {
                int groupid;
                int groupid1;
                int whichN;
                float zshift;
                float zlower;
                float dr;
                int every_n_frame;
                string contact_angle_filename;
                string fitting_function;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "contact_angle_filename") {
                        contact_angle_filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "zshift") {
                        zshift = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "zlower") {
                        zlower = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "every_n_frame") {
                        every_n_frame = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "fitting_function") {
                        fitting_function = analysis_opt[argid+1];
                    }
                }
                analysis.push_back(new ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,contact_angle_filename,zshift,dr,zlower,fitting_function,every_n_frame));




        } else if (analysis_opt[0] == "trajectory_extract") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                float dr;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }
                }
                analysis.push_back(new ANALYSIS_TRAJECTORY_EXTRACT(system,sels[groupid],vector1d,vector2d,voidf,filename));




        } else if (analysis_opt[0] == "msd") {
		int groupid;
	 	int whichN;
                int dtmax;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "dtmax") {
			dtmax = stoi(analysis_opt[argid+1]);
		    }
		}
        	analysis.push_back(new ANALYSIS_MSD(system,sels[groupid],vector1d,vector2d,voidf,filename,dtmax));

        } else if (analysis_opt[0] == "functionalize") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                float dr;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "dist_crit") {
                        dist_crit = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_FUNCTIONALIZE(system,sels[groupid],sels[groupid1],vector1d,vector2d,voidf,filename,dist_crit,dr));
        } else if (analysis_opt[0] == "patch_no_order") {
                int groupid;
                int groupid1;
                int whichN;
                float dist_crit;
                float dr;
                float cellsizex;
                float cellsizey;
                float cellsizez;
                int crosslink_in_cluster;
                int pre_crosslinking;
                string input_cluster;
                string output_cluster;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "group1") {
                        groupid1 = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "input_cluster") {
                        input_cluster = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "output_cluster") {
                        output_cluster = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "dist_crit") {
                        dist_crit = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "dr") {
                        dr = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "cellsizex") {
                        cellsizex = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "cellsizey") {
                        cellsizey = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "cellsizez") {
                        cellsizez = stof(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "crosslink_in_cluster") {
                        crosslink_in_cluster = stoi(analysis_opt[argid+1]);
                    } else if (analysis_opt[argid] == "pre_crosslinking") {
                        pre_crosslinking = stoi(analysis_opt[argid+1]);
                    }
                }
                analysis.push_back(new ANALYSIS_PATCH_NO_ORDER(system,sels[groupid],sels[groupid1],sels,vector1d,vector2d,voidf, input_cluster, output_cluster,filename,dist_crit,dr,cellsizex,cellsizey,cellsizez,crosslink_in_cluster,pre_crosslinking));
} else if (analysis_opt[0] == "passivation") {
                int groupid;
                int groupid1;
                int whichN;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                     if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }
                }
                analysis.push_back(new ANALYSIS_PASSIVATION(system,sels,vector1d,vector2d,voidf, filename));
        } else if (analysis_opt[0] == "msd_region") {
		int groupid;
	 	int whichN;
                int dtmax;
                float low;
                float high;
                int msd_xflag = 0;
                int msd_yflag = 0;
                int msd_zflag = 0;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "dtmax") {
			dtmax = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "x") {
                        msd_xflag = 1;
			low = stof(analysis_opt[argid+1]);
			high = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "y") {
                        msd_yflag = 1;
			low = stof(analysis_opt[argid+1]);
			high = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "z") {
                        msd_zflag = 1;
			low = stof(analysis_opt[argid+1]);
			high = stoi(analysis_opt[argid+2]);
		    }
		}
        	analysis.push_back(new ANALYSIS_MSD_REGION(system,sels[groupid],vector1d,vector2d,voidf,filename,dtmax,msd_xflag,msd_yflag,msd_zflag,low,high));

        } else if (analysis_opt[0] == "random_walk") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
		float fringe;
	        float bl;
                int wall_rw;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "wall_rw") {
                        wall_rw = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "fringe") {
                        fringe = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_RANDOM_WALK(system,sels[groupid],vector1d,vector2d,voidf,bl,filename,wall_rw,monomers,fringe));
        } else if (analysis_opt[0] == "straight_chain") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
		float fringe;
	        float bl;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "fringe") {
                        fringe = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_STRAIGHT_CHAIN(system,sels[groupid],vector1d,vector2d,voidf,bl,filename,monomers,fringe));
        } else if (analysis_opt[0] == "packing") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
	        float lx;
	        float ly;
	        float lz;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "lx") {
                        lx = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "ly") {
                        ly = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "lz") {
                        lz = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_PACKING(system,sels[groupid],vector1d,vector2d,voidf,filename,lx,ly,lz));
        } else if (analysis_opt[0] == "rotate") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
	        float bl;
		int nbinsangle;
		int axisid;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "nbinsangle") {
                        nbinsangle = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "axisid") {
                        axisid = stoi(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_ROTATE(system,sels[groupid],vector1d,vector2d,voidf,bl,name1,name2,name3,name4,name_ref1,name_ref2,filename,nbinsangle,axisid));
        } else {
   error1.error_exit("Wrong analysis option, please check the input script!") ;
        }
    }
    return analysis;
}


ANALYSIS_POINTERS::~ANALYSIS_POINTERS()
{
    for (unsigned i = 0; i < analysis.size(); i++) {
	analysis[i] = NULL;
    }
    analysis.clear() ;
    
    for (unsigned i = 0; i < sels.size(); i++) {
        sels[i] = NULL;
    }
    sels.clear();

    system = NULL;
    input = NULL;

}

