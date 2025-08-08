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
#include <numeric>

#include <math.h>
#include "analysis_contact_angle.hpp"
#include "analysis_contact_angle_density_profile.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE(PSF *system, GROUP *sel1, GROUP *sel2, int vector1d, int vector2d, int voidf, string filename, string contact_angle_filename, float zshift, float dr, float zlower, string fitting_function, int every_n_frame)
{
    cout << "test initializing droplet density profile: " << endl;
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;


    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->contact_angle_filename = contact_angle_filename;
    contact_angle_file = new ofstream (contact_angle_filename.c_str());
    this->dr = dr;
    this->zlower = zlower;
    this->iframe = 0;

    this->zshift = -zshift;

    this->fitting_function = fitting_function;
    this->every_n_frame = every_n_frame;
    this->nframes = 0;
}

void ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::init() {
}


void ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::polyfit(	const std::vector<double> &t,	const std::vector<double> &v,	std::vector<double> &coeff,	int order)
{
	// Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
	Eigen::MatrixXd T(t.size(), order + 1);
    Eigen::MatrixXd Tp(order+1, t.size());
    Eigen::MatrixXd T2(order+1,order+1);
	Eigen::VectorXd V = Eigen::VectorXd::Map(&v.front(), v.size());
    Eigen::VectorXd Vp;
	Eigen::VectorXd result;

	// check to make sure inputs are correct
	assert(t.size() == v.size());
	assert(t.size() >= order + 1);
	// Populate the matrix
	for(size_t i = 0 ; i < t.size(); ++i)
	{
		for(size_t j = 0; j < order + 1; ++j)
		{
			T(i, j) = pow(t.at(i), j);
		}
	}
	//std::cout<<T<<std::endl;

	// Solve for linear least square fit
//    Tp = T.transpose();
//    T2 = Tp*T;

//    Vp = Tp*V;
	
	result  = T.householderQr().solve(V);
	//result  = T2.householderQr().solve(Vp);

 //   result = T2.inverse()*Vp;
    cout << "size_of_result: " << result.size() << endl;


	coeff.resize(order+1);
	for (int k = 0; k < order+1; k++)
	{
		coeff[k] = result[k];
	}

}

float ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::find_ellipse_x(float yan, vector<double> coeffs) {
    float xan;

    xan = sqrt((1.0 - coeffs[1]*yan -coeffs[2]*yan*yan)/coeffs[0]);

    return xan;
}

vector<float> ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::find_ellipse_y(vector<float> xan, vector<double> coeffs) {
    vector<float> yan_array(xan.size());
    double A = coeffs[0];
    double C = coeffs[1];
    double B = coeffs[2];
    for (int i = 0; i < xan.size(); i++) {
        float x = xan[i];
        yan_array[i] = sqrt((1.0 - A*x*x + C*C*0.25/B)/B) - C*0.5/B;
        float inside_sq = (1.0 - A*x*x + C*C*0.25/B)/B;
        if (inside_sq < 0.0) {
           // cout << "i: " << i << " y: " << x << " inside_sq:" << inside_sq << endl;
        } 
    }
    return yan_array;
}

vector<float> ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::find_ellipse_y_neg(vector<float> xan, vector<double> coeffs) {
    vector<float> yan_array(xan.size());
    double A = coeffs[0];
    double C = coeffs[1];
    double B = coeffs[2];
    for (int i = 0; i < xan.size(); i++) {
        float x = xan[i];
        yan_array[i] = -sqrt((1.0 - A*x*x + C*C*0.25/B)/B) - C*0.5/B;
    }
    return yan_array;

}


vector<float> ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::find_parabola_y(vector<float> xan, vector<double> coeffs) {
    vector<float> yan_array(xan.size());
    double A = coeffs[0];
    double C = coeffs[1];
    double B = coeffs[2];

    cout << "xan_size: " << xan.size() << endl;

    for (int i = 0; i < xan.size(); i++) {
        float x = xan[i];
        yan_array[i] = -sqrt((1.0 - A*x*x + C*C*0.25/B)/B) - C*0.5/B;
     //   cout << "xan.size: " << xan.size() << " i: " << i << " x: " << x << " yan_array: " << yan_array[i] << endl; 
    }

    return yan_array;

}


void ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::output_density(vector<vector<float>> density_yz) {
    for (int iybin=0; iybin < this->ybins; iybin++) {
        for (int izbin=0; izbin < this->zbins; izbin++) {
            float y = this->dr * float(iybin) - this->yshift;
            float z = this->dr * float(izbin) - this->zshift;
            if (density_yz[iybin][izbin] > 0.1*this->density_bulk) *this->density_file << y << " " << z << " " << density_yz[iybin][izbin] * this->every_n_frame << endl;
        }
    }
}

void ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::output_contour(vector<double> y_contour_points, vector<double> z_contour_points) {
            for (int icontour = 0; icontour < y_contour_points.size(); icontour++) {
	            *file_contour_last << y_contour_points[icontour] << " ";
	        }
            *file_contour_last << endl;

            for (int icontour = 0; icontour < z_contour_points.size(); icontour++) {
	            *file_contour_last << z_contour_points[icontour] << " ";
	        }
            *file_contour_last << endl;
}

void ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::output_density_solute(vector<vector<float>> density_yz_solute, vector<double> coeff) {
    density_z_solute.clear();
    density_z_solute.resize(this->zbins,0.0);

    int izbinp;
    float droplet_edge = 5.0;

    float zshift_solute_density = 0.5 * system->box_first_frame[2];
    for (int izbin=0; izbin < this->zbins; izbin++) {
        float density_z = 0.0;
        float npoints = 0.0;

        izbinp = izbin + this->zbins/2 - int(this->zshift/this->dr);

        for (int iybin=0; iybin < this->ybins; iybin++) {
            float y = this->dr * float(iybin) - this->yshift;
            float z = this->dr * float(izbin) - this->zshift;
            //float z = this->dr * float(izbin) + zshift_solute_density;
            if (density_yz[iybin][izbin] > 0.1*this->density_bulk) *this->density_file_solute << y << " " << z << " " << density_yz_solute[iybin][izbin] * this->every_n_frame << endl;

            float yan = find_ellipse_x(z,coeff);

            if (abs(y) < abs(yan) - droplet_edge && izbinp < this->zbins) {
                density_z_solute[izbinp] += density_yz_solute[iybin][izbin] / system->box_first_frame[0];
                npoints += 1.0;
            }
        }
        if (npoints >0.01 && izbinp < this->zbins) density_z_solute[izbinp] /= npoints / 1000.0;
    }
}


vector<float> ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::compute_vector() {
    float x,y,z;
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    float dist2;
    int iybin,izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    float dr = this->dr;
    vector<float> ytemp;
    vector<float> ztemp;
    vector<float> ytemp_sel2;
    vector<float> ztemp_sel2;
    this->iframe += 1;
    this->nframes += 1;
    float R_dense = 3.0;
    float volume_bulk = 0.5*PI*R_dense*R_dense;
    vector<float> frame_angle;
    vector<double> y_contour_points;
    vector<double> z_contour_points;



    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    if (sel2->NATOM == 0) error1.error_exit("ERROR: sel2 doesn't contain any atoms!");

    if (this->iframe == 1){
        vector<double> box;
        box.resize(3);
        box[0] = system->box_first_frame[0];
        box[1] = system->box_first_frame[1]; 
        box[2] = system->box_first_frame[2];

        this->yshift = box[1] * 0.5;
        this->ybins = int(box[1]/dr);
        this->zbins = int(box[2]/dr);
        this->density_yz.resize(this->ybins,vector<float>(this->zbins));
        this->density_yz_solute.resize(this->ybins,vector<float>(this->zbins));
        this->density_z_solute.resize(this->zbins);
    }



    for (auto &segment:sel1->segments_ind) {
      for (int ind : segment) {
	      x = system->x[ind];
	      y = system->y[ind];
          z = system->z[ind];
        
          ytemp.push_back(y);
          ztemp.push_back(z);

        }
      }


    for (auto &segment:sel2->segments_ind) {
      for (int ind : segment) {
	      x = system->x[ind];
	      y = system->y[ind];
          z = system->z[ind];
        
          ytemp_sel2.push_back(y);
          ztemp_sel2.push_back(z);

        }
      }





    float y_com = accumulate(begin(ytemp), end(ytemp), 0.0);
    float z_com = accumulate(begin(ztemp), end(ztemp), 0.0);
    y_com = y_com / float(ytemp.size());
    z_com = z_com / float(ztemp.size());

    this->z_com_frames += z_com;


    for (int i=0; i<ytemp.size(); i++) {
        y = ytemp[i];
        z = ztemp[i];
        iybin = int((y - y_com + yshift)/this->dr);//adjust the y position such that the com of the droplet is in the center of y-axis. Note: ycom can be close to 0 to begin with.
        izbin = int((z + this->zshift)/this->dr);
        if (izbin >=0 && izbin < this->zbins && iybin >=0 && iybin < ybins) {
            density_yz[iybin][izbin] += 1.0;
            if (z > -this->zshift && (z+zshift)*(z+zshift) + (y - y_com)*(y - y_com)< R_dense*R_dense ) {
                this->density_bulk += 1.0;
            }
        }
    }


    for (int i=0; i<ytemp_sel2.size(); i++) {
        y = ytemp_sel2[i];
        z = ztemp_sel2[i];
        iybin = int((y - y_com + yshift)/this->dr);
        izbin = int((z + this->zshift)/this->dr);
        if (izbin >=0 && izbin < this->zbins && iybin >=0 && iybin < ybins) {
            density_yz_solute[iybin][izbin] += 1.0;
            //if (z > -this->zshift && (z+zshift)*(z+zshift) + (y - y_com)*(y - y_com)< R_dense*R_dense ) {
            //    this->density_bulk += 1.0;
            //}
        }
    }



    if (this->iframe % this->every_n_frame == 0 ) {

        this->density_bulk /= (volume_bulk * float(this->every_n_frame));
        this->z_com_frames /= float(this->every_n_frame);

        for (int iybin=0; iybin < this->ybins; iybin++) {
            for (int izbin=0; izbin < this->zbins; izbin++) {
                this->density_yz[iybin][izbin] /= ( this->dr * this->dr * float(this->every_n_frame ) );
                this->density_yz_solute[iybin][izbin] /= ( this->dr * this->dr * float(this->every_n_frame ) );
                float current_dens = this->density_yz[iybin][izbin]; 
                if (current_dens > 0.05*this->density_bulk && current_dens < 0.1*this->density_bulk) {
                    float y_contour = this->dr * float(iybin) - yshift;
                    float z_contour = this->dr * float(izbin) - zshift;
                    if (z_contour > this->zlower) {
                        y_contour_points.push_back(y_contour);
                        z_contour_points.push_back(z_contour);

                    }
                }
            }
        }


        cout << "this->iframe: " << this->iframe << " " <<  this->every_n_frame * (system->nframes_tot/this->every_n_frame) << endl; 
        //output_density(density_yz);

        cout << "remainder: " << this->iframe % this->every_n_frame << " npoints:" << y_contour_points.size() << " bulk_density: " << density_bulk << endl; 
           // cout << "y_com: " << y_com << " npoints: " << y_contour_points.size() << endl; //Meng debug

            for (int icontour = 0; icontour < y_contour_points.size(); icontour++) {
	            *file_temp << y_contour_points[icontour] << " ";
	        }
            *file_temp << endl;

            for (int icontour = 0; icontour < z_contour_points.size(); icontour++) {
	            *file_temp << z_contour_points[icontour] << " ";
	        }
            *file_temp << endl;

            vector<double> y_contour_points2(y_contour_points.size());
            vector<double> z_contour_points2(z_contour_points.size());
            for (int icontour = 0; icontour < y_contour_points.size(); icontour++) {
	            y_contour_points2[icontour] = y_contour_points[icontour]*y_contour_points[icontour]; 
	            z_contour_points2[icontour] = z_contour_points[icontour]*z_contour_points[icontour]; 
	        }

        int fitting_order;

        double small_double = 0.00000001;
        std::vector<double> coeff_temp;

        int order_general = 3;

        if (z_contour_points.size() >= order_general) {

            coeff.clear();
            if (this->fitting_function == "ellipse") {

                fitting_order = 2;
	            polyfit(z_contour_points, y_contour_points2, coeff_temp, fitting_order);
                coeff.resize(fitting_order+1,0.0);

                if (coeff_temp[0]>-small_double && coeff_temp[0] < small_double ) {
                    error1.error_exit("coeff0 is 0!"); 
                }

                for (int icoeff = 1; icoeff < coeff_temp.size(); icoeff++) {
                    coeff[icoeff] = -coeff_temp[icoeff] / coeff_temp[0];
            //coeff[icoeff] /= -coeff[0];
                }

            } else if (this->fitting_function == "sphere") {
                vector<double> y2_plus_z2 = vector_sum(y_contour_points2, z_contour_points2);
                //vector<double> y2_plus_z2 =  y_contour_points2 + z_contour_points2; 
                fitting_order = 1;
	            polyfit(z_contour_points, y2_plus_z2, coeff_temp, fitting_order);
                coeff.resize(fitting_order+2,0.0);

                if (coeff_temp[0]>-small_double && coeff_temp[0] < small_double ) {
                    error1.error_exit("coeff0 is 0!"); 
                }

                for (int icoeff = 1; icoeff < coeff_temp.size(); icoeff++) {
                    coeff[icoeff] = -coeff_temp[icoeff] / coeff_temp[0];
                //coeff[icoeff] /= -coeff[0];
                }
                coeff[2] = 1.0/coeff_temp[0];

            } else {
                error1.error_exit("please specify the fitting_function. Options: 1. ellipse; 2. sphere");
            }

            coeff[0] = 1.0/ coeff_temp[0]; 
     //   float swap_coeff = coeff[1];
     //   coeff[1] = coeff[2];
     //   coeff[2] = swap_coeff;



            cout << "coeffs linear: ";
            for (int icoeff = 0; icoeff < coeff_temp.size(); icoeff++) {
                cout << coeff_temp[icoeff] << " ";
            }
            cout << endl;

            cout << "coeffs: ";
            for (int icoeff = 0; icoeff < coeff.size(); icoeff++) {
                cout << coeff[icoeff] << " ";
            }
            cout << endl;

            float yan = -find_ellipse_x(-this->zshift,coeff);
          //  float inside_sq = sqrt((1.0 + coeff[1]*coeff[1]*0.25/coeff[2])/coeff[0]);

          //  cout << "yan: " << yan << " inside_sq: " << inside_sq << endl;

            vector<float> yan_array = this->linspace(yan,-yan,0.1);
            vector<float> zan_array(yan_array.size());

            cout << "yan_size: " << yan_array.size() <<  " first_yan_array: "<< yan_array[0] << endl;

            if (coeff[0] > 0.0 && coeff[2] > 0.0 ) {
                cout << "A>0, B>0" << endl;
                zan_array = find_ellipse_y(yan_array,coeff);
            } else if (coeff[0] > 0.0 && coeff[2] < 0.0) {
                cout << "A>0, B<0" << endl;
                zan_array = find_parabola_y(yan_array,coeff);
            }

            float slope = (zan_array[1] - zan_array[0] )/(yan_array[1] - yan_array[0] );

            this->contact_angle = atan(slope) * 180.0 / PI;

            if (coeff[0] > 0.0 && coeff[2] > 0.0 && coeff[1] < 0.0 ) {
                this->contact_angle = 180.0 - this->contact_angle;
            }

            this->nframes = 0;

            //*contact_angle_file << float(this->iframe) << " " << this->contact_angle << " " << this->z_com_frames << endl;
            *contact_angle_file << float(this->iframe) << " " << this->contact_angle << " " << z_com << endl;
            this->z_com_frames = 0.0;

        } 

/*
    //fitting benchmark begin

	// time value
	std::vector<double> time = {0, 0.0192341804504395, 0.0394501686096191,  0.059575080871582, 0.0790810585021973, 0.0792751312255859, 0.0987141132354736,  0.119336366653442,  0.138712167739868,  0.159000158309937,  0.178890228271484,   0.19960618019104,  0.219112157821655,   0.23919415473938,  0.259442090988159,  0.279186248779297,  0.299112319946289,  0.319219350814819,  0.339494228363037,  0.339675188064575,  0.359552145004272,   0.37941837310791,  0.399189233779907,  0.419828176498413,  0.439810276031494,  0.459331274032593,  0.479461193084717,  0.499663114547729,  0.519809246063232,  0.539092063903809,  0.559118270874023,  0.579315185546875,  0.598889112472534,  0.619685173034668,  0.638863086700439,  0.639052152633667,  0.658920288085938,  0.679149150848389,  0.699787139892578,   0.71905517578125,   0.73898720741272,  0.739143371582031,  0.758654117584229,  0.779210329055786,  0.799195289611816,  0.819046258926392,  0.839539289474487,   0.85923433303833,   0.87903618812561,  0.899263143539429,  0.919251203536987,  0.939138174057007,  0.959244251251221,  0.979074239730835,  0.998935222625732,   1.01904726028442,    1.0387852191925,   1.03895926475525,   1.05906510353088,   1.07873225212097,   1.09908628463745,   1.11907029151917,   1.13899827003479,   1.15879201889038};
	// velocity value
	std::vector<double> velocity = {1.8, 1.86, 2.03, 2.08, 2.14, 2.14, 2.25, 2.36, 2.42, 2.59,  2.7, 2.81, 2.87, 3.04, 3.15, 3.26, 3.32, 3.43, 3.54, 3.54,  3.6, 3.71, 3.83, 3.94, 4.11, 4.22, 4.33, 4.44, 4.56, 4.67, 4.78, 4.84, 4.84, 4.89, 4.89, 4.89, 4.95, 5.01, 5.06, 5.06, 5.06, 5.06, 5.01, 5.06, 5.12, 5.18, 5.18, 5.23, 5.23, 5.23, 5.29, 5.34, 5.29,  5.4,  5.4, 5.46, 5.51, 5.51, 5.51, 5.46,  5.4, 5.34, 5.34, 5.34};

	// placeholder for storing polynomial coefficient
        std::vector<double> coeff;
	polyfit(time, velocity, coeff, 3);

	std::vector<double> fitted_velocity;
	std::cout<< "Printing fitted values" << std::endl;
	for(int p = 0; p < time.size(); ++ p)
	{
		double vfitted = coeff[0] + coeff[1]*time.at(p) + coeff[2]*(pow(time.at(p), 2)) +coeff[3]*(pow(time.at(p), 3)) ;
		std::cout<< vfitted<<", ";
		fitted_velocity.push_back(vfitted);
	}
	std::cout<<std::endl;

// Fitting benchmark end

*/
        if (this->iframe == this->every_n_frame * (system->nframes_tot/this->every_n_frame)) {
            output_density(density_yz);
            output_contour(y_contour_points, z_contour_points);
        }


        frame_angle.push_back(float(this->iframe));
        frame_angle.push_back(this->contact_angle);
        frame_angle.push_back(z_com);

        output_density_solute(density_yz_solute,this->coeff);

        density_yz.clear();
        density_yz_solute.clear();

        this->density_yz.resize(this->ybins,vector<float>(this->zbins,0.0));
        this->density_yz_solute.resize(this->ybins,vector<float>(this->zbins,0.0));
        this->density_bulk = 0.0;

    }

    return density_z_solute;
}


ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE::~ANALYSIS_CONTACT_ANGLE_DENSITY_PROFILE()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
    this->density_yz.clear();
    this->density_yz_solute.clear();
    this->density_z_solute.clear();
    this->density_file->close();
}

