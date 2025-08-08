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
//***************** use of this class must follow anglezs_rings **************************

#include <iostream>
#include <vector>
#include "psf.hpp"
#include "atom.hpp"
#include "group.hpp"
#include "analysis.hpp"
#include <Eigen/Dense>
#include <Eigen/QR>
#ifndef ANALYSIS_CONTACT_ANGLE_HPP
#define	ANALYSIS_CONTACT_ANGLE_HPP

using namespace std;

class  ANALYSIS_CONTACT_ANGLE : public ANALYSIS
{

private:
    //private attributes
    //no private methods
    
public:
   float contact_angle = 0.0;
    // no public attributes
    // public methods
//    GROUP *sel1;
 //   int whichN;
    //
    float yshift;
    float zshift;
    int ybins;
    int zbins;
    int zlower;
    int every_n_frame;
    int nframes;

    int iframe;
    float z_com_frames = 0.0;

    float density_bulk = 0.0;
    //vector<double> coeff;
    
    ofstream *file_temp = new ofstream ("contour.dat");
    ofstream *file_contour_last = new ofstream ("contour_last.dat");
    ofstream *contact_angle_file;
    ofstream *density_file = new ofstream("density_2D.dat");

    vector<vector<float>> density_yz;
    string contact_angle_filename;
    string fitting_function;

    ANALYSIS_CONTACT_ANGLE(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, string contact_angle_filename, float zshift, float dr, float zlower, string fitting_function, int every_n_frame); //constructor
    
    void init();

    vector<float> compute_vector();
    float find_ellipse_x(float yan, vector<double> coeffs);
    vector<float> find_ellipse_y(vector<float> xan, vector<double> coeffs);
    vector<float> find_ellipse_y_neg(vector<float> xan, vector<double> coeffs);
    vector<float> find_parabola_y(vector<float> xan, vector<double> coeffs);
    void output_density(vector<vector<float>> density_yz);
    void output_contour(vector<double> y_contour, vector<double> z_contour);

    void polyfit(	const std::vector<double> &t, const std::vector<double> &v, std::vector<double> &coeff,	int order);

    ~ANALYSIS_CONTACT_ANGLE();

};

#endif	/* DCD_R_HPP */

