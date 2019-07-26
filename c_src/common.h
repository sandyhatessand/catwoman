/* The batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/* Must be defined in the C file that includes this header. */
inline double intensity(double x, double* args);

inline double area(double d, double x, double R, double theta)
{
	/*
	Returns area of overlapping circles with radii x and R; separated by a distance d
	*/
	
	int sit = 0;
	double A = (x*x-R*R-d*d)/(2*d);
	double B = A*(2*cos(theta)*cos(theta) - 1) - 2*cos(theta)*sin(fabs(theta))*sqrt(R*R - A*A);	
	double x_minpos = -d*cos(theta)+sqrt(x*x-d*d*sin(theta)*sin(theta));
	double x_minneg = -d*cos(theta)-sqrt(x*x-d*d*sin(theta)*sin(theta));
	double a1 = -d*cos(theta) + sqrt(x*x - d*d*sin(theta)*sin(theta));
	double b1 = 0;
	double a2 = A*cos(theta)+sin(fabs(theta))*sqrt(R*R - A*A);
	double b2 = sqrt(R*R*cos(theta)*cos(theta) - A*B);
	double a3 = -d*cos(theta) + x;
	double b3 = d*sin(fabs(theta));

	if (fabs(d * sin(theta)) < x) {
		if (fabs(x_minpos)<=R && fabs(x_minneg)>R) {
			if (b3>=b2){
				sit = 1;    //(a-1)
			} else if ((b3<b2) && (b3>b1) && (a2>=a1)){
				sit = 2;   //(a-2)
			} else if ((b3<b2) && (b3>b1) && (a2<a1)){
				sit = 3;    //(a-3)
			} else if ((b3<=b1) && (a2<a1)){
				sit = 4;  //(b)
			}
		} else if (fabs(x_minpos)<=R && fabs(x_minneg)<=R){
			if (A*A>R*R)||(((A*sin(fabs(theta))+cos(theta)*sqrt(R*R-A*A))<0) && ((A*sin(fabs(theta))-cos(theta)*sqrt(R*R-A*A))<0)) {
				sit = 5;  //(c-1)
			} else if (((A*sin(fabs(theta))+cos(theta)*sqrt(R*R-A*A))>=0)&& ((A*sin(fabs(theta))-cos(theta)*sqrt(R*R-A*A))>=0)) {
				sit = 7;  //(c-3)
			}
		} else if (fabs(x_minpos)>R && fabs(x_minneg)>R){
			if (((A*sin(fabs(theta))+cos(theta)*sqrt(R*R-A*A))>=0) && ((A*sin(fabs(theta))-cos(theta)*sqrt(R*R-A*A))>=0)) {
				sit = 6;   //(c-2)
			}
		}
	} else if (d>=x+R) {
		sit = 8;  //no intersection
	} else if ((d<x+R) && (d+x>=R)) {
		sit = 9;  //intersection of two circles
	} else if ((d<x+R) && (d+x<R)) {
		sit = 10;  //circle completely inside semi-circle
	}
			
	switch(sit) {
		case 1 :   //(a-1)
			double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(x*x-(a1+d*cos(theta))*(a1+d*cos(theta)))+asin((a1+d*cos(theta))/x);
			double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(x*x-(a2+d*cos(theta))*(a2+d*cos(theta)))+asin((a2+d*cos(theta))/x);			
			double delta_f = f_a1 - f_a2;
			double g_a2 = acos(a2/R) - (a2/R)*sqrt(1-(a2/R)*(a2/R));
			
			return M_PI*R*R/2 - (x*x/2)*delta_f - d*sin(fabs(theta))*(a2-a1) - (R*R/2)*g_a2;
			break;
		case 2 :    //(a-2) 	
			double h_a1 = acos(a1/R) + (a1/R)*sqrt(1-(a1/R)*(a1/R));
			double p_a1 = asin(a1/R)+(a1/R)*sqrt(1-a1*a1/(R*R));
			double p_a2 = asin(a2/R)+(a2/R)*sqrt(1-a2*a2/(R*R));
			double delta_p = p_a2 - p_a1;
			double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));
			
			return (R*R/2)*(h_a1+delta_p)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha));
			break;
		case 3 :   //(a-3)
			double h_a2 = acos(a2/R) + (a2/R)*sqrt(1-(a2/R)*(a2/R));
			double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));

			return (R*R/2)*h_a2-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha));
			break;
		case 4 :     //(b)
			double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(x*x-(a1+d*cos(theta))*(a1+d*cos(theta)))+asin((a1+d*cos(theta))/x);
                        double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(x*x-(a2+d*cos(theta))*(a2+d*cos(theta)))+asin((a2+d*cos(theta))/x);
                        double delta_f = f_a1 - f_a2;
			double h_a2 = acos(a2/R) + (a2/R)*sqrt(1-(a2/R)*(a2/R));
			
			return (R*R/2)*h_a2 + (x*x/2)*delta_f + d*sin(fabs(theta))*(a2-a1);
			break
		case 5 :    //(c-1)
			double w = sqrt(x*x-d*d*sin(theta)*sin(theta));
			double z = d*sin(fabs(theta));

			return x*x*acos(z/x) - w*z;
			break;
		case 6 :  //(c-2)
			double u = (d*d+x*x-R*R)/(2*d*x);
			double v = (d*d+R*R-x*x)/(2*d*R);
			double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
			double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);

			return A_int - M_PI*R*R/2;
			break;
		case 7 : //(c-3)
			double w = sqrt(x*x-d*d*sin(theta)*sin(theta));
                        double z = d*sin(fabs(theta));
			double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);
			double A1 = x*x*acos(z/x) - w*z;

			return A_int - A1;
			break;
		case 8 : //no intersection
			return 0;
			break;
		case 9 : //intersection of two circles
			double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);

			return A_int;
			break;
		case 10 : //circle completely inside semi-circle
			return M_PI*x*x;
			break;

		default : 
		        break;
	}


void calc_limb_darkening(double* f_array, double* d_array, int N, double rprs, double fac, int nthreads, double* intensity_args, double phi, double b, double min_i)
{
	/*
		This function takes an array of sky distances (d_array) of length N, computes stellar intensity by calling intensity with
		intensity_args, and puts the results in f_array.  To use this function, include this file in a .c file and implement
		the intensity function within that .c file.

		The proper way of implementing this function is to accept a pointer to the intensity function.  Unfortunately, few
		compilers that implement OpenACC support function pointers, so this approach is not yet possible.
	*/

	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	#if defined (_OPENACC)
	#pragma acc parallel loop copyout(f_array[:N]) present(intensity_args)
	#elif defined (_OPENMP)
	#pragma omp parallel for
	#endif
	for(int i = 0; i < N; i++)
	{
		double d = d_array[i];
		double x_in = MAX(d - rprs, 0.);					//lower bound for integration
		double x_out = MIN(d + rprs, 1.0);					//upper bound for integration

		if(x_in >= 1.) f_array[i] = 1.0;					//flux = 1. if the planet is not transiting
		else if(x_out - x_in < 1.e-7) f_array[i] = 1.0;				//pathological case	
		else
		{
			double delta = 0.;						//variable to store the integrated intensity, \int I dA
			double x = x_in;						//starting radius for integration
			double dx = fac*acos(x); 					//initial step size

			x += dx;						//first step

			double A_i = 0.;						//initial area
			
			double theta;
			
			if (phi >= 0. && i < min_i){
				theta = asin(b/d) + phi;
			} else if (phi >= 0 && i >= min_i) {
				theta = asin(b/d) - phi;
	                } else if (phi <= 0 && i < min_i){
				theta = asin(b/d) - phi; 
			} else if (phi <= 0 && i >= min_i){
				theta = asin(b/d) + phi;
			}

			while(x < x_out)
			{
				double A_f = area(d, x, rprs, theta);				//calculates area of overlapping circles
				double I = intensity(x - dx/2., intensity_args); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dx = fac*acos(x);  				//updating step size
				x = x + dx;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dx = x_out - x + dx;  					//calculating change in radius for last step  FIXME
			x = x_out;						//final radius for integration
			double A_f = area(d, x, rprs, theta);					//area for last integration step
			double I = intensity(x - dx/2., intensity_args); 		//intensity at the midpoint
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step

			f_array[i] = 1.0 - delta;	//flux equals 1 - \int I dA
		}
	}
}
