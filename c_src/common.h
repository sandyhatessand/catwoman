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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Must be defined in the C file that includes this header. */
inline double intensity(double x, double* args);

double area(double d, double x, double R, double theta);


double area(double d, double x, double R, double theta)
{
        /*
        Returns area of an overlapping semi-circle with radii R and circle x; separated by a distance d
        */ 
/*	double arg1 = (d*d + x*x - R*R)/(2.*d*x);
	double arg2 = (d*d + R*R - x*x)/(2.*d*R);
	double arg3 = MAX((-d + x + R)*(d + x - R)*(d - x + R)*(d + x + R), 0.);
	if(x <= R - d){
		return M_PI*x*x;							//planet completely overlaps stellar circle
	} else if(x >= R + d){
		return M_PI*R*R;						//stellar circle completely overlaps planet
	} else {
		return x*x*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3);	
	}
*/	
	//define and initialise variables
//  	printf("d = %.20f\n x = %.20f\n R = %.20f\n theta = %.20f\n", d,x,R,theta);     
	int sit = 12;
        double A = 1, x_minpos=1,x_minneg=1, a1=1,b1=1,a2=1,b2=1,a3=1,b3=1;

	A = (x*x-R*R-d*d)/(2*d);
        x_minpos = (-1)*d*cos(theta)+sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
	x_minneg = (-1)*d*cos(theta)-sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
//	printf("x_minppos = %.20f\n x_minneg = %.20f\n", x_minpos, x_minneg);
	a1 = (-1)*d*cos(theta) + sqrt(fabs(x*x - d*d*sin(theta)*sin(theta)));
        b1 = 0.0;
        b2 = -A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A));  //theta should be left raw here (no fabs)
	a2 = b2*tan(theta)+A/cos(theta);      //theta should be left raw here (no fabs)
        a3 = (-1)*d*cos(theta) + x;
        b3 = d*sin(theta);   //theta should be left raw here (no fabs)
	double lim = pow(10,-7);
//	printf("lim = %.20f\n a1 = %.20f\n a2 = %.20f\n b1 = %.20f\n b2 = %.20f\n a3 = %.20f\n b3 = %.20f\n", lim,a1,a2,b1,b2,a3,b3);   
       	
	if ((fabs(d * sin(theta)) < x)&&(fabs(d*sin(theta)-x)>lim)) {       
               	if ((fabs(x_minpos)<R||(fabs(x_minpos-R)<lim)) && fabs(x_minneg)>R && fabs(x_minneg-R)>=lim) {  //implies circle is to the left of semi-circle
                        if (theta>0||(fabs(theta)<lim)){
				if (b3>b2||(fabs(b3-b2)<lim)){
                               		sit = 1;    //(a-1)
                       		} else if ((b3<b2)&&(fabs(b3-b2)>=lim) && (b3>b1||(fabs(b3-b1)<lim)) && (a2>a1||(fabs(a1-a2)<lim))){
                               		sit = 2;   //(a-2)
                       		} else if ((b3<b2)&&(fabs(b3-b2)>=lim) && (b3>b1||(fabs(b3-b1)<lim)) && (a2<a1)&&(fabs(a1-a2)>=lim)){
                               		sit = 3;    //(a-3)
				}
			} else if (theta<0||(fabs(theta)<lim)){ 
                       		if ((b3<b1||(fabs(b3-b1)<lim)) && (a2<a1||(fabs(a1-a2)<lim))){
                               		sit = 4;  //(b)
                	       	} 
			}
               	}else if ((fabs(x_minpos)<R||(fabs(x_minpos-R)<lim)) && (fabs(x_minneg)<R||(fabs(x_minneg-R)<lim))){

                       	if ((((-A*sin((theta))+cos(theta)*sqrt(fabs(R*R-A*A)))>0)||(fabs(-A*sin((theta))+cos(theta)*sqrt(fabs(R*R-A*A)))<lim)) && (((-A*sin((theta))-cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin((theta))-cos(theta)*sqrt(fabs(R*R-A*A))))<lim)) {

                               	sit = 7;  //(c-3) not sure if there should be negative sign in front of A...

                       	} else if ((A*A>R*R)||(((A*sin((-theta))+cos(theta)*sqrt(fabs(R*R-A*A)))<0) && ((A*sin(((-theta)))-cos(theta)*sqrt(fabs(R*R-A*A)))<0))) {
				sit = 5; //(c-1)
			}
               	}else if (fabs(x_minpos)>R && fabs(x_minneg)>R){

			if ((((A*sin(-(theta))+cos(theta)*sqrt(R*R-A*A))>0)||(fabs(A*sin(-(theta))+cos(theta)*sqrt(R*R-A*A))<lim)) && (((A*sin(-(theta))-cos(theta)*sqrt(R*R-A*A))>=0)||(fabs(A*sin(-(theta))-cos(theta)*sqrt(R*R-A*A))<lim))) {

                               	sit = 6;   //(c-2)                

			} else if ((A*A>R*R)||(((A*sin((-theta))+cos(theta)*sqrt(fabs(R*R-A*A)))<0) && ((A*sin(((-theta)))-cos(theta)*sqrt(fabs(R*R-A*A)))<0))) {                         
				sit = 11; //semi-circle completely inside circle
			} else if ((d>x+R)||(fabs(d-(x+R))<lim)) {
				sit = 8;  //no intersection
			}
		}
        } else if ((d>x+R)||(fabs(d-(x+R))<lim)) {    //this is wrong
               	sit = 8;  //no intersection
       	} else if ((d+x>R)||(fabs(d+x-R)<lim)) {
       		sit = 9;  //intersection of two circles
        } else if (d+x<R) {
                sit = 10;  //circle completely inside semi-circle
	}
	
//	printf("sit =  %d\n", sit);	

	switch(sit) {
                case 1:{   //(a-1)
			double h_a2 = asin(a2/R) + (a2/R)*sqrt(1-(a2/R)*(a2/R));
			double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+asin((a1+d*cos(theta))/x);
                        double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+asin((a2+d*cos(theta))/x);
                        double delta_f = f_a1 - f_a2;
//                      printf("delta_f = %.20f\n  x*x-(a1+d*cos(theta))*(a1+d*cos(theta)) = %.20f\n x*x-(a2+d*cos(theta))*(a2+d*cos(theta)) = %.20f\n",delta_f,x*x-(a1+d*cos(theta))*(a1+d*cos(theta)),x*x-(a2+d*cos(theta))*(a2+d*cos(theta)));
			double ret = (M_PI*R*R/2 - (x*x/2)*delta_f - d*sin(fabs(theta))*(a2-a1) - (R*R/2)*(M_PI/2-h_a2));
			if (isnormal(ret)||ret == 0){
				return ret;
			} else{
				double c1 = 1.0;
	                     	double c2 = 0.0;
                        	double c3 = 0.0;
                        	double c4 = 1.0;
                        	double c5 = 1.0;
				double c6 = 1.0;

				if (((a2/R)>0 && fabs((a2/R)-1)>=lim)||((a2/R)<0 && fabs((a2/R)+1)>=lim)){
					c1 = asin(a2/R);
				} else if ((a2/R)>0 && fabs((a2/R)-1)<lim){
					c1 = M_PI/2;
				} else if ((a2/R)<0 && fabs((a2/R)+1)<lim){
					c1 = -M_PI/2;
				}

				if (fabs(1-(a2/R)*(a2/R))>=lim){
					c2 = sqrt(1-(a2/R)*(a2/R));
				} else {
					c2 = 0.0;
				}

			/*	if ((x*x-(a1+d*cos(theta))*(a1+d*cos(theta)))>=0){
					c3 = sqrt(x*x-(a1+d*cos(theta))*(a1+d*cos(theta)));
				} else{
					c3 = 0.0;
				}
                        */
				if ((((a1+d*cos(theta))/x)>0 && fabs(((a1+d*cos(theta))/x)-1)>=lim)||(((a1+d*cos(theta))/x)<0 && fabs(((a1+d*cos(theta))/x)+1)>=lim)){
                                        c4 = asin((a1+d*cos(theta))/x);
                                } else if (((a1+d*cos(theta))/x)>0 && fabs(((a1+d*cos(theta))/x)-1)<lim){
                                        c4 = M_PI/2;
                                } else if (((a1+d*cos(theta))/x)<0 && fabs(((a1+d*cos(theta))/x)+1)<lim){
                                        c4 = -M_PI/2;
                                }

	                /*	if ((x*x-(a2+d*cos(theta))*(a2+d*cos(theta)))>=0){
                                        c5 = sqrt(x*x-(a2+d*cos(theta))*(a2+d*cos(theta)));
                                } else{
                                        c5 = 0.0;
	                	}		
			*/	 
				if ((((a2+d*cos(theta))/x)>0 && fabs(((a2+d*cos(theta))/x)-1)>=lim)||(((a2+d*cos(theta))/x)<0 && fabs(((a2+d*cos(theta))/x)+1)>=lim)){
                                        c6 = asin((a2+d*cos(theta))/x);
                                } else if (((a2+d*cos(theta))/x)>0 && fabs(((a2+d*cos(theta))/x)-1)<lim){
                                        c6 = M_PI/2;
                                } else if (((a2+d*cos(theta))/x)<0 && fabs(((a2+d*cos(theta))/x)+1)<lim){
                                        c6 = -M_PI/2;
                                }

				h_a2 = c1 + (a2/R)*c2;
				f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+c4;
	                        f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+c6; //asin((a2+d*cos(theta))/x);
        	                delta_f = f_a1 - f_a2;
				ret = (M_PI*R*R/2 - (x*x/2)*delta_f - d*sin(fabs(theta))*(a2-a1) - (R*R/2)*(M_PI/2-h_a2));
				
				return ret;
			}
                        break; }

                case 2 :{    //(a-2)
                        double c;
			if ((fabs(a2/R)-1)<lim){
                                c = 1.0;
                        } else {
                                c = a2/R;
                        }
			double h_a2 = asin(c) + (c)*sqrt(1-(c)*(c));
                        double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));

                        return ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));
                        }break;
                case 3 :{   //(a-3)
                        double h_a2 = asin(a2/R) + (a2/R)*sqrt(1-(a2/R)*(a2/R));
                        double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));

                        return ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));
                        }break;
                case 4 :{     //(b) need to swap sign of theta
                        double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(x*x-(a1+d*cos(theta))*(a1+d*cos(theta)))+asin((a1+d*cos(theta))/x);
                        double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(x*x-(a2+d*cos(theta))*(a2+d*cos(theta)))+asin((a2+d*cos(theta))/x);
                        double delta_f = f_a1 - f_a2;
                        double h_a2 = asin(a2/R) + (a2/R)*sqrt(1-(a2/R)*(a2/R));
                        return ((R*R/2)*(h_a2+M_PI/2) + (x*x/2)*delta_f + d*sin(fabs(theta))*(a2-a1));
                        }break;
                case 5 :{    //(c-1)
                        double s = sqrt(x*x-d*d*sin(theta)*sin(theta));
                        double z = d*sin(fabs(theta));

                        return (x*x*acos(z/x) - s*z);
                        }break;
                case 6 :{  //(c-2)
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);

                        return (A_int - M_PI*R*R/2);
                        }break;
                case 7 :{ //(c-3)
                        double s = sqrt(x*x-d*d*sin(theta)*sin(theta));
                        double z = d*sin(fabs(theta));
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);
                        double A1 = x*x*acos(z/x) - s*z;

                        return (A_int - A1);
                        }break;
                case 8 :{ //no intersection
                        return 0;
                        }break;
                case 9 :{ //intersection of two circles
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(w);

                        return A_int;
                        }break;
                case 10 :{ //circle completely inside semi-circle
                        return M_PI*x*x;
                        }break;

        	case 11: { //semi-circle completely inside circle
			return M_PI*R*R/2;
			}break;
		
	        default :
                        return 0;
                        break;
        }

} 


void calc_limb_darkening(double* f_array, double* d_array, int N, double rprs, double fac, int nthreads, double* intensity_args, double phi, double b, double mini)
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
		printf("i = %d ******************** \n",i);
		double d = d_array[i];
		double x = 0.1;
		double x_in = MAX(d - rprs, 0.);					//lower bound for integration
		double x_out = MIN(d + rprs, 1.0);					//upper bound for integration
		if(x_in >= 1.) f_array[i] = 1.0;					//flux = 1. if the planet is not transiting
		else if(x_out - x_in < 1.e-7) f_array[i] = 1.0;				//pathological case	
		else
		{
			double delta = 0.;						//variable to store the integrated intensity, \int I dA
			x = x_in;						//starting radius for integration
			double dx = fac*acos(x); 					//initial step size

			x += dx;						//first step
			double A_i = 0.;						//initial area
		
			double theta = 0.;
			
			/* if (i <= mini){                                     //finding theta for specific values of d, b and phi.
				if ((phi>=0) && (phi<=acos(b/d))) {
					theta = phi + asin(b/d);
				} else if ((phi>acos(b/d)) && (phi<= M_PI/2)) {
					theta = M_PI - asin(b/d) - phi;
				} else if ((phi<=0) && (phi>=-asin(b/d))) {
					theta = asin(b/d) + phi;
				} else if ((phi<-asin(b/d)) && (phi>=-M_PI/2)) {
					theta = -phi -asin(b/d);
				}
			} else {				
				if ((phi >= -M_PI/2) && (phi < -acos(b/d))) {
                                        theta = M_PI+phi - asin(b/d);
                                } else if ((phi<=M_PI/2) && (phi > asin(b/d))) {
                                        theta = - asin(b/d) + phi;
                                } else if ((phi>0) && (phi<=asin(b/d))) {
                                        theta = asin(b/d) - phi;
                                } else if ((phi<=0) && (phi>=-acos(b/d))) {
                                        theta = -phi +asin(b/d);
				}
			} */
		
			if (i<=mini) {	//if planet to left of star
				theta = -phi;
			} else {       // if planet to right of star
				theta = phi; 
			}
			
			while(x < x_out)
			{
				double A_f = area(d, x, rprs, theta);				//calculates area of overlapping circles
				printf("Area = %f\n\n", A_f);
				double I = intensity(x - dx/2., intensity_args); 	//intensity at the midpoint
				delta += (A_f - A_i)*I;				//increase in transit depth for this integration step
				dx = fac*acos(x);  				//updating step size
				x = x + dx;					//stepping to next element
				A_i = A_f;					//storing area
			}
			dx = x_out - x + dx;  					//calculating change in radius for last step  FIXME
			x = x_out;						//final radius for integration
			double A_f = area(d, x, rprs, theta);					//area for last integration step
			printf("Area = %f\n\n", A_f);
			double I = intensity(x - dx/2., intensity_args); 		//intensity at the midpoint
			delta += (A_f - A_i)*I;					//increase in transit depth for this integration step
			f_array[i] = 1.0 - delta;	//flux equals 1 - \int I dA
		}
	}
}
