#include <iostream>
#include <fstream>
#include "capd/capdlib.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>

using namespace capd;
using namespace std;

// Constant values
const interval C = pow(3,2./3)*pow(2,-1./3);
const float eps = pow(10,-5);

// Angle of discontinuity of the Melnikov function
const interval bad_angle = sqrt(2)/3;

// Bound constants
const interval k = 0.001;
const interval K = 100;


// TOL: distance to the bad_angle, the closest to the bad angle the worst precision we can get in capd
const interval TOL = 0.45;

// num_angles: Number of equidistant angles between 0 and 2pi
const int num_angles = 10000;


// Integral I1.
interval df(interval angles, interval t)
{
	interval fpart = power(t,interval(2./3))*cos(angles-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(angles-t),interval(-3./2));
	interval spart = interval(3.0)*C*power(t,interval(4./3))*sin(angles-t)*sin(angles-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(angles-t),interval(-5./2));
	
	return fpart - spart;
}

interval part(interval J,int n,int i) // i-th sub interval of J out of n
{
	return J.left()+i*(J.right()-J.left())/n+(J-J.left())/n;
}


// Integral from k to K (capd algorithm)
interval I(interval k, interval K, interval angles, int n)
{
	interval sum=0.0;
	
	for(int i=0;i<n;i++)
	{
		interval t=part(intervalHull(k,K),n,i);
		sum = sum + df(angles,t)*K/n;
	}
	return sum;
}

// Computation of I1
interval CAPD(interval angles, float dist)
{
	// Number of subintervals for each iteration of the algorithm. The closer to the singularity, the better has to be the precision
	// and therefore more subintervals are needed.
	int n=1000000/dist; 

	// A toy example to check the execution
	// int n = 100000;

	// The upper and lower bound for the derivative 
	interval dlower_bound = power(k,interval(5./3))*interval(3./5)* power(1+power(C,interval(2.0))*power(k,interval(4./3)) -2*C*power(k,interval(2./3)),interval(-3./2))*interval(-1.0,1.0);
	interval dupper_bound = interval(3.0)*power(K,interval(-1./3))*power(C-power(K,interval(-2./3)),-3)*interval(-1.0,1.0);
	
	return dlower_bound + I(k,K,angles,n) + dupper_bound;

}

int main(int argc, char* argv[])
{
	try
	{		
		// Write data in file
		ofstream file;
		ofstream file_angle;
		
		file.open("data_der.txt");
		file_angle.open("angles_der.txt");
		
		
		for (int i = 0; i <= num_angles; i++)
		{	
			
			
			// Interval to compute the integral. The constant "eps" is written to ensure the adjacent intervals intersect.
			interval angles = interval(atan(1)*8*i/num_angles, atan(1)*8*(i+1)/num_angles + eps);
		
			if (abs(angles-bad_angle) > TOL){
				float dist = abs(atan(1)*8*i/num_angles-sqrt(2)/3);
				file_angle << angles << endl;
				
				// I1 is obtained through the CAPD algorithm.
				interval I1 = C*CAPD(angles,dist);
				
				// I2 is obtained explicitly.
				interval I2 = power(interval(2.0)/C,interval(1./2)) * sin(angles)*interval(tgamma(2./3)/2) - power(interval(2.0)/C,interval(1./2)) * cos(angles)*interval(4*atan(1)/tgamma(1./3));
				
				interval result = I1 - I2;
				file << result << endl;
				
			}
		}
		file_angle.close();
		file.close();

	}
  	catch(exception& e)
    	{
      		cout << "\n\nException caught: " << e.what() << endl;
    	}
    	
    return 0;
}
