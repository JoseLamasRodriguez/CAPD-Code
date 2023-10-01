#include <iostream>
#include <fstream>
#include "capd/capdlib.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>

using namespace capd;
using namespace std;

// CONSTANT VALUES

const interval C = pow(3,2./3)*pow(2,-1./3);
const interval bad_angle = sqrt(2)/3;

// Value of K in J1 and J2 
const interval K = 10000.0;

// A: variable theta of Lemma 6.1 
const interval A = interval(4./3);


// Integral J1
interval J1(interval t)
{
	interval fpart = power(t,interval(2./3))*cos(-A-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(-A-t),interval(-3./2));
	interval spart = interval(3.0)*C*power(t,interval(4./3))*sin(-A-t)*sin(-A-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(-A-t),interval(-5./2));
	
	return fpart - spart;
}

// Integral J2
interval J2(interval t)
{

	interval integrand = cos(A+t)*power(t,interval(-4./3));
	return integrand;
}


interval part(interval J,int n,int i) // i-th sub interval of J out of n
{
	return J.left()+i*(J.right()-J.left())/n+(J-J.left())/n;
}

// Integral J1 from 4/3 to K (capd algorithm)
interval I_J1(interval k, interval K, int n)
{
	interval sum=0.0;
	
	for(int i=0;i<n;i++)
	{
		interval t=part(intervalHull(k,K),n,i);
		sum = sum + J1(t)*K/n;
		
	}
	return sum;
}

// Integral J2 from 4/3 to K (capd algorithm)
interval I_J2(interval k, interval K, int n)
{
	interval sum=0.0;
	
	for(int i=0;i<n;i++)
	{
		interval t=part(intervalHull(k,K),n,i);
		sum = sum + J2(t)*K/n;
		
	}
	return sum;
}

// Computation of J1
interval CAPD_J1()
{

	int n=100000000;
	
	// Upper bound of the derivative
	interval upper_bound= 3*power(K,interval(-1./3))*power(C-power(K,interval(-2./3)),-3)*interval(-1.0,1.0);
	
	// In this case we have to do only from 4/3 to inf, so we can avoid the lower bound (since it is greater than 0 and the bad angle)
	return I_J1(A,K,n) + upper_bound;

}


// Computation of J2
interval CAPD_J2()
{

	int n= 100000000;
	
	// Upper bound of the integral
	interval upper_bound = 3*power(K,interval(-1./3));
	
	cout << "Upper bound J2 " << upper_bound << endl;
	return I_J2(A,K,n) + upper_bound;
}
int main(int argc, char* argv[])
{
	
	try
	{	
		
		interval J1 = C*CAPD_J1();
		interval J2 = power(interval(2.0)/C,interval(1./2)) * (cos(8./3)/power(interval(4./3),interval(1./3)) - 1./3 * CAPD_J2());
		
		interval result = J1 + J2;
		
		cout << "The value of the derivative is contained in " << result<< endl;

		
	}
  	catch(exception& e)
    	{
      		cout << "\n\nException caught: " << e.what() << endl;
    	}
    	
    return 0;
}
