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

const interval delta  = 1e-10;
// A: variable theta of Lemma 6.1 
const interval A = pow(2,1./2)/3*power(delta,interval(3.0));


/////////////////////////
const float deltaf = 1e-10;
const float Cf = pow(3,2./3)*pow(2,-1./3);
const float Af = pow(2,1./2)/3*pow(deltaf,3.0);


// integral J1
interval J1(interval t)
{
	interval fpart = power(t,interval(2./3))*cos(A-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(A-t),interval(-3./2));
	interval spart = interval(3.0)*C*power(t,interval(4./3))*sin(A-t)*sin(A-t)*power(1+power(C,2.0)*power(t,interval(4./3))-2*C*power(t,interval(2./3))*cos(A-t),interval(-5./2));
	
	return fpart - spart;
}


interval part(interval J,int n,int i) // i-th sub interval of J out of n
{
	return J.left()+i*(J.right()-J.left())/n+(J-J.left())/n;
}

// integral J1 from sqrt(2/3)delta^3 to K (capd algorithm)
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



// Computation of J1
interval CAPD_J1()
{

	int n=100000000;
	
	// Upper bound of the derivative
	interval upper_bound= 3*power(K,interval(-1./3))*power(C-power(K,interval(-2./3)),-3)*interval(-1.0,1.0);
	
	// In this case we have to do only from sqrt(2)/3*delta^3 to inf, so we can avoid the lower bound (since it is greater than 0 and the bad angle)
	return I_J1(A,K,n) + upper_bound;

}


int main(int argc, char* argv[])
{
	
	try
	{	

		interval J1 = C*CAPD_J1();
		
		cout << "Value J1 " << J1 << endl;
		
		float I2_aux= pow(2.0/Cf,1./2)*tgamma(2./3)*(sin(Af)/2+pow(3.0,1./2)/2*cos(Af));
		
		//float I2_aux = power(interval(2.0)/C,interval(1./2)) * tgamma(2./3)*(sin(A)/2 + power(interval(3.0),interval(1./2))/2 * cos(A));
		
		interval J2 = interval(I2_aux-pow(deltaf,2.0),I2_aux);
		
		cout << "Value J2 " << J2 << endl;
		
		interval result = J1 - J2;
		
		cout << "The value of the derivative is contained in " << result<< endl;

		
	}
  	catch(exception& e)
    	{
      		cout << "\n\nException caught: " << e.what() << endl;
    	}
    	
    return 0;
}
