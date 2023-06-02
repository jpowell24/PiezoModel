// #define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <iostream>

//#include <C:/Users/Jackson/Code/PiezoModel/cpp_code/lapack-3.11/OPENBlas>


// NOTE: I went into the file def_lapack.hpp and commented 
// out a bunch of things that almost certainly destroy the code
// ctrl+f for "//" to find which

using namespace std;
using namespace arma;

int main()
{
    //initialize the random generator 
    //Create a 5 x 5 random matrix and printing it
	
	mat A = randu<mat>(5,5); //random matrix of size 5x5 (generared by this syntax)
	mat B = inv(A);	         //inverse of matrix
    
	cout << "Matrix A::\n"<<endl;
    cout<< A << endl;
	cout << "Matrix B(Inverse)::\n"<<endl;
  	
	cout << B << endl;

	cout << inv(A) << endl;

	return 0;
}