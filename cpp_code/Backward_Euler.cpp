#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>


using namespace std;
using namespace Eigen;

const char *path1="../data_files/Diffusion_test.csv";


int main(){
	cout << "begin" << endl;
	ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

	double delta_t = 0.01/1000; //in seconds
	double delta_x = 0.1/10; //in meters

	double total_t = 0.001;
	double time_index = total_t/delta_t;

	double initial_C = 0.1;

	int size = 100; 

	
	MatrixXd A = MatrixXd::Zero(size, size);
	
	// cout << "Here is the matrix A:\n" << A << "\n" << endl;

	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			if((i < size - 1) && (j < size - 1) && (i > 0) && (j > 0)){
				if(i == j){
					A(i,j) = -2;
					A(i,j + 1) = 1;
					A(i,j - 1) = 1;
				}
			}
			if((i == 0)){
				if(i == j){
					A(i,j) = -1;
					A(i,j + 1) = 1;
				}
			}
			if((i == size - 1)){
				if(i == j){
					A(i,j) = -1;
					A(i,j - 1) = 1;
				}
			}
		}
	}

	A = (delta_t / pow(delta_x,2))*A;

	// cout << A << endl;

	VectorXd F1(size);
	for(int i = 0; i < size; i++){
		F1[i] = initial_C;
	}
	F1[50] = 10;


	VectorXd F2 = VectorXd::Zero(size);

	// cout << F1 << " and " << endl;

	// cout << F1 << "\n" << endl;

	MatrixXd I(size, size);
	I.setIdentity(size, size);

	for(double time = 0; time < time_index; time ++){

		F2 = (I - A).inverse()*F1;

		// cout << F2 << endl;

		myfile << "Time " << time << ",";

		for (int i = 0; i < size - 1; i++)
    	{
        	if(time < time_index - 1) myfile << F1[i] << ",";
        	if(time == time_index - 1) myfile << F1[i] << ",";
    	}

		myfile << "\n";

		F1 = F2; 

		cout << time << endl;
	}

	myfile.close();
	cout << "end" << endl;

	return(0);
}