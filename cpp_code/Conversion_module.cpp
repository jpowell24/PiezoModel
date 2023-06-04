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

VectorXd vec_std_to_Eigen(vector<double> convert_me){
    int size = convert_me.size(); 
    VectorXd A(size);
    for(int i = 0; i < size; i++){
        A[i] = convert_me[i];
    }
    return(A);
}

vector<double> vec_Eigen_to_std(VectorXd convert_me){ 
    int size = convert_me.size();
    vector<double> A(size);
    for(int i = 0; i < size; i++){
        A[i] = convert_me[i];
    }
    return(A);
}

int generator(int x){
    vector<double> test;
    VectorXd test2(3);

    for(int i = 0; i < 3; i++){
        test2[i] = i;
    }

    test = vec_Eigen_to_std(test2);

    for(int i = 0; i < test.size(); i++){
        cout << test[i] << endl;
    }

    return(0);
}

int main(){
	cout << "begin" << endl;
    generator(0);
	cout << "end" << endl;
	return(0);
}