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

double delta_t = 0.1; 
double delta_x = 0.1; 

bool rows_eigen_A_initialized; 
bool cols_eigen_A_initialized; 
bool eigen_A_initialized; 

// MatrixXd euler_A;
MatrixXd rows_euler_A;
MatrixXd cols_euler_A;



const char *path1="../data_files/Diffusion_test.csv";

MatrixXd Pairwise_distances(MatrixXd A, MatrixXd B){
    MatrixXd Distances(A.cols(),B.cols());
    for(int i = 0; i < A.cols(); i++){
        for(int j = 0; j < B.cols(); j++){
            Distances(i,j) = sqrt((A.col(i) - B.col(j)).array().pow(2).sum());
        }
    }
    return(Distances);
}

VectorXd vec_std_to_Eigen(vector<double> convert_me){
    int size = convert_me.size(); 
    VectorXd A(size);
    for(int i = 0; i < size; i++){
        A[i] = convert_me[i];
    }
    return(A);
}

MatrixXd mat_std_to_Eigen(vector<vector<double>> convert_me){
    int rows = convert_me.size(); 
    int cols = convert_me[0].size(); 
    MatrixXd A(rows,cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            A(i,j) = convert_me[i][j];
        }
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

vector<vector<double>> mat_Eigen_to_std(MatrixXd convert_me){
    int rows = convert_me.rows(); 
    int cols = convert_me.cols(); 
    vector<double> B;
    vector<vector<double>> A;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            B.push_back(convert_me(i,j));
        }
        A.push_back(B);
        B.clear();
    }
    return(A);
}

MatrixXd euler_A_maker(int size){

    MatrixXd euler_A_temp = MatrixXd::Zero(size, size);

        for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if((i < size - 1) && (j < size - 1) && (i > 0) && (j > 0)){
                if(i == j){
                    euler_A_temp(i,j) = -2;
                    euler_A_temp(i,j + 1) = 1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
            if((i == 0)){
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j + 1) = 1;
                }
            }
            if((i == size - 1)){
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
        }
    }
    eigen_A_initialized = true; 

    return(euler_A_temp);
}

VectorXd backward_euler(VectorXd Diffuse_me, MatrixXd euler_A){
    int size = Diffuse_me.size();

    euler_A = (delta_t / pow(delta_x,2))*euler_A;
    VectorXd F2 = VectorXd::Zero(size);

    MatrixXd I(size, size);
    I.setIdentity(size, size);

    F2 = (I - euler_A).inverse()*Diffuse_me;

    return(F2);
}

MatrixXd call_backward_euler(MatrixXd t1){


    MatrixXd thalfway(t1.rows(),t1.cols());
    MatrixXd t2(t1.rows(),t1.cols()); 

    MatrixXd A_cols; 
    A_cols = euler_A_maker(t1.rows());

    MatrixXd A_rows; 
    A_rows = euler_A_maker(t1.cols());
    
    for(int i = 0; i < t1.cols(); i++){
        thalfway.col(i) = backward_euler(t1.col(i), A_cols);
    }

    for(int i = 0; i < thalfway.rows(); i++){
        t2.row(i) = backward_euler(thalfway.row(i), A_rows);
    }

    return(t2);
}

int generator(int x){
    MatrixXd A(20,15);
    A <<    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,100,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;

    cout << "here" << endl;

    A = call_backward_euler(A);

    cout << A << endl;

    return(0);
}

int main(){
	cout << "begin" << endl;

    generator(0);

	cout << "end" << endl;
}