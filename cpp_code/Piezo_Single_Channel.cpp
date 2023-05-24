#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/Piezo_Single_Channel.csv";

default_random_engine generator;
// normal_distribution<double> stochastic_opening(0,4);


normal_distribution<double> stiffness(0.2,0.01);
normal_distribution<double> pressure(0,0.01);
normal_distribution<double> pressure2(0,5);
normal_distribution<double> voltage(-70,1);

vector<double> vec_P_Substrate;
vector<double> vec_P_Pressure;
vector<double> vec_P_Voltage;
vector<double> vec_P_Total; 

double Reset_vecs(double i){

    return(0);
}

double Piezo_P_Pressure(double i){
    double F_inf;
    F_inf = 1/(exp((30 - i)/7) + 1);
    vec_P_Pressure.push_back(F_inf);
    return(F_inf);
}

double Piezo_P_Substrate(double i){
    double S_inf;
    S_inf = (1/(0.25*pow(2*M_PI,0.5))*exp(-0.5*pow((i - 0.7)/0.25,2)))/1.6;
    vec_P_Substrate.push_back(S_inf);
    return(S_inf);
}

double Piezo_P_Voltage(double i){
    double V_inf;
    V_inf = 1/(exp((100 - i)/20) + 1);
    vec_P_Voltage.push_back(V_inf);
    return(V_inf);
}

double Piezo_Channel(int time, double pressure_temp){

    return(0);
}


double output_file(double x){
    double p = Piezo_P_Pressure(30);

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> distrib({1 - p, p});

    for(int i = 0; i <= 100; i++){
        cout << distrib(gen) << endl;
        x += distrib(gen);
    }

    cout << x << endl;

    return(0);
}


int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}