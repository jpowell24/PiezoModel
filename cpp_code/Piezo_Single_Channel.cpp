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


normal_distribution<double> stiffness(0.5,0.01);
normal_distribution<double> pressure(30,0.01);
normal_distribution<double> pressure2(0,5);
normal_distribution<double> voltage(-70,1);

vector<double> vec_P_Substrate;
vector<double> vec_P_Pressure;
vector<double> vec_P_Voltage;
vector<double> vec_P_Total; 

vector<int> vec_Channel;

vector<int> vec_Channel0;
vector<int> vec_Channel1;
vector<int> vec_Channel2;
vector<int> vec_Channel3;
vector<int> vec_Channel4;
vector<int> vec_Channel5;
vector<int> vec_Channel6;
vector<int> vec_Channel7;
vector<int> vec_Channel8;
vector<int> vec_Channel9;

double delta_T = 1; 

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

double P_open(double x){
    double P_open; 
    double Pressure_input = 30; 
    double Substrate_input = 0.5; 
    double Voltage_input = -70;

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);

    double P_total = P_P*P_S + P_V;

    vec_P_Total.push_back(P_total);

    P_open = 1/(exp((0.5 - P_total)/0.1) + 1) - 0.00669;
    
    return(P_open);
}

double P_close(double x){
    double P_close;

    double tau_inact = exp(-delta_T/500);
    double tau_open = exp(-delta_T/12); 

    P_close = tau_open;

    return(P_close);
}

double Piezo_Channel(int x){
    double p;
    vec_Channel.push_back(0);

    for(int time = 0; time < 100; time++){

        if(vec_Channel[time] == 0){
            p = P_open(0);
        }
        else if(vec_Channel[time] == 1){
            p = P_close(0);
        }
        else{
            cout << "You have an error" << endl;
        }

        random_device rd;
        mt19937 gen(rd());
        discrete_distribution<> distrib({1 - p, p});

        if(x == 0){
            vec_Channel0.push_back(distrib(gen));
        }
        else if(x == 1){
            vec_Channel1.push_back(distrib(gen));
        }
        else if(x == 2){
            vec_Channel2.push_back(distrib(gen));
        }
        else if(x == 3){
            vec_Channel3.push_back(distrib(gen));
        }
        else if(x == 4){
            vec_Channel4.push_back(distrib(gen));
        }
        else if(x == 5){
            vec_Channel5.push_back(distrib(gen));
        }
        else if(x == 6){
            vec_Channel6.push_back(distrib(gen));
        }
        else if(x == 7){
            vec_Channel7.push_back(distrib(gen));
        }
        else if(x == 8){
            vec_Channel8.push_back(distrib(gen));
        }
        else if(x == 9){
            vec_Channel9.push_back(distrib(gen));
        }
        else{
            cout << "There is an error" << endl;
        }


    }

    return(0);
}

double output_file(double x){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    for(int i = 0; i < 10; i++){
        Piezo_Channel(i);
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_Channel.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_Open;

    myfile << "Piezo_Open\n";

    for (int i = 0; i < max_size - 1; i++)
    {
        Bool_Open = (vec_Channel.size() > i) ? true : false;

        if(Bool_Open) myfile << vec_Channel[i];

        myfile << "\n";
    }

    myfile.close();

    return(0);
}


int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}