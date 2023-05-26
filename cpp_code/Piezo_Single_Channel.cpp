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

vector<double> vec_open1;
vector<double> vec_inactive;
vector<double> vec_inactive_held;
vector<double> vec_closed;

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

    vec_P_Substrate.clear();
    vec_P_Pressure.clear();
    vec_P_Voltage.clear();
    vec_P_Total.clear();
    vec_Channel.clear();

    vec_open1.clear();
    vec_inactive.clear();
    vec_inactive_held.clear();
    vec_closed.clear();

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
    double Pressure_input = x; 
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

double Piezo_Channel(int x, double pressure){
    double p;

    cout << pressure << endl;

    vec_Channel.push_back(0);
    vec_open1.push_back(0);
    vec_inactive_held.push_back(0);
    vec_inactive.push_back(0);
    vec_closed.push_back(1);

    random_device rd;
    mt19937 gen(rd());

    //cout << "I have been called " << x << " times" << endl;
    int time = 0; 

    for(double counter = 0; counter < 800; counter+=delta_T){

        double Pressure_input = pressure; 
        double Substrate_input = 0.5; 
        double Voltage_input = -70; 

        if(counter < 200 || counter > 600){
            Pressure_input = 0;
        }
        else{
            Pressure_input = pressure;
        }
        //cout << "I have been called " << x << " times" << endl;

        // if(vec_Channel[time] == 0){
        //     p = P_open(pressure);
        // }
        // else if(vec_Channel[time] == 1){
        //     p = P_close(pressure);
        // }
        // else{
        //     cout << "You have an error" << endl;
        // }

        double P_P = Piezo_P_Pressure(Pressure_input);
        double P_S = Piezo_P_Substrate(Substrate_input);
        double P_V = Piezo_P_Voltage(Voltage_input);

        double P_total = P_P*P_S + P_V;

        vec_P_Total.push_back(P_total);

        double P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1) - 0.00669;

        double tau_inact = 0.999;
        double tau_open = 0.9; 
        
        vec_open1.push_back(tau_open*vec_open1[time] + (1 - exp(-delta_T*P_opening_temp))*vec_closed[time]);
        vec_inactive_held.push_back(((1 - exp(-delta_T*P_P))*vec_inactive_held[time]) + (P_P*vec_open1[time])*(1-tau_open));
        vec_inactive.push_back(tau_inact*vec_inactive[time] + ((1 - P_P)*(1-tau_open)*vec_open1[time]) + ((exp(-delta_T*P_P))*vec_inactive_held[time]));
        vec_closed.push_back((exp(-delta_T*P_opening_temp))*vec_closed[time] + (vec_inactive[time] - tau_inact*vec_inactive[time]));

        p = vec_open1[time];

        discrete_distribution<> distrib({1 - p, p});

        int temp = distrib(gen);

        if(x == 0){
            vec_Channel0.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 1){
            vec_Channel1.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 2){
            vec_Channel2.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 3){
            vec_Channel3.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 4){
            vec_Channel4.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 5){
            vec_Channel5.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 6){
            vec_Channel6.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 7){
            vec_Channel7.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 8){
            vec_Channel8.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else if(x == 9){
            vec_Channel9.push_back(temp);
            vec_Channel.push_back(temp);
        }
        else{
            cout << "There is an error" << endl;
        }
        
        time++;
    }

    return(0);
}

double output_file(double x){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vec_Channel0.push_back(0);
    vec_Channel1.push_back(0);
    vec_Channel2.push_back(0);
    vec_Channel3.push_back(0);
    vec_Channel4.push_back(0);
    vec_Channel5.push_back(0);
    vec_Channel6.push_back(0);
    vec_Channel7.push_back(0);
    vec_Channel8.push_back(0);
    vec_Channel9.push_back(0);

    for(int i = 0; i < 10; i++){
        Piezo_Channel(i, (10+5*i));
        Reset_vecs(0);
        // cout << i << endl;
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_Channel0.size());
    sizes.insert(sizes.begin(),vec_Channel1.size());
    sizes.insert(sizes.begin(),vec_Channel2.size());
    sizes.insert(sizes.begin(),vec_Channel3.size());
    sizes.insert(sizes.begin(),vec_Channel4.size());
    sizes.insert(sizes.begin(),vec_Channel5.size());
    sizes.insert(sizes.begin(),vec_Channel6.size());
    sizes.insert(sizes.begin(),vec_Channel7.size());
    sizes.insert(sizes.begin(),vec_Channel8.size());
    sizes.insert(sizes.begin(),vec_Channel9.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_Open0;
    bool Bool_Open1;
    bool Bool_Open2;
    bool Bool_Open3;
    bool Bool_Open4;
    bool Bool_Open5;
    bool Bool_Open6;
    bool Bool_Open7;
    bool Bool_Open8;
    bool Bool_Open9;


    myfile << "Piezo_Open0,Piezo_Open1,Piezo_Open2,Piezo_Open3,Piezo_Open4,Piezo_Open5,Piezo_Open6,Piezo_Open7,Piezo_open8,Piezo_open9\n";

    for (int i = 0; i < max_size - 1; i++)
    {
        Bool_Open0 = (vec_Channel0.size() > i) ? true : false;
        Bool_Open1 = (vec_Channel1.size() > i) ? true : false;
        Bool_Open2 = (vec_Channel2.size() > i) ? true : false;
        Bool_Open3 = (vec_Channel3.size() > i) ? true : false;
        Bool_Open4 = (vec_Channel4.size() > i) ? true : false;
        Bool_Open5 = (vec_Channel5.size() > i) ? true : false;
        Bool_Open6 = (vec_Channel6.size() > i) ? true : false;
        Bool_Open7 = (vec_Channel7.size() > i) ? true : false;
        Bool_Open8 = (vec_Channel8.size() > i) ? true : false;
        Bool_Open9 = (vec_Channel9.size() > i) ? true : false;

        if(Bool_Open0) myfile << vec_Channel0[i]  << ",";
        if(!Bool_Open0) myfile << ",";
        if(Bool_Open1) myfile << vec_Channel1[i]  << ",";
        if(!Bool_Open1) myfile << ",";
        if(Bool_Open2) myfile << vec_Channel2[i]  << ",";
        if(!Bool_Open2) myfile << ",";
        if(Bool_Open3) myfile << vec_Channel3[i]  << ",";
        if(!Bool_Open3) myfile << ",";
        if(Bool_Open4) myfile << vec_Channel4[i]  << ",";
        if(!Bool_Open4) myfile << ",";
        if(Bool_Open5) myfile << vec_Channel5[i]  << ",";
        if(!Bool_Open5) myfile << ",";
        if(Bool_Open6) myfile << vec_Channel6[i]  << ",";
        if(!Bool_Open6) myfile << ",";
        if(Bool_Open7) myfile << vec_Channel7[i]  << ",";
        if(!Bool_Open7) myfile << ",";
        if(Bool_Open8) myfile << vec_Channel8[i]  << ",";
        if(!Bool_Open8) myfile << ",";
        if(Bool_Open9) myfile << vec_Channel9[i];

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