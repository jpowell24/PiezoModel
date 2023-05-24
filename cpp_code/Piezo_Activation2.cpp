#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/Piezo_Activation.csv";

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

vector<double> vec_open1;
vector<double> vec_open2;
vector<double> vec_inactive;
vector<double> vec_inactive_held;
vector<double> vec_inactive2;
vector<double> vec_open3;
vector<double> vec_closed;
vector<double> vec_closed2;
vector<double> vec_current;

vector<double> vec_current10;
vector<double> vec_current20;
vector<double> vec_current30;
vector<double> vec_current40;
vector<double> vec_current50;
vector<double> vec_current60;

vector<double> vec_P_total10;
vector<double> vec_P_total20;
vector<double> vec_P_total30;
vector<double> vec_P_total40;
vector<double> vec_P_total50;
vector<double> vec_P_total60;

double open_local, inactive, closed, tau_inact, tau_open3, tau_open, tau_inact2, tau_open2; 

double local_N_Piezo = 100;

double Reset_vecs(double i){

    vec_open1.clear();
    vec_open2.clear();
    vec_inactive.clear();
    vec_inactive2.clear();
    vec_inactive_held.clear();
    vec_open3.clear();
    vec_closed.clear();
    vec_closed2.clear();
    vec_P_Substrate.clear();
    vec_P_Pressure.clear();
    vec_P_Voltage.clear();
    vec_P_Total.clear(); 

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

    double stochastic_stiffness = abs(stiffness(generator));
    double stochastic_pressure = abs(pressure(generator));
    double stochastic_voltage = voltage(generator);

    //cout << stochastic_stiffness << " " << stochastic_pressure << " " << stochastic_voltage << endl;

    double P_opening_temp;
    double Pressure_input, Substrate_input, Voltage_input; 

    Pressure_input = stochastic_pressure;
    Substrate_input = stochastic_stiffness;
    Voltage_input = stochastic_voltage;

    if(time > 100 && time < 600){
        Pressure_input = pressure_temp + ((pressure_temp/100)*pressure2(generator));
        Pressure_input = abs(Pressure_input);
    }

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);

    double P_total = P_P*P_S + P_V;

    vec_P_Total.push_back(P_total);

    P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1) - 0.00669;

    tau_inact = 0.999;
    tau_open = 0.9; 

    double delta_T = 1;

    vec_open1.push_back(tau_open*vec_open1[time] + (1 - exp(-delta_T*P_opening_temp))*vec_closed[time]);
    vec_inactive_held.push_back(((1 - exp(-delta_T*P_P))*vec_inactive_held[time]) + (P_P*vec_open1[time])*(1-tau_open));
    vec_inactive.push_back(tau_inact*vec_inactive[time] + ((1 - P_P)*(1-tau_open)*vec_open1[time]) + ((exp(-delta_T*P_P))*vec_inactive_held[time]));
    vec_closed.push_back((exp(-delta_T*P_opening_temp))*vec_closed[time] + (vec_inactive[time] - tau_inact*vec_inactive[time]));

    tau_inact2 = 0.99;
    tau_open2 = 0.98; 

    // vec_open2.push_back(tau_open2*vec_open2[time] + P_P*vec_closed2[time]);
    // vec_inactive2.push_back(tau_inact2*vec_inactive2[time] + (1 - tau_open2)*vec_open2[time]);
    // vec_closed2.push_back((1 - P_P)*vec_closed2[time] + (1 - tau_inact2)*vec_inactive2[time]);
    
    if(time < 599){
        vec_open2.push_back(P_P*vec_closed2[0]);
    }
    else{
        vec_open2.push_back(tau_open2*vec_open2[time]);
    }
    

    //cout << inverse_P_total << " x " << inverse_P_total << "  = " << inverse_P_total*vec_closed2[time] << endl;
    //cout << vec_closed2[time] << " and " << vec_open2[time] << endl;

    // vec_open2[time] = 0;

    // 0.95 = e^-1/10

    //open_local = tau_open*vec_open1[time] + P_opening_temp*vec_closed[time];
    // cout << open_local << " and " << P_opening_temp*vec_closed[time] << endl;

    // cout << P_opening_temp << endl;
    // open_local = vec_num_open[time] + P_opening_temp*closed;

    // if(open_local >= local_N_Piezo){
    //     open_local = local_N_Piezo;
    // }

    if(pressure_temp == 10){
        vec_current10.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total10.push_back(-0.25*vec_open2[time]);
    }
    else if(pressure_temp == 20){
        vec_current20.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total20.push_back(-0.25*vec_open2[time]);
    }
    else if(pressure_temp == 30){
        vec_current30.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total30.push_back(-0.25*vec_open2[time]);
    }
    else if(pressure_temp == 40){
        vec_current40.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total40.push_back(-0.25*vec_open2[time]);
    }
    else if(pressure_temp == 50){
        vec_current50.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total50.push_back(-0.25*vec_open2[time]);
    }
    else{
        vec_current60.push_back(-15*vec_open1[time] + -2*vec_inactive_held[time] + -0.25*vec_open2[time]);
        vec_P_total60.push_back(-0.25*vec_open2[time]);
    }

    return(0);
}

double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    for(double x = 10; x <= 60; x+=10){
        vec_open1.push_back(0);
        vec_open2.push_back(0);
        vec_inactive.push_back(0);
        vec_inactive2.push_back(0);
        vec_inactive_held.push_back(0);
        vec_open3.push_back(0);
        vec_closed.push_back(1);
        vec_closed2.push_back(1);

        for(double i = 0; i <= 800; i+= 1){
            Piezo_Channel(i, x);
        }

        if(x < 60){
            Reset_vecs(0);
        }
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_open1.size());
    sizes.insert(sizes.begin(),vec_open2.size());
    sizes.insert(sizes.begin(),vec_P_Pressure.size());
    sizes.insert(sizes.begin(),vec_P_Substrate.size());
    sizes.insert(sizes.begin(),vec_P_Voltage.size());
    sizes.insert(sizes.begin(),vec_closed.size());
    sizes.insert(sizes.begin(),vec_inactive.size());
    sizes.insert(sizes.begin(),vec_current.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open1;
    bool Bool_P_open2;
    bool Bool_P_open3;
    bool Bool_open1;
    bool Bool_open2;
    bool Bool_open3;
    bool Bool_P_active;
    bool Bool_P_inactive;
    bool Bool_current;

    bool Bool_current10;
    bool Bool_current20;
    bool Bool_current30;
    bool Bool_current40;
    bool Bool_current50;
    bool Bool_current60;

    bool Bool_P_total10;
    bool Bool_P_total20;
    bool Bool_P_total30;
    bool Bool_P_total40;
    bool Bool_P_total50;
    bool Bool_P_total60;


    //cout << "Break point 4" << endl;

    myfile << "Piezo_Open1,Piezo_Inactive1,Piezo_Inactive_held,Piezo_Closed1,Piezo_Open2,Piezo_Closed2,Substrate,Voltage,Current,10,20,30,40,50,60,P10,P20,P30,P40,P50,P60\n";

    for (int i = 0; i < max_size - 1; i++)
    {
        //cout << "Break point 5" << endl;
        Bool_open1 = (vec_open1.size() > i) ? true : false;
        Bool_open2 = (vec_open2.size() > i) ? true : false;
        Bool_open3 = (vec_open3.size() > i) ? true : false;
        Bool_P_open1 = (vec_P_Pressure.size() > i) ? true : false;
        Bool_P_open2 = (vec_P_Substrate.size() > i) ? true : false;
        Bool_P_open3 = (vec_P_Voltage.size() > i) ? true : false;
        Bool_P_active = (vec_closed.size() > i) ? true : false;
        Bool_P_inactive = (vec_inactive.size() > i) ? true : false;
        Bool_current = (vec_current.size() > i) ? true : false;

        Bool_current10 = (vec_current10.size() > i) ? true : false;
        Bool_current20 = (vec_current20.size() > i) ? true : false;
        Bool_current30 = (vec_current30.size() > i) ? true : false;
        Bool_current40 = (vec_current40.size() > i) ? true : false;
        Bool_current50 = (vec_current50.size() > i) ? true : false;
        Bool_current60 = (vec_current60.size() > i) ? true : false;

        Bool_P_total10 = (vec_P_total10.size() > i) ? true : false;
        Bool_P_total20 = (vec_P_total20.size() > i) ? true : false;
        Bool_P_total30 = (vec_P_total30.size() > i) ? true : false;
        Bool_P_total40 = (vec_P_total40.size() > i) ? true : false;
        Bool_P_total50 = (vec_P_total50.size() > i) ? true : false;
        Bool_P_total60 = (vec_P_total60.size() > i) ? true : false;


        //cout << "Break point 6" << endl;

        if(Bool_open1) myfile << vec_open1[i] << ",";
        if(!Bool_open1) myfile << ",";
        if(Bool_open1) myfile << vec_inactive[i] << ",";
        if(!Bool_open1) myfile << ",";
        if(Bool_P_inactive) myfile << vec_inactive_held[i] << ",";
        if(!Bool_P_inactive) myfile << ",";
        if(Bool_P_active) myfile << vec_closed[i] << ",";
        if(!Bool_P_active) myfile <<",";
        if(Bool_open2) myfile << vec_open2[i] << ",";
        if(!Bool_open2) myfile << ",";
        if(Bool_P_open1) myfile << vec_closed2[i] << ",";
        if(!Bool_P_open1) myfile << ",";
        if(Bool_P_open2) myfile << vec_P_Substrate[i] << ",";
        if(!Bool_P_open2) myfile << ",";
        if(Bool_P_open3) myfile << vec_P_Voltage[i] << ",";
        if(!Bool_P_open3) myfile << ",";
        if(Bool_current) myfile << vec_current[i]  << ",";
        if(!Bool_current) myfile << ",";

        if(Bool_current10) myfile << vec_current10[i]  << ",";
        if(!Bool_current10) myfile << ",";
        if(Bool_current20) myfile << vec_current20[i]  << ",";
        if(!Bool_current20) myfile << ",";
        if(Bool_current30) myfile << vec_current30[i]  << ",";
        if(!Bool_current30) myfile << ",";
        if(Bool_current40) myfile << vec_current40[i]  << ",";
        if(!Bool_current40) myfile << ",";
        if(Bool_current50) myfile << vec_current50[i]  << ",";
        if(!Bool_current50) myfile << ",";
        if(Bool_current60) myfile << vec_current60[i] << ",";
        if(!Bool_current60) myfile << ",";

        if(Bool_P_total10) myfile << vec_P_total10[i]  << ",";
        if(!Bool_P_total10) myfile << ",";
        if(Bool_P_total20) myfile << vec_P_total20[i]  << ",";
        if(!Bool_P_total20) myfile << ",";
        if(Bool_P_total30) myfile << vec_P_total30[i]  << ",";
        if(!Bool_P_total30) myfile << ",";
        if(Bool_P_total40) myfile << vec_P_total40[i]  << ",";
        if(!Bool_P_total40) myfile << ",";
        if(Bool_P_total50) myfile << vec_P_total50[i]  << ",";
        if(!Bool_P_total50) myfile << ",";
        if(Bool_P_total60) myfile << vec_P_total60[i];


        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }

    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}