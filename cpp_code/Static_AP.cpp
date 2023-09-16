#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

#include "Calcium_Dynamics.h"

using namespace std;

const char *path1="../data_files/Piezo_Channel.csv";
const char *path2="../data_files/static_wt_output.csv";


default_random_engine generator;
normal_distribution<double> stochastic_opening(0,4);

int reset_vecs(int x){
    vec_n.clear();
    vec_m.clear();
    vec_h.clear();
    vec_inf_n.clear();
    vec_inf_m.clear();
    vec_inf_h.clear();
    vec_Na_I.clear();
    vec_K_I.clear();
    vec_L_I.clear();

    return(0);
}

double dynamical_h(double V){
    //for some reason I decided to add in double h, so that we can store this value local to the method
    //and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    a_h = 0.07*exp(-V/20);
    b_h = ((1)/(exp((30-V)/10) + 1));
    h_inf = a_h/(a_h + b_h);
    tau_h = 1/(a_h + b_h);
    //the if statements are added to expunge NaN from the data
    //storing the values in vectors so that they can be easily written to a csv file
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);
    return(0);
}

double dynamical_n(double V){
    a_n = 0.01*((10-V)/(exp((10-V)/10) - 1));
    if(V == 10){
        a_n = 0.1; // This is the Taylor approx value for when divide by 0
    }
    b_n = 0.125*exp(-V/80);
    n_inf = a_n/(a_n + b_n);
    tau_n = 1/(a_n + b_n);
    //cout << tau_n << endl;
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    //cout << V << endl;
    return (0);
}

double dynamical_m(double V){
    a_m = 0.1*((25 - V)/(exp((25-V)/10) - 1));
    if (V == 25){
        a_m = 1; // This is the Taylor approx value for when divide by 0
    }
    b_m = 4*exp(-V/18);
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    //cout << tau_m << endl;
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return (0);
}


double Static_WT_AP(double local_g_k){

        ofstream create_file(path1);
        ofstream myfile;
        myfile.open(path1);
    
        // if(static_ap_counter >= 5000){
        //   current = 5;
        // }

        //cout << "Break point 1" << endl;

        double current = 0; 

        dynamical_m(vec_V[static_ap_counter]);
        dynamical_h(vec_V[static_ap_counter]);
        dynamical_n(vec_V[static_ap_counter]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[static_ap_counter] + delta_T*((vec_inf_n[static_ap_counter] - vec_n[static_ap_counter])/vec_tau_n[static_ap_counter]));
        vec_m.push_back(vec_m[static_ap_counter] + delta_T*((vec_inf_m[static_ap_counter] - vec_m[static_ap_counter])/vec_tau_m[static_ap_counter]));
        vec_h.push_back(vec_h[static_ap_counter] + delta_T*((vec_inf_h[static_ap_counter] - vec_h[static_ap_counter])/vec_tau_h[static_ap_counter]));

        //cout << "Break point 3" << endl;

        K_I_temp = (local_g_k*pow(vec_n[static_ap_counter+1],4)*((vec_V[static_ap_counter]) - E_k));
        Na_I_temp = (g_Na*pow(vec_m[static_ap_counter+1],3)*pow(vec_h[static_ap_counter+1],1)*((vec_V[static_ap_counter]) - E_Na));
        L_I_temp = (g_l*((vec_V[static_ap_counter]) - E_l));

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp)/C_m;

        if(local_g_k > 20){
            vec_V.push_back(vec_V[static_ap_counter] + delta_T*V_dt);
        }
        else{
            vec_V2.push_back(vec_V[static_ap_counter] + delta_T*V_dt);
        }

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);

        static_ap_counter++; 

        myfile.close();

    return(0);
}

double output_WT_Static_AP(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    vec_V.push_back(V_start);
    vec_n.push_back(dynamical_n(0));
    vec_m.push_back(dynamical_m(0));
    vec_h.push_back(dynamical_h(0));

    for(double i = 0; i < 1000; i++){
        Static_WT_AP(36);
    }

    reset_vecs(0);
    static_ap_counter = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(dynamical_n(0));
    vec_m.push_back(dynamical_m(0));
    vec_h.push_back(dynamical_h(0));

    for(double i = 0; i < 1000; i++){
        Static_WT_AP(0);
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_V2.size());
    sizes.insert(sizes.begin(),vec_K_I.size());
    sizes.insert(sizes.begin(),vec_Na_I.size());
    sizes.insert(sizes.begin(),vec_L_I.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool V; 
    bool V_2;
    bool K_I; 
    bool Na_I; 
    bool L_I; 
    
    cout << "Break point 4" << endl;

    // myfile << "V,K_I,Na_I,L_I\n";
    // for (int i = 0; i < max_size; i++)
    // {
    //     //cout << "Break point 5" << endl;
    //     V = (vec_V.size() > i) ? true : false;
    //     K_I = (vec_K_I.size() > i) ? true : false;
    //     Na_I = (vec_Na_I.size() > i) ? true : false;
    //     L_I = (vec_L_I.size() > i) ? true : false;

    //     //cout << "Break point 6" << endl;

    //     if(V) myfile << vec_V[i] - 70 <<"," ;
    //     if(!V) myfile <<"," ;
    //     if(K_I) myfile << vec_K_I[i] << ",";
    //     if(!K_I) myfile << ",";
    //     if(Na_I) myfile << vec_Na_I[i] << ",";
    //     if(!Na_I) myfile <<",";
    //     if(L_I) myfile << vec_L_I[i];

    //     //if(!dn_V_0) myfile << vec_tiny_N[i];

    //     myfile << "\n";
    //     //cout << "Break point 6" << endl;
    // }

    myfile << "V,V2\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        V = (vec_V.size() > i) ? true : false;
        V_2 = (vec_V2.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(V) myfile << vec_V[i] - 70 <<"," ;
        if(!V) myfile << "," ;
        if(V_2) myfile << vec_V2[i] - 70;

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }


    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_WT_Static_AP(0);
  cout << "End" << endl;
}
