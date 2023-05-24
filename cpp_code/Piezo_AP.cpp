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


double Static_WT_AP(double Ca_input){

    //vec_Cl_I.push_back(0);
    
        current = 10000000*Ca_input;
        //cout << x << endl; 

        // if(static_ap_counter >= 5000){
        //   current = 5;
        // }

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[static_ap_counter]);
        dynamical_h(vec_V[static_ap_counter]);
        dynamical_n(vec_V[static_ap_counter]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[static_ap_counter] + delta_T*((vec_inf_n[static_ap_counter] - vec_n[static_ap_counter])/vec_tau_n[static_ap_counter]));
        vec_m.push_back(vec_m[static_ap_counter] + delta_T*((vec_inf_m[static_ap_counter] - vec_m[static_ap_counter])/vec_tau_m[static_ap_counter]));
        vec_h.push_back(vec_h[static_ap_counter] + delta_T*((vec_inf_h[static_ap_counter] - vec_h[static_ap_counter])/vec_tau_h[static_ap_counter]));

        //cout << "Break point 3" << endl;

        K_I_temp = (g_k*pow(vec_n[static_ap_counter+1],4)*((vec_V[static_ap_counter]) - E_k));
        Na_I_temp = (g_Na*pow(vec_m[static_ap_counter+1],3)*pow(vec_h[static_ap_counter+1],1)*((vec_V[static_ap_counter]) - E_Na));
        L_I_temp = (g_l*((vec_V[static_ap_counter]) - E_l));

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp)/C_m;
        vec_V.push_back(vec_V[static_ap_counter] + delta_T*V_dt);

        vec_VWT.push_back(vec_V[static_ap_counter] + delta_T*V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);

        //cout << V << endl;
        //cout << V_temp << endl;
        static_ap_counter++; 
        //cout << x << endl;
    return(0);
}

double output_WT_Static_AP(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    vector<int> sizes;


    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_K_I.size());
    sizes.insert(sizes.begin(),vec_Na_I.size());
    sizes.insert(sizes.begin(),vec_L_I.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool V; 
    bool K_I; 
    bool Na_I; 
    bool L_I; 
    
    //cout << "Break point 4" << endl;

    myfile << "V,K_I,Na_I,L_I\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        V = (vec_V.size() > i) ? true : false;
        K_I = (vec_K_I.size() > i) ? true : false;
        Na_I = (vec_Na_I.size() > i) ? true : false;
        L_I = (vec_L_I.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(V) myfile << vec_V[i] <<"," ;
        if(!V) myfile <<"," ;
        if(K_I) myfile << vec_K_I[i] << ",";
        if(!K_I) myfile << ",";
        if(Na_I) myfile << vec_Na_I[i] << ",";
        if(!Na_I) myfile <<",";
        if(L_I) myfile << vec_L_I[i];

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    myfile.close();
    return (x);
}

double PotentialE(double out, double in, int Z) { //calculated in VOLTS NOT MILLIVOLTS
  double E = (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "\n THIS IS E in V: " << E << endl;
  return (E);
}

double Compute_J_on(double C_cyt){
  double scaling_factor = 100000; 
  double local_C_cyt = C_cyt*scaling_factor;
  //cout << "Compute_J_on active" << endl;
  k_buff_bind = 0.600; //for the buffer BAPTA in mM
  k_buff_unbind = 0.100;

  J_on = k_buff_bind*local_C_cyt*buff_unbound;
  J_off = k_buff_unbind*buff_bound; 
  buff_c_dT = k_buff_bind*local_C_cyt*buff_unbound - k_buff_unbind*buff_bound;
  vec_buff_bound.push_back(vec_buff_bound[buff_counter] + buff_c_dT);

  buff_diff = J_off - J_on;
  buff_diff = buff_diff*pow(scaling_factor,-1);
  vec_J_on.push_back(buff_diff);

  buff_counter++;
  return(buff_diff*0.000001); //*0.000001
}

double Compute_J_serca(double serc_local){
  //cout << "Compute_J_serca active" << endl;
  double local_C_cyt = 1000*serc_local;

  J_serca = v_serca*(pow(local_C_cyt,2)/(pow(local_C_cyt,2) + pow(K_p,2)));
  if(J_serca != J_serca){
    J_serca = 0;
  }
  vec_J_serca.push_back(J_serca*0.01); //this and ryr are scaled weirdly, I don't know why this works better - otherwise Piezo will dominate
  return(J_serca*0.01);
}

double Compute_J_ryr(double ryr_local){ // I am almost certain that there is something wrong with the kinetic equations that go beyond the paper
  //cout << "Compute_J_ryr active1" << endl;

  double local_C_cyt = 1000*ryr_local; // this is here in case we want to scale 

  w_inf = ((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b))/((1/K_c) + (K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)); 
  //cout << (K_a/pow(local_C_cyt,4)) << endl;

  //cout << "Compute_J_ryr active2" << endl;

  tau_w = w_inf/K_d;
  w_dt = (w_inf - vec_w[w_counter])/tau_w;

  //cout << "Compute_J_ryr active3" << endl;

  vec_w.push_back(vec_w[w_counter] + w_dt);
  P_open = (vec_w[w_counter + 1]*((1 + pow(local_C_cyt,3))/K_b))/((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b));
  J_ryr = (v_rel*P_open + v_leak)*(C_er - local_C_cyt); //I do not like the C_er - local_c_cyt

  if(J_ryr != J_ryr){
    J_ryr = 0;
  }

  w_counter++;

  vec_J_ryr.push_back(J_ryr*0.01);
  return(J_ryr*0.01);
}

double Piezo_Channel(double potential){
  //cout << "Piezo_Channel active" << endl;
  int stochastic = stochastic_opening(generator);
  int open_temp;
  double P_open_temp;

  open_temp = vec_num_open[open_counter] + abs(stochastic);

  // cout << P_open_temp << endl;
  // cout << exp((delta_T*0.001)/T_Piezo) << endl;

  if (open_temp >= N_Piezo_channels){
    open_temp = N_Piezo_channels;
  }

  int local_tau = 1.6/delta_T;
  //cout << local_tau << endl;
  //cout << open_counter % local_tau << endl;
  if (!(open_counter % local_tau)){

    //cout << open_counter << endl;
    open_temp = open_temp*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
  }

   if ((open_counter >= 50000) && (open_counter <= 50010)){
     open_temp = N_Piezo_channels; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
   }


  //cout << open_temp << endl;

  vec_num_open.push_back(open_temp);  
  vec_num_closed.push_back(1 - open_temp);
    
  G_Piezo_total = open_temp*G_Piezo_single;

  Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
  vec_Piezo_current.push_back(Piezo_current);

  open_counter++;
  return(Piezo_current);
}

double Calcium_concentration(double time_range, double delta_T){
  vec_V.push_back(V_start);
  vec_n.push_back(dynamical_n(0));
  vec_m.push_back(dynamical_m(0));
  vec_h.push_back(dynamical_h(0));

  E_Ca = PotentialE(0.0024, 0.0000001, 2);
  vec_w.push_back(0);
  vec_num_open.push_back(0);
  vec_num_closed.push_back(N_Piezo_channels);
  vec_buff_bound.push_back(0);
  vec_Ca_conc.push_back(0.00000012); //supposed to be 120nM
  // cone_circumference = 2*M_PI*cone_radius;
  // cone_cross_area = M_PI*pow(cone_radius,2);

  // Ca_c_dT = D_diff_Ca*Ca_c_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + Compute_J_ryr(C_cyt) - Compute_J_serca(C_cyt) + Compute_J_on(C_cyt);

  double scaling_factor = 1;

  for(double i = 0; i <= time_range; i += delta_T){
    C_cyt = vec_Ca_conc[Ca_counter]; 

    //cout << C_cyt << endl;

    Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca) + scaling_factor*Compute_J_ryr(C_cyt) - scaling_factor*Compute_J_serca(C_cyt) + Compute_J_on(C_cyt));
    vec_Ca_conc.push_back(vec_Ca_conc[Ca_counter] + Ca_c_dT);
    Ca_counter++;

    Static_WT_AP(vec_Ca_conc[Ca_counter]);
    
  }
  return(0);
}

double voltage_output(double x)
{

    Calcium_concentration(1000, delta_T);
    output_WT_Static_AP(0);

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_Ca_conc.size());
    sizes.insert(sizes.begin(),vec_Piezo_current.size());
    sizes.insert(sizes.begin(),vec_J_on.size());
    sizes.insert(sizes.begin(),vec_J_ryr.size());
    sizes.insert(sizes.begin(),vec_J_serca.size());
    sizes.insert(sizes.begin(),vec_num_open.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool bool_Ca_conc;
    bool bool_Piezo_current;
    bool bool_J_on;
    bool bool_J_ryr;
    bool bool_J_serca;
    bool bool_num_open;

    //cout << "Break point 4" << endl;

    myfile << "Ca_concentration,Buffering,Piezo_current,J_ryr,J_serca,Piezo_open\n";
    for (int i = 0; i < max_size; i++)
    {
      if(i >= 1000){
        //cout << "Break point 5" << endl;
        bool_Ca_conc = (vec_Ca_conc.size() > i) ? true : false;
        bool_Piezo_current = (vec_Piezo_current.size() > i) ? true : false;
        bool_J_on = (vec_J_on.size() > i) ? true : false;
        bool_J_ryr = (vec_J_ryr.size() > i) ? true : false;
        bool_J_serca = (vec_J_serca.size() > i) ? true : false;
        bool_num_open = (vec_num_open.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(bool_Ca_conc) myfile << 1000000000*vec_Ca_conc[i] << ",";
        if(!bool_Ca_conc) myfile << ",";
        if(bool_J_on) myfile << 1*vec_J_on[i] << ",";
        if(!bool_J_on) myfile << ",";
        if(bool_Piezo_current) myfile << 1*vec_Piezo_current[i] << ",";
        if(!bool_Piezo_current) myfile << ",";
        if(bool_J_ryr) myfile << 1*vec_J_ryr[i] << ",";
        if(!bool_J_ryr) myfile << ",";
        if(bool_J_serca) myfile << -1*1*vec_J_serca[i]<< ",";;
        if(!bool_J_serca) myfile << ",";
        if(bool_num_open) myfile << vec_num_open[i];


        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
      }
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  voltage_output(0);
  cout << "End" << endl;
}
