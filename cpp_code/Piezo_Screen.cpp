#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/Piezo_Channel.csv";

//CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double log_convert = 2.303; // to convert from ln to log10
double R_constant = 8.1345;
double F = 96485.3321;
double body_temp = 310.15;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//General use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double delta_T = 0.001; // this is in ms
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// calcium concentration: 2.4mM outside, 100nM inside https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553253/#:~:text=Extracellular%20calcium%2C%20and%20particularly%20the,calcium%20of%201.1%E2%80%931.4%20mM
double E_Ca; // 131.373 --> this is for humans, i.e., body temp of 310K etc. Unsure what it is for Drosophila

// Piezo Kinetics %%%%%%%%%%%%%%%%%%%%%%%%%
double G_Piezo_single = 0.000000000030; // something like 30 pS
double G_Piezo_total;

vector<double> temp_vec, vec_channel_number, vec_time;

double PotentialE(double out, double in, int Z) {
  double E = 1000 * (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "Calculated E_Ca in mV: " << E << endl;
  return (E);
}

double Piezo_screen(double channel_number, double time_max){
  E_Ca = PotentialE(0.0024, 0.0000001, 2); //Calculates potential of Ca, taking it to be 100nM inside, 2.4mM outside, also this outputs in millivolts
  double Charge_temp = 0;
  double Charge_temp_dt = 0;
  double Conc_temp = 0;
  for(double x = 0; x <= channel_number; x += 10){ //scans through the number of channels
    for(double y = 0; y <= time_max; y += delta_T){ // scans through time is in milliseconds, so step with 1/1000 of ms seems good
      Charge_temp_dt = x*(E_Ca/1000)*G_Piezo_single*delta_T*0.001; // the 0.0000001 was added to convert from seconds to 0.001 ms, E_ca is divided by 1000 to convert mV to V
      Charge_temp += Charge_temp_dt; 
      
      Conc_temp = (Charge_temp/F)/((4/3)*M_PI*pow(0.00175,3)); //the radius, in this case, is in centimeters. I do not know why, but this seems to work better

      if((Conc_temp >= 0.0000000022) && (Conc_temp <= 0.0000000023)){
        vec_time.push_back(y); //stores the time that it took to reach 2.3nM
        vec_channel_number.push_back(x); //stores # of channels
        //cout << "occured" << endl;
      }
    }
    Charge_temp = 0; //resets after each time scan 
    Conc_temp = 0; 
    Charge_temp_dt = 0;
  }
  return(0);
}

double voltage_output(double x)
{

    Piezo_screen(1000, 10);

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_time.size());
    sizes.insert(sizes.begin(),vec_channel_number.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector size: " << max_size << endl;

    bool time_bool;
    bool channel_number_bool;

    //cout << "Break point 4" << endl;

    myfile << "Time,Channel_number\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        time_bool = (vec_time.size() > i) ? true : false;
        channel_number_bool = (vec_channel_number.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(time_bool) myfile << vec_time[i] << ",";
        if(!time_bool) myfile << ",";
        if(channel_number_bool) myfile << vec_channel_number[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  voltage_output(0);
  cout << "End" << endl;
}