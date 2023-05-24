#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

int x_max = 50;
int y_max = 50;
double divs = (x_max + 1)*(y_max + 1);
double mols_divs = 0.0000000012/divs;

int time_max = 10;
double delta_T = 0.01;
int time_max_calc = time_max/delta_T; 
double divide = (y_max + 1) * (x_max + 1);
double size_scale = 1/divide; 

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 2d specific thingsc
vector<double> vec_x;
vector<double> vec_y; 
vector<vector<double> > vec_coords;
vector<vector<vector<double> > > vec_time(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
vector<double> vec_average;
//vector<vector<vector<double> > > vec_time;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double log_convert = 2.303; // to convert from ln to log10
double R_constant = 8.1345;
double F = 96485.3321;
double body_temp = 310.15;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//General use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// calcium concentration: 2.4mM outside, 100nM inside https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553253/#:~:text=Extracellular%20calcium%2C%20and%20particularly%20the,calcium%20of%201.1%E2%80%931.4%20mM
double Ca_in, Ca_out;
double E_Ca; // 131.373 --> this is for humans, i.e., body temp of 310K etc. Unsure what it is for Drosophila
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Piezo Kinetics %%%%%%%%%%%%%%%%%%%%%%%%%
double G_Piezo_single = 0.000000000030; 
double G_Piezo_total;
int N_Piezo_channels = 10000;
double p_open = 0; 
vector<vector<vector<double> > > vec_num_open(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
double p_closed = 1;
vector<vector<vector<double> > > vec_num_closed(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
double Piezo_current;
vector<vector<vector<double> > > vec_Piezo_current(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
int open_counter = 0;
vector<double> vec_P_Substrate;
vector<double> vec_P_Pressure;
vector<double> vec_P_Voltage;
double avg_piezo_temp;
vector<double> vec_Piezo_average;
// Piezo1 = 29pS https://www.sciencedirect.com/science/article/pii/S0968000416301505
// The decay rate, according to this paper, is 1ms
// 25-30pS https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/dvdy.401
// Time constant = 6.2 ± 0.3 ms https://www.nature.com/articles/nature10812
// This paper also says conductance is 30pS and that Drosophila is closer to 3.3pS....
// Piezo1 (16ish ms) vs Piezo2 (7ish ms) inactiation kinetics https://www.science.org/doi/full/10.1126/science.1193270?casa_token=WnOGIrz8PbAAAAAA%3AdA422PZ9tigaVHLD6IxvFlWpQer0KOFm_jvCrohrhgnSzzdEPOXhZo_koo-cZeWbLeNSTFPtSJ3bdtA
// Force gating for Piezo1 is about –30 mmHg or 1.4 mN/m https://www.sciencedirect.com/science/article/pii/S0968000421000220#bb0270
// Piezo1 half-maximal activation T50 = 2.7 ± 0.1 mN/m https://elifesciences.org/articles/12088
// Piezo1 P50 of –28.1 ± 2.8 and –31.2 ± 3.5 mmHg https://www.science.org/doi/full/10.1126/science.1193270?casa_token=gKZRiU1R_vAAAAAA%3AAuhBFZP-QjXckn7G9aBL6A_ZJfCjsRrUIqoXCbBA1887i29wPWtzlMBfdwShr45kBM7Pj-N4NYVlfEQ
// Piezo1 P50 of -60mV voltage gating https://www.nature.com/articles/s41467-018-03502-7
double T_Piezo = 0.016; //16 ms tau

// Piezo dimensions %%%%%%%%%%%%%%%%%%%%%%%
// Blade = 200A https://www.nature.com/articles/nature25453
// SA, then, is 17320.51 A^2, or 173.20 nm^2, or 72658.95
// If we assume the growth cone of a sphere is 1um, then the SA is 12.57 um^2
// I think that means that 72658.95 Piezo proteins can be on the surface

// Overall Definitions %%%%%%%%%%%%%%%%%%%%%
vector<double> vec_Ca_conc; //calcium concentration
double Ca_c_dT; //derivative of concentration of calcium per time in the cytoplasm
double D_diff_Ca; //coefficient of diffusion
double Ca_c_dT_dT; //2nd derivative of Ca concentration per unit time
double J_ipr; //IPR is the IP3 receptor, and is supposedly the biggest calcium leaker
// Note: it does not seem that Piezo upregulating PLC/IP3 has been shown
double cone_radius = 0.0002; //this is in ... centimeters ... for some reason
// update from later Jackson (whose brain hath grown). cm is best because you can convert to mL easily, which is relevant for molar quantity
// also, most constant values, like capacitance are solved in terms of cm
double cone_circumference; //circumference of the growth cone
double cone_cross_area; //cross sectional area of the growth cone
double J_in; //flux of Ca2+ into the cell 
double J_pm; //flux pumped out of the cell
double J_ryr; //flux from RyR receptor
double J_serca; //flux from serca pump, maybe includes good SERCA model: https://link.springer.com/article/10.1140/epjp/s13360-023-03691-1
double J_on; //function of binding of Ca2+ to buffers, Page 309 of mathematical physiology seems good
double J_off; //function of unbinding of Ca2+ from buffers
// Ca_c_dT = D_c*Ca_c_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + J_ryr - J_serca - J_on + J_off; 

// RyR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Model 1
// Very good reference for calcium dynamics: https://www.frontiersin.org/articles/10.3389/fphys.2012.00114/full#F1
double vol_D; //dyadic space volume
double vol_ER; //ER volume
double g_ryr; //RyR channel conductance
double C_er = 0.0001; //concentration inside the ER; units of M; in this case, we are assuming that C_er is relatively constant
double C_cyt; //concentration in the cytoplasm
double N_ryr; //stochastic number of RyR channels
// J_ryr = N_ryr*(g_ryr/vol_D)*(C_er - C_cyt); //This kind of decribes the local movement due to gradient

// RyR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Model 2
// From this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4648987/pdf/10827_2015_Article_579.pdf
double v_rel = 0.005; // no clue what this is, just some constant, units of ms^-1
double v_leak = 0.00012; // still, no clue, just some constant, units of ms^-1
// c_cyt and c_er are still used 
double P_open; // probability/proportion that channels are open
double w; // dynamical variable -- I don't want to talk about it
double w_inf; // infinite version of w; w_inf is a state-function and does not need a vector
double w_dt; // w derivative with respect to time
double tau_w; // time constant of w
double K_a = 0.0192; // kinetics a, units of uM^4
double K_b = 0.2573; // kinetics b, units of uM^3
double K_c = 0.0571; // kinetics c, unitless
double K_d = 0.0001; // kinetics d, units of ms^-1
int w_counter = 0;

vector<vector<vector<double> > > vec_J_ryr(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
vector<vector<vector<double> > > vec_w(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
// J_ryr = (v_rel*P_open + v_leak)(C_er - C_cyt);
// P_open = (w*((1 + C_cyt^3)/K_b))/((K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b));
// w_inf = ((K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b))/((1/K_c) + (K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b));
// w_dt = (w_inf - w)/tau_w; 
// tau_w = w_inf/K_d;
// w = w + w_dt;

// SERCA Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Most SERCA models seem to simply be the Hill function
double v_serca = 0.12; // some constant, units of muM / ms
double K_p = 0.3; // kinetics p, units of uM
vector<vector<vector<double> > > vec_J_serca(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));

// TECHNICALLY THiS MODEL IS OUTDATED, CHECK PAGE 284 OF MATHEMATICAL PHYS
// J_serca = v_serca*((C_cyt^2)/(C_cyt^2 + K_p^2));

// IPR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Taking this mostly from Mathematical Physiology page 293. This is a simplified IPR model
// It is difficult to know if this will be sufficient 
// I am afraid that the addition of IPR will completely dominate the code, and it is unknown if this is true to the growth cone

// Buffering Definitions %%%%%%%%%%%%%%%%%
// Reference includes a list of models published by year: https://www.frontiersin.org/articles/10.3389/fncom.2018.00014/full
double buff_unbound = 0.02*size_scale; //concentration of unbound buffer, which we are taking to be b_total
double buff_bound = 0.02*size_scale; //concentration of bound buffer
double buff_total = buff_unbound + buff_bound;
vector<vector<vector<double> > > vec_buff_bound(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
vector<vector<vector<double> > > vec_buff_unbound(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
int buff_counter = 0;
double k_buff_bind; //binding affinity/Kon of buffer
double k_buff_unbind; //unbinding affinity/Koff of buffer
double buff_c_dT; //derivative of concentration of buffer with respect to time
double buff_c_dT_dT; //second derivative of concentration of buffer with respect to time, we do not really ned this unless we are measuring the diffusion of the buffer
double buff_diff; //store diff between jon and joff
vector<double> vec_buff_average; // using to store the average value of buffer bound over the whole area
vector<double> vec_ubuff_average; 
vector<vector<vector<double> > > vec_J_on(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1)));
// J_on = k_buff_bind*C_cyt*buff_unbound;
// J_off = k_buff_unbind*buff_bound; 
// buff_c_dT = k_buff_bind*C_cyt*buff_unbound - k_buff_unbind*buff_bound;
// buff_bound = buff_bound + buff_c_dT;