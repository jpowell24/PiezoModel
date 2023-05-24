#include "2d_Calcium_Dynamics.h"

const char *path1="../data_files/2d_Piezo_Channel.csv";
const char *path2="../data_files/2d_Piezo_Channel_avg.csv";

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,0.25);

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    vec_coords.clear();
    vec_time.clear();
    vec_num_open.clear();
    vec_num_closed.clear();
    vec_Piezo_current.clear();
    vec_buff_bound.clear();
    vec_w.clear();

    return(0);
}

vector<vector<double> > fill_2dvecs(int x, int y, double value){
    // this is an optional method to fill 3D vectors with some value
    // this may be useful if we want an editable size in the future
    for(double i = 0; i <= x; i++){
        for(double j = 0; j <= y; j++){ 
            vec_y.push_back(value);
        }
    vec_coords.push_back(vec_y);
    }

    return(vec_coords);

    reset_vecs(0);
}

double PotentialE(double out, double in, int Z) { //calculated in VOLTS NOT MILLIVOLTS
  double E = (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //cout << "THIS IS E in V: " << E << endl;
  return (E);
}

double Compute_J_diffusion(int time, int x, int y) { 
    double R = 1;
    double total_diffusion, x_diffusion, y_diffusion;
    if((x >= 1) && (x <= x_max - 1)){
        x_diffusion = (1/R)*(vec_time[time][y][x + 1] - (2 * vec_time[time][y][x]) + vec_time[time][y][x - 1]);
    }
    if((y >= 1) && (y <= y_max - 1)){
        y_diffusion = (1/R)*(vec_time[time][y + 1][x] - (2 * vec_time[time][y][x]) + vec_time[time][y - 1][x]);
    }
    
    if((y == 0)){
        y_diffusion = (1/R)*(vec_time[time][y + 1][x] - (vec_time[time][y][x]));
    }
    if((y == y_max)){
        y_diffusion = (1/R)*(vec_time[time][y - 1][x] - (vec_time[time][y][x]));
    }

    if((x == 0)){
        x_diffusion = (1/R)*(vec_time[time][y][x + 1] - (vec_time[time][y][x]));
    }
    if((x == x_max)){
        x_diffusion = (1/R)*(vec_time[time][y][x - 1] - (vec_time[time][y][x]));
    }

    total_diffusion = x_diffusion + y_diffusion;
  return(total_diffusion);
}


double Compute_J_on(double C_cyt, int time, int x, int y){
  double scaling_factor = 10000000; 
  double local_C_cyt = C_cyt*scaling_factor;

  k_buff_bind = 0.600; //for the buffer BAPTA in mM
  k_buff_unbind = 0.100;

  J_on = k_buff_bind*local_C_cyt*buff_unbound;
  J_off = k_buff_unbind*buff_bound; 
  buff_c_dT = k_buff_bind*local_C_cyt*buff_unbound - k_buff_unbind*buff_bound;
  vec_buff_bound[time + 1][x][y] = (vec_buff_bound[time][x][y] + buff_c_dT);

  buff_diff = J_off - J_on;
  buff_diff = buff_diff*pow(scaling_factor,-1);
  //vec_J_on.push_back(buff_diff);

  buff_counter++;
  return(buff_diff*0.00000001); //*0.000001
}

double Compute_J_serca(double serc_local, int time, int x, int y){
  double local_C_cyt = 100000*serc_local;

  J_serca = v_serca*(pow(local_C_cyt,2)/(pow(local_C_cyt,2) + pow(K_p,2)));
  if(J_serca != J_serca){
    J_serca = 0;
  }
  
  // vec_J_serca.push_back(J_serca*0.0001); //this and ryr are scaled weirdly, I don't know why this works better - otherwise Piezo will dominate
  
  return(J_serca*10);
}

double Compute_J_ryr(double ryr_local, int time, int x, int y){ // I am almost certain that there is something wrong with the kinetic equations that go beyond the paper

  double local_C_cyt = 100000000*ryr_local; // this is here in case we want to scale 

  w_inf = ((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b))/((1/K_c) + (K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)); 

  tau_w = w_inf/K_d;
  w_dt = (w_inf - vec_w[time][x][y])/tau_w;

  vec_w[time + 1][x][y] = (vec_w[time][x][y] + w_dt);
  P_open = (vec_w[time + 1][x][y]*((1 + pow(local_C_cyt,3))/K_b))/((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b));
  J_ryr = (v_rel*P_open + v_leak)*(C_er - local_C_cyt); //I do not like the C_er - local_c_cyt

  if(J_ryr != J_ryr){
    J_ryr = 0;
  }

  w_counter++;

    double scale_down = 1/100000000;

  return(J_ryr*scale_down);
}

double Piezo_Channel(double potential, int time, int x, int y){

    int stochastic = stochastic_opening(generator);
    double open_temp;
    double P_open_temp;

    open_temp = vec_num_open[time][x][y] + abs(stochastic);

    if (open_temp >= N_Piezo_channels){
        open_temp = N_Piezo_channels;
    }

    int local_tau = 1.6/delta_T;
    if (!(open_counter % local_tau)){
        open_temp = open_temp*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
    }

    vec_num_open[time + 1][x][y] = open_temp;
    vec_num_closed[time + 1][x][y] = 1 - open_temp;
    
    G_Piezo_total = open_temp*G_Piezo_single;
    Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
    vec_Piezo_current[time + 1][x][y] = Piezo_current;

    return(Piezo_current*0.1);
}

double Calcium_concentration(double x){

    cout << "  High" << endl;

    E_Ca = PotentialE(0.0024, 0.00000012, 2);
    
    double divs = (x_max + 1)*(y_max + 1);

    double mols_divs = 0.0000000000001/divs;

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = mols_divs;
            vec_num_closed[0][i][j] = N_Piezo_channels;
            vec_Piezo_current[0][i][j] = 0;
            vec_buff_bound[0][i][j] = 0;
            vec_w[0][i][j] = 0;
        }
    }

    //cout << "break point 2" << endl;

    double scaling_factor = 1;
    double avg_temp;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){
        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                C_cyt = vec_time[time_temp][i][j];

                //Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca, time_temp, i, j) + scaling_factor*Compute_J_ryr(C_cyt, time_temp, i, j) - Compute_J_serca(C_cyt, time_temp, i, j) + Compute_J_on(C_cyt, time_temp, i, j));
                
                // scaling_factor*Piezo_Channel(E_Ca, time_temp, i, j) +
                Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca, time_temp, i, j) + Compute_J_diffusion(time_temp, j, i) + Compute_J_on(C_cyt, time_temp, i, j) + scaling_factor*Compute_J_ryr(C_cyt, time_temp, i, j) - Compute_J_serca(C_cyt, time_temp, i, j));
                
                vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;

                // what will be better is to solve for the number of moles in each cube, then use that to calculate overall concentration
                avg_temp += vec_time[time_temp][i][j];
                //cout << avg_temp << endl;
            }
        }
        vec_average.push_back(avg_temp/(x_max*y_max));
        //cout << avg_temp/(x_max*y_max) << endl;
        avg_temp = 0;
    }

    cout << "    to Low" << endl;

    return(0);
}


double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    Calcium_concentration(0);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_time.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector length: " << max_size << endl;

    for(int time = 0; time <= time_max; time++){
        for (int x = 0; x <= x_max; x++){
            for (int y = 0; y <= y_max; y++){  
                if(y < y_max) {
                    myfile << vec_time[time][x][y] << ",";
                }  
                else {
                    myfile << vec_time[time][x][y];
                }
            }
        myfile << "\n";
        }
    myfile << "\n";
    }

    reset_vecs(0);

    myfile.close();

    return (x);
}

double output_avg_file(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_average.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector length: " << max_size << endl;

    bool bool_average;

    //cout << "Break point 4" << endl;

    myfile << "Average\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        bool_average = (vec_average.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(bool_average) myfile << vec_average[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    reset_vecs(0);

    myfile.close();

    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  output_avg_file(0);
  cout << "End" << endl;
}
