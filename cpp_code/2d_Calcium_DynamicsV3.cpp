#include "2d_Calcium_Dynamics.h"

const char *path1="../data_files/2d_Piezo_Channel.csv";
const char *path2="../data_files/2d_Piezo_Channel_avg.csv";

default_random_engine generator; //important: must be outside of the loop/method that calls it

normal_distribution<double> stiffness(0.7,0.1);
normal_distribution<double> pressure(30,5);
normal_distribution<double> voltage(-70,10);

int reset_vecs(int x){
    vec_time.clear();
    vec_buff_bound.clear();
    vec_w.clear();

    return(0);
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

    // NOTE: diffusion is not scaled to the size of the cell! this is very important!! 
    // in the current state, a huge cell will have an over-inflated Ca2+ storage
    // This is because it is more difficult to diffuse through the greater divisions


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

    //cout << (x_diffusion + y_diffusion) << " " << (divs/100)*(x_diffusion + y_diffusion) << endl;

    total_diffusion = 1000*(divs/100)*(x_diffusion + y_diffusion);

  return(total_diffusion);
}

double Piezo_P_Pressure(double i){
    double F_inf;
    F_inf = 1/(exp((30 - i)/6) + 1);
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

double Piezo_Channel(double potential, int time, int x, int y, int loc){

    double stochastic_stiffness = stiffness(generator);
    double stochastic_pressure = pressure(generator);
    double stochastic_voltage = voltage(generator);

    //cout << stochastic_stiffness << " " << stochastic_pressure << " " << stochastic_voltage << endl;

    int local_N_Piezo, open_local, closed_inactive_local, closed_active_local; 

    double P_opening_temp;
    double Pressure_input, Substrate_input, Voltage_input; 

    Pressure_input = stochastic_pressure;
    Substrate_input = stochastic_stiffness;
    Voltage_input = stochastic_voltage;

    local_N_Piezo = N_Piezo_channels/(2*x_max + 2*y_max);

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);

    double P_total = P_P*P_S + P_V;

    P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1);

    open_local = P_opening_temp*local_N_Piezo;

    // cout << P_opening_temp << endl;
    //open_local = vec_num_open[time][x][y] + P_opening_temp*closed_active_local;

    if(open_local >= local_N_Piezo){
        open_local = local_N_Piezo;
    }

    int local_tau = 1.6/delta_T;
    
    if(!(open_counter % local_tau)){
        closed_inactive_local = vec_num_open[time][x][y]*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
    } //this needs work 

    open_local = open_local - closed_inactive_local;

    if(loc == 2){
        loc = 1;
    }

    open_local = open_local*loc;
    //cout << open_temp << endl;

    vec_num_open[time + 1][x][y] = open_local;
    vec_num_closed[time + 1][x][y] = closed_inactive_local;
    
    G_Piezo_total = open_local*G_Piezo_single;

    //cout << G_Piezo_total << endl;

    Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
    vec_Piezo_current[time + 1][x][y] = Piezo_current;

    //cout << loc << " equals " << Piezo_current << endl;

    return(Piezo_current);
}

double Compute_J_on(double C_cyt, int time, int x, int y){
    buff_bound = vec_buff_bound[time][x][y]; 
    buff_unbound = vec_buff_unbound[time][x][y]; 

    double local_C_cyt = C_cyt;
    //double local_C_cyt = 1.1;

    k_buff_bind = 0.600; //for the buffer BAPTA in mM
    k_buff_unbind = 0.100;

    J_on = (local_C_cyt/mols_divs)*k_buff_bind*buff_unbound;
    J_off = (mols_divs/local_C_cyt)*k_buff_unbind*buff_bound;
    buff_c_dT = (local_C_cyt/mols_divs)*k_buff_bind*buff_unbound - (mols_divs/local_C_cyt)*k_buff_unbind*buff_bound;
    
    //buff_c_dT = -(k_buff_unbind*buff_bound);

    vec_buff_bound[time + 1][x][y] = (vec_buff_bound[time][x][y] + delta_T*buff_c_dT);
    vec_buff_unbound[time + 1][x][y] = (vec_buff_unbound[time][x][y] - delta_T*buff_c_dT);

    if(vec_buff_unbound[time + 1][x][y] < 0){
        vec_buff_unbound[time + 1][x][y] = 0;
    }
    else if(vec_buff_unbound[time + 1][x][y] > buff_total){
        vec_buff_unbound[time + 1][x][y] = buff_total;
    }
    if(vec_buff_bound[time + 1][x][y] < 0){
        vec_buff_bound[time + 1][x][y] = 0;
    }
    else if(vec_buff_bound[time + 1][x][y] > buff_total){
        vec_buff_bound[time + 1][x][y] = buff_total;
    }

    buff_diff = delta_T*(J_off - J_on);
    //buff_diff = J_off;
    //buff_diff = buff_diff*pow(scaling_factor,-1);
    //vec_J_on.push_back(buff_diff);
    //cout << time << " " << buff_diff << endl;
    return(buff_diff); //*0.000001
}

// maybe we can just do a simple PID controller for regulating the ER concentration
// and set the efflux and influx rates to the SERCA pump
// i guess RyR will need to be independent, obviously

double Compute_efflux(double C_cyt, int time, int x, int y, int loc){
    // efflux can only occur at the barriers, which makes things a bit weirder

    double efflux;
    double proportional_error; 
    double deriative_error; 
    double integral_error; 

    if(C_cyt > mols_divs){
        efflux = -1000*loc*(C_cyt - mols_divs);
    }
    else{
        efflux = 0;
    }

    //efflux = -0.01*C_cyt; 
    
    return(efflux);
}


double Calcium_concentration(double x){

    cout << "  High" << endl;

    E_Ca = PotentialE(0.0024, 0.00000012, 2);

    cout << divs << endl;
    cout << mols_divs << endl;

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = mols_divs;
            vec_buff_bound[0][i][j] = buff_bound;
            vec_buff_unbound[0][i][j] = buff_unbound;
            vec_num_closed[0][i][j] = N_Piezo_channels;
            vec_num_open[0][i][j] = 0;
            vec_Piezo_current[0][i][j] = 0;
        }
    }

    double scaling_factor;
    double divide = (y_max + 1)*(x_max + 1);
    scaling_factor = 1/divide; 
    double avg_temp, avg_buff_temp, avg_ubuff_temp;

    for(int time_temp = 0; time_temp <= time_max_calc; time_temp++){
        //cout << "break point 4 " << time_temp << endl;
        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                int location;

        
                if(((i == 0) || (i == x_max)) && ((j == 0) || (j == y_max))){
                    location = 2;
                    //cout << true;
                }
                else if(((i == 0) || (i == x_max))){
                    location = 1;
                    //cout << true;
                }
                else if((j == 0) || (j == y_max)){
                    location = 1;
                    //cout << true;
                }
                else{
                    location = 0;
                    //cout << false;
                }

                C_cyt = vec_time[time_temp][i][j];

                Ca_c_dT = delta_T*(scaling_factor*Compute_J_diffusion(time_temp, j, i) + scaling_factor*Compute_J_on(C_cyt, time_temp, i, j) + scaling_factor*Compute_efflux(C_cyt, time_temp, i, j, location) + Piezo_Channel(E_Ca, time_temp, i, j, location));
                //vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;
                vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;
                

                // if((vec_time[time_temp][i][j] + Ca_c_dT) >= 0){
                //     vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;
                // }
                // else{
                //     vec_time[time_temp + 1][i][j] = 0; 
                //     cout << "something is wrong, horribly wrong" << endl;
                // }

                avg_temp += vec_time[time_temp][i][j];
                avg_buff_temp += vec_buff_bound[time_temp][i][j];
                avg_ubuff_temp += vec_buff_unbound[time_temp][i][j];
                avg_piezo_temp += vec_Piezo_current[time_temp][i][j];
                // avg_buff_temp += 1;

            }
        }
        vec_average.push_back(avg_temp/(divide));
        vec_buff_average.push_back(avg_buff_temp/(divide));
        vec_ubuff_average.push_back(avg_ubuff_temp/(divide));
        vec_Piezo_average.push_back(avg_piezo_temp/(divide));
        //cout << avg_temp/(x_max*y_max) << endl;
        avg_temp = 0;
        avg_buff_temp = 0;
        avg_ubuff_temp = 0;
        avg_piezo_temp = 0;
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

    for(int time = 0; time <= time_max_calc; time++){
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
    sizes.insert(sizes.begin(),vec_buff_average.size());
    sizes.insert(sizes.begin(),vec_ubuff_average.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector length: " << max_size << endl;

    bool bool_average;
    bool bool_buff_average;
    bool bool_ubuff_average;
    bool bool_Piezo_average; 

    //cout << "Break point 4" << endl;

    myfile << "Average,buff_average,ubuff_average,Piezo_avg\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        bool_average = (vec_average.size() > i) ? true : false;
        bool_buff_average = (vec_buff_average.size() > i) ? true : false;
        bool_ubuff_average = (vec_ubuff_average.size() > i) ? true : false;
        bool_Piezo_average = (vec_Piezo_average.size() > i) ? true : false;

        //cout << "Break point 6" << endl;
        
        if(bool_average) myfile << vec_average[i] << ",";
        if(!bool_average) myfile << ",";
        if(bool_buff_average) myfile << vec_buff_average[i] << ",";
        if(!bool_buff_average) myfile << ",";
        if(bool_ubuff_average) myfile << vec_ubuff_average[i] << ",";
        if(!bool_ubuff_average) myfile << ",";
        if(bool_Piezo_average) myfile << vec_Piezo_average[i];
        

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
