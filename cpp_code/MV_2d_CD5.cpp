#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// Global Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Growth cone dimensions %%%%%%%%%%%%%%%%%%
int x_max = 15; // Number of divisions in the X direction
int y_max = 15; // Number of divisions in the Y direction
double cell_size = 0.000003; // 3um 

double divs = (x_max + 1) * (y_max + 1); // Finds number of divisions

double delta_X = cell_size/x_max; // Finds delta_X in um

double time_max = 500; //in s
double delta_T = 2; // in s
double time_max_calc = time_max/delta_T; // Finds total # of time steps

// vec_time is a 3d vector, with the coordinates vec_time[time][x][y]
vector<vector<vector<double> > > vec_time(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1))); 
vector<double> vec_average; // Holds the average concentration of the growth cone
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double log_convert = 2.303; // to convert from ln to log10
double R_constant = 8.1345;
double F = 96485.3321;
double body_temp = 310.15;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Piezo Initializaiton %%%%%%%%%%%%%%%%%%%%
vector<double> vec_open1; // first open state
vector<double> vec_inactive; // inactive state
vector<double> vec_inactive_held; // held open state
vector<double> vec_closed; // closed state

// Eigen initialization
MatrixXd A_cols, A_rows; // Both are the ... 1, -2, 1, ... matrices for rows & cols

default_random_engine generator; // Random generator used to add stochasticity
normal_distribution<double> stiffness(0.5,0.03);
normal_distribution<double> pressure(0,1);
normal_distribution<double> voltage(-70,10);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MatrixXd mat_std_to_Eigen(vector<vector<double>> convert_me){
    // Converts from std class to Eigen class
    int rows = convert_me.size(); 
    int cols = convert_me[0].size(); 
    MatrixXd A(rows,cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            A(i,j) = convert_me[i][j];
        }
    }
    return(A);
}

vector<vector<double>> mat_Eigen_to_std(MatrixXd convert_me){
    // Converts from Eigen class to std class
    int rows = convert_me.rows(); 
    int cols = convert_me.cols(); 
    vector<double> B;
    vector<vector<double>> A;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            B.push_back(convert_me(i,j));
        }
        A.push_back(B);
        B.clear();
    }
    return(A);
}

MatrixXd euler_A_maker(int size){
    MatrixXd euler_A_temp = MatrixXd::Zero(size, size);

    // Fills the diagonal with ... 1, -2, 1 ...
    // If at a corner, fills with -1, 1, ... 
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if((i < size - 1) && (j < size - 1) && (i > 0) && (j > 0)){ // All along the diagonal
                if(i == j){
                    euler_A_temp(i,j) = -2;
                    euler_A_temp(i,j + 1) = 1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
            if((i == 0)){ // At top left corner
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j + 1) = 1;
                }
            }
            if((i == size - 1)){ // At bottom right corner
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
        }
    }

    return(euler_A_temp);
}

VectorXd backward_euler(VectorXd Diffuse_me, MatrixXd euler_A){ // Performs backward Euler on vectors inserted
    int size = Diffuse_me.size();

    euler_A = (delta_T / pow(delta_X,2))*euler_A;
    VectorXd F2 = VectorXd::Zero(size); // F2 is where our diffused vector will be saved/returned

    MatrixXd I(size, size);
    I.setIdentity(size, size);

    F2 = (I - euler_A).inverse()*Diffuse_me; 

    return(F2);
}

vector<vector<double>> Compute_J_diffusion(vector<vector<double>> A){ // Controller to solve diffusion
    MatrixXd t1; // Our un-diffused matrix
    vector<vector<double>> B;
    bool euler_A_made; // Bool to ensure we only initialize matrix A once
    bool euler_B_made; // Bool to ensure we only initialize matrix B once

    t1 = mat_std_to_Eigen(A); // Converts 2d vector (in std class) to MatrixXd (in Eigen class) 

    MatrixXd thalfway(t1.rows(),t1.cols());
    MatrixXd t2(t1.rows(),t1.cols()); 

    if(!euler_A_made){ // Makes ... 1, -2, 1, ... matrix A for rows
        A_cols = euler_A_maker(t1.rows());
        euler_A_made = true; 
    }

    if(!euler_B_made){ // Makes ... 1, -2, 1, ... matrix B for rows
        A_rows = euler_A_maker(t1.cols());
        euler_B_made = true; 
    }
    
    // Note: backward_euler() is called in order to run Backward Euler for each vector in rows, then in cols
    for(int i = 0; i < t1.cols(); i++){ // Diffusion along rows
        thalfway.col(i) = backward_euler(t1.col(i), A_cols);
    }

    for(int i = 0; i < thalfway.rows(); i++){ // Diffusion along cols
        t2.row(i) = backward_euler(thalfway.row(i), A_rows);
    }

    B = mat_Eigen_to_std(t2); // Converts MatrixXd (in Eigen class) to 2d vector (in std class)

    return(B);
}

double Probability_of_Piezo(int time_counter, double Pressure_input, double Substrate_input, double Voltage_input){
    double p; // local probability of Piezo, returned at the end

    if(time_counter == 0){ // initializes all of the Piezo channels to be closed at time(0)
        // These vectors are declared globally
        vec_open1.push_back(0);
        vec_inactive_held.push_back(0);
        vec_inactive.push_back(0);
        vec_closed.push_back(1);
    }

    double P_P = 1/(exp((30 - Pressure_input)/6) + 1); // probability of opening due to pressure 
    double P_S = (1/(0.25*pow(2*M_PI,0.5))*exp(-0.5*pow((Substrate_input - 0.7)/0.25,2)))/1.6; // probability of opening due to substrate stiffness
    double P_V = 1/(exp((100 - Voltage_input)/20) + 1); // probability of opening due to voltage 

    double P_total = P_P*P_S + P_V; // used as input for total probability 

    double P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1) - 0.00669; // total opening probability

    double tau_inact = 0.999; // time constant for inactive state
    double tau_open = 0.9; // time constant for open state (when not held open)
    
    // solves & tracks the Markov/state transition for Piezo in vectors
    vec_open1.push_back(tau_open*vec_open1[time_counter] + (1 - exp(-delta_T*P_opening_temp))*vec_closed[time_counter]);
    vec_inactive_held.push_back(((1 - exp(-delta_T*P_P))*vec_inactive_held[time_counter]) + (P_P*vec_open1[time_counter])*(1-tau_open));
    vec_inactive.push_back(tau_inact*vec_inactive[time_counter] + ((1 - P_P)*(1-tau_open)*vec_open1[time_counter]) + ((exp(-delta_T*P_P))*vec_inactive_held[time_counter]));
    vec_closed.push_back((exp(-delta_T*P_opening_temp))*vec_closed[time_counter] + (vec_inactive[time_counter] - tau_inact*vec_inactive[time_counter]));

    p = vec_open1[time_counter] + vec_inactive_held[time_counter]; // probability of having any open form

    return(p);
}

double Compute_J_efflux(double C_cyt, int loc){
    // efflux is calculated by magnitude of difference between edge piece and desired concentration of 120nm
    double efflux;

    if(C_cyt > 0.0000000012){
        efflux = 0.5*loc*(C_cyt - 0.0000000012); // loc is used to ensure efflux happens only at the edges
    }
    else{
        efflux = 0; //if not at edge, efflux is 0
    }

    return(efflux);
}

double Compute_J_on(double C_cyt){ // This method calculates buffering 
    double buff_unbound = 0.000001; // Concentration of unbound buffer, which we are taking to be b_total
    double buff_bound = 0.000007; // Concentration of bound buffer
    double J_on; // Function of binding of Ca2+ to buffers, Page 309 of mathematical physiology seems good
    double J_off; // Function of unbinding of Ca2+ from buffers
    double k_buff_bind = 0.600; // Binding affinity/Kon for the buffer BAPTA in mM
    double k_buff_unbind = 0.100; // Unbinding affinity/Koff of buffer

    J_on = delta_T*(C_cyt/0.0000000012)*k_buff_bind*C_cyt*buff_unbound;
    J_off = delta_T*(0.0000000012/C_cyt)*k_buff_unbind*buff_bound; //multiplying by delta_T to scale with time

    double buff_diff; //store difference between J_on and J_off
    buff_diff = J_off - J_on;

    return(buff_diff); 
}

double Write_Out_Average(double x){ // This method is used to output the values stored in the averages vector
    const char *path2="../data_files/2d_Piezo_Channel_avg.csv";
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2); // Defines the .csv file used to store averages

    myfile << "Average\n";

    // Writes out averages, starting at time(100) to allow the matrix to reach equilibrium
    for (int i = 100; i < vec_average.size(); i++){
        myfile << vec_average[i];
        myfile << "\n";
    }

    myfile.close(); // Closes file

    return (0);
}

double Model_Growth_Cone(int x){ // This method calls all other methods used in this program
    const char *path1="../data_files/2d_Piezo_Channel.csv";
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1); // Defines .csv file used to store growth cone matrix values

    random_device rd;
    mt19937 gen(rd()); // Random number generator used to decide if Piezo is open or closed

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = 0.00000012; // Fills 2d vector with concentration of 120nM at time(0)
        }
    }

    vec_average.push_back(0.00000012); // First value of average vector is 120nM

    int counter = 0; // counter counts up with time in integer values, i.e., rather than with delta_T it is 1 every time
    for(double time_temp = 0; time_temp <= time_max; time_temp+=delta_T){
        cout << time_temp << endl; // prints the time, to see progress in the code

        double P_Piezo; // this is the probability of Piezo being open, and is solved by calling Probability_of_Piezo()

        if(counter % 100 < 30){ // used to simulate force being applied at regular intervals
            P_Piezo = Probability_of_Piezo(counter, 60, stiffness(generator), voltage(generator));
        }
        else{
            P_Piezo = Probability_of_Piezo(counter, abs(pressure(generator)), abs(stiffness(generator)), voltage(generator));
        }

        double avg_temp = 0; // used to average the concentration in the growth cone

        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                int location = 0; // used to determine location in growth cone (either corner, edge, or center)
                if(((i == 0) || (i == x_max)) && ((j == 0) || (j == y_max))){
                    location = 1;
                }
                else if(((i == 0) || (i == x_max))){
                    location = 1;
                }
                else if((j == 0) || (j == y_max)){
                    location = 1;
                }
                else{
                    location = 0; // All center positions have value = 0
                }

                // Uses P_Piezo as probability to generate a 1 (representing open) or 0 (representing closed)
                // This is done for each edge position to quantize opening vs. closing as opposed to a decimal value
                discrete_distribution<> distrib({1 - P_Piezo, P_Piezo});
                int On_Off_rand = distrib(gen); 

                double Piezo_current;
                double G_Piezo_single = 0.000000000030; //Conductance of single channel
                int N_Piezo_channels = 100000; // Total number of Piezo channels in entire growth cone
                Piezo_current = On_Off_rand*location*G_Piezo_single*(N_Piezo_channels/(2*x_max + 2*y_max)); // Solves current generated Piezo

                // Total Ca2+ for each position:
                vec_time[1][i][j] = vec_time[0][i][j] + Compute_J_on(vec_time[0][i][j]) + Piezo_current - Compute_J_efflux(vec_time[0][i][j], location);
                
                avg_temp += vec_time[1][i][j]; // Sums the total Ca2+, which is then used to compute the average
            }
        }

        avg_temp = avg_temp/((x_max + 1) * (y_max + 1)); // Computes the average 

        vec_average.push_back(avg_temp); // Stores the averages

        vec_time[1] = Compute_J_diffusion(vec_time[1]); // Diffuses all the Ca2+ within the growth cone

        
        if(time_temp > 100){ // Writes out the growth cone matrix (vec_time) to .csv file
        // Writing out begins after time 100, which is in order to allow the matrix to reach equilibrium
            for (int x = 0; x <= x_max; x++){ // Iterates through rows & columns 
                for (int y = 0; y <= y_max; y++){  
                    if(y < y_max) {
                        myfile << vec_time[1][x][y] << ",";
                    }  
                    else {
                        myfile << vec_time[1][x][y];
                    }
                }
                myfile << "\n";
            }
            myfile << "\n";
        }

        vec_time[0] = vec_time[1]; // Sets the output of this code to be the input for the next recursion
        counter++; // Updates with each delta_T
    }

    myfile.close(); // Closes file

    Write_Out_Average(0); // Writes out average values of the growth cone

    return(0);
}

int main(void) {

    // NOTE: This code uses two classes of matrices/vectors. That is: the std:: class 
    // (which is built into C++) and the Eigen:: class (from the package Eigen), listed
    // in the includes. Because the code was previously built using the std class, I 
    // decided that rather than rewriting all of the code into the Eigen:: class, it 
    // was easiest to build mini-modules to convert between the two as needed, such as 
    // when calculating diffusion. These modules include 'mat_std_to_Eigen()', etc. 

    cout << "Begin" << endl;

    Model_Growth_Cone(0); // Runs all of the algorithms and writes growth cone matrix

    cout << "End" << endl;
}
