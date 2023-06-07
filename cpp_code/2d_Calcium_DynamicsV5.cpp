#include "2d_Calcium_Dynamics.h"

const char *path1="../data_files/2d_Piezo_Channel.csv";
const char *path2="../data_files/2d_Piezo_Channel_avg.csv";

default_random_engine generator; //important: must be outside of the loop/method that calls it

normal_distribution<double> stiffness(0.7,0.1);
normal_distribution<double> pressure(30,5);
normal_distribution<double> voltage(-70,10);

MatrixXd Pairwise_distances(MatrixXd A, MatrixXd B){
    MatrixXd Distances(A.cols(),B.cols());
    for(int i = 0; i < A.cols(); i++){
        for(int j = 0; j < B.cols(); j++){
            Distances(i,j) = sqrt((A.col(i) - B.col(j)).array().pow(2).sum());
        }
    }
    return(Distances);
}

VectorXd vec_std_to_Eigen(vector<double> convert_me){
    int size = convert_me.size(); 
    VectorXd A(size);
    for(int i = 0; i < size; i++){
        A[i] = convert_me[i];
    }
    return(A);
}

MatrixXd mat_std_to_Eigen(vector<vector<double>> convert_me){
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

vector<double> vec_Eigen_to_std(VectorXd convert_me){ 
    int size = convert_me.size();
    vector<double> A(size);
    for(int i = 0; i < size; i++){
        A[i] = convert_me[i];
    }
    return(A);
}

vector<vector<double>> mat_Eigen_to_std(MatrixXd convert_me){
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

        for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if((i < size - 1) && (j < size - 1) && (i > 0) && (j > 0)){
                if(i == j){
                    euler_A_temp(i,j) = -2;
                    euler_A_temp(i,j + 1) = 1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
            if((i == 0)){
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j + 1) = 1;
                }
            }
            if((i == size - 1)){
                if(i == j){
                    euler_A_temp(i,j) = -1;
                    euler_A_temp(i,j - 1) = 1;
                }
            }
        }
    }

    return(euler_A_temp);
}

VectorXd backward_euler(VectorXd Diffuse_me, MatrixXd euler_A){
    int size = Diffuse_me.size();

    euler_A = (delta_T / pow(size_scale,2))*euler_A;
    VectorXd F2 = VectorXd::Zero(size);

    MatrixXd I(size, size);
    I.setIdentity(size, size);

    F2 = (I - euler_A).inverse()*Diffuse_me;

    return(F2);
}

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
    return (E);
}

vector<vector<double>> Compute_J_diffusion(vector<vector<double>> A){
    MatrixXd t1; 
    vector<vector<double>> B;
    bool euler_A_made; 
    bool euler_B_made;

    t1 = mat_std_to_Eigen(A);

    MatrixXd thalfway(t1.rows(),t1.cols());
    MatrixXd t2(t1.rows(),t1.cols()); 

    if(!euler_A_made){
        A_cols = euler_A_maker(t1.rows());
        euler_A_made = true; 
    }

    if(!euler_B_made){
        A_rows = euler_A_maker(t1.cols());
        euler_B_made = true; 
    }
    
    for(int i = 0; i < t1.cols(); i++){
        thalfway.col(i) = backward_euler(t1.col(i), A_cols);
    }

    for(int i = 0; i < thalfway.rows(); i++){
        t2.row(i) = backward_euler(thalfway.row(i), A_rows);
    }

    B = mat_Eigen_to_std(t2);

    return(B);
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

double J_Piezo(int time_counter, double Pressure_input, double Substrate_input, double Voltage_input){
    double p;

    if(time_counter == 0){
        vec_open1.push_back(0);
        vec_inactive_held.push_back(0);
        vec_inactive.push_back(0);
        vec_closed.push_back(1);
    }

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);

    double P_total = P_P*P_S + P_V;

    vec_P_Total.push_back(P_total);

    double P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1) - 0.00669;

    double tau_inact = 0.999;
    double tau_open = 0.9; 
        
    vec_open1.push_back(tau_open*vec_open1[time_counter] + (1 - exp(-delta_T*P_opening_temp))*vec_closed[time_counter]);
    vec_inactive_held.push_back(((1 - exp(-delta_T*P_P))*vec_inactive_held[time_counter]) + (P_P*vec_open1[time_counter])*(1-tau_open));
    vec_inactive.push_back(tau_inact*vec_inactive[time_counter] + ((1 - P_P)*(1-tau_open)*vec_open1[time_counter]) + ((exp(-delta_T*P_P))*vec_inactive_held[time_counter]));
    vec_closed.push_back((exp(-delta_T*P_opening_temp))*vec_closed[time_counter] + (vec_inactive[time_counter] - tau_inact*vec_inactive[time_counter]));

    p = vec_open1[time_counter] + vec_inactive_held[time_counter];

    return(p);
}

double Compute_efflux(double C_cyt, int loc){
    double efflux;

    if(C_cyt > mols_divs){
        efflux = 0.5*loc*(C_cyt - mols_divs);
    }
    else{
        efflux = 0;
    }

    return(efflux);
}

double Compute_J_on(double C_cyt){
    k_buff_bind = 0.600; //for the buffer BAPTA in mM
    k_buff_unbind = 0.100;

    J_on = (C_cyt/mols_divs)*k_buff_bind*C_cyt*buff_unbound;
    J_off = (mols_divs/C_cyt)*k_buff_unbind*buff_bound; 

    buff_diff = J_off - J_on;

    return(buff_diff*0.000001); //*0.000001
}

double main2(int x){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = mols_divs;
            vec_num_closed[0][i][j] = N_Piezo_channels;
            vec_num_open[0][i][j] = 0;
            vec_Piezo_current[0][i][j] = 0;
            vec_buff_bound[0][i][j] = buff_bound;
            vec_buff_unbound[0][i][j] = buff_unbound;
        }
    }
    vec_average.push_back(mols_divs);

    for(int time_temp = 0; time_temp <= time_max_calc; time_temp++){
        cout << time_temp << endl;

        double P_Piezo;

        if(time_temp % 100 < 9 && time_temp > 50){
            P_Piezo = J_Piezo(time_temp, 60, 0.5, -70);
        }
        else{
            P_Piezo = J_Piezo(time_temp, 10, 0.5, -70);
        }

        double avg_temp = 0; 

        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                int location = 0; 
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
                    location = 0;
                }

                double P_Piezo_temp;
                if(location != 0){
                    P_Piezo_temp = P_Piezo;
                }
                else{
                    P_Piezo_temp = 0; 
                }

                Piezo_current = P_Piezo_temp*G_Piezo_single*(N_Piezo_channels/(2*x_max + 2*y_max));

                vec_time[1][i][j] = vec_time[0][i][j] + Compute_J_on(vec_time[0][i][j]) + Piezo_current - Compute_efflux(vec_time[0][i][j], location);
                
                avg_temp += vec_time[1][i][j];
            }
        }

        avg_temp = avg_temp/(x_max * y_max);

        vec_average.push_back(avg_temp);

        vec_time[1] = Compute_J_diffusion(vec_time[1]);

        for (int x = 0; x <= x_max; x++){
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

        vec_time[0] = vec_time[1];
    }

    return(0);
}

double output_avg(double x)
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

    myfile << "Average\n";

    for (int i = 0; i < max_size; i++)
    {
        bool_average = (vec_average.size() > i) ? true : false;
        
        if(bool_average) myfile << vec_average[i];

        myfile << "\n";
    }

    myfile.close();

    return (0);
}

int main(void) {
    cout << "Begin" << endl;

    main2(0);
    output_avg(0);

    cout << "End" << endl;
}
