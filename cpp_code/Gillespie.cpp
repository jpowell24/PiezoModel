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
int time_max_calc = time_max/delta_T; // Finds total # of time steps

// vec_time is a 3d vector, with the coordinates vec_time[time][x][y]
vector<vector<vector<double> > > vec_time(time_max_calc + 2, vector<vector<double> >(y_max + 1, vector<double>(x_max + 1))); 
vector<double> vec_average; // Holds the average concentration of the growth cone
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Gill initialization %%%%%%%%%%%%%%%%%%%%%
int N = 1000; 
vector<vector<double> > Gillespie_Datatype(N, vector<double>(4));

// Data type key: 
// 0         1          2         3          4        
 // index     x          y         state      time

// States: 
// 0         1          2         3          
// Open1     Open2      Inact     Closed

int dimensions = 2*(x_max + 1) + 2*(y_max + 1);
int scale = dimensions-4;
vector<vector<double> > positions_key(scale, vector<double>(2));
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double initialize_Gill(int x){
    vector<double> temp; 

    int counter = 0; 

    for(int i = 0; i <= y_max-1; i++){
        temp.push_back(0);
        temp.push_back(i);
        positions_key[counter] = temp;
        counter++;
        temp.clear();
    }
    for(int i = 0; i <= y_max; i++){
        temp.push_back(x_max);
        temp.push_back(i);
        positions_key[counter] = temp;
        counter++;
        temp.clear();
    }
    for(int i = 1; i <= x_max-1; i++){ // x_max - 1 prevents double counting of corners
        temp.push_back(i);
        temp.push_back(0);
        positions_key[counter] = temp;
        counter++;
        temp.clear();
    }
    for(int i = 0; i <= x_max-1; i++){
        temp.push_back(i);
        temp.push_back(y_max);
        positions_key[counter] = temp;
        counter++;
        temp.clear();
    }
  
    for(int i = 0; i < Gillespie_Datatype.size(); i++){
        int temp = rand() % scale;
        Gillespie_Datatype[i][0] = i; 
        Gillespie_Datatype[i][1] = positions_key[temp][0];
        Gillespie_Datatype[i][2] = positions_key[temp][1];
        Gillespie_Datatype[i][3] = 0;
        Gillespie_Datatype[i][4] = 0;
    }

    for(int i = 0; i < Gillespie_Datatype.size(); i++){
        double P_x = 0.3; // Probability
        double temp_rand = (float) rand()/RAND_MAX; // Generate random number
        double temp_scaled_rand = -pow(P_x,-1)*log(temp_rand); // Generate random time
        Gillespie_Datatype[i][4] = round(temp_scaled_rand/delta_T)*delta_T; // Scale random time to nearest delta_T
    }

    return(0);
}

double Gillespie(int x){ 
    const char *path1="../data_files/Gillespie.csv";
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1); // Defines .csv file used to store growth cone matrix values

    srand(time(0)); // Ensures the random number generator is seeded differently each time

    initialize_Gill(0);

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = 0; // Fills 2d vector with 0
        }
    }


    for(int time = 0; time < 20; time+=delta_T){

        for(int i = 0; i < Gillespie_Datatype.size(); i++){
            
            if(Gillespie_Datatype[i][3] == 0 || Gillespie_Datatype[i][3] == 1){
                vec_time[0][Gillespie_Datatype[i][1]][Gillespie_Datatype[i][2]] = vec_time[0][Gillespie_Datatype[i][1]][Gillespie_Datatype[i][2]] + 1; 
            }

            if(Gillespie_Datatype[i][4] < delta_T){
                double P_x;
                if(Gillespie_Datatype[i][3] == 0){
                    P_x = 0.3;
                    double P_1 = 0.2; 
                    double P_2 = 0.8; 
                    double controller = (float) rand()/RAND_MAX;

                    if(controller > P_1){
                        Gillespie_Datatype[i][3] = 2; 
                    }
                    else{
                        Gillespie_Datatype[i][3] = 1; 
                    }
                }   
                else if(Gillespie_Datatype[i][3] == 1){
                    P_x = 0.3;
                    Gillespie_Datatype[i][3] = 2; 
                }   
                else if(Gillespie_Datatype[i][3] == 2){
                    P_x = 0.3;
                    Gillespie_Datatype[i][3] = 3; 
                }   
                else if(Gillespie_Datatype[i][3] == 3){
                    P_x = 0.3;
                    Gillespie_Datatype[i][3] = 0; 
                }   
                else{
                    cout << "We have a serious problem" << endl;
                }

                double temp_rand = (float) rand()/RAND_MAX; // Generate random number
                double temp_scaled_rand = -pow(P_x,-1)*log(temp_rand); // Generate random time
                Gillespie_Datatype[i][4] = round(temp_scaled_rand/delta_T)*delta_T; // Scale random time to nearest delta_T
            }

            Gillespie_Datatype[i][4] =  Gillespie_Datatype[i][4] - delta_T; 


        }

        // for (int x = 0; x < Gillespie_Datatype.size(); x++){ // Iterates through rows & columns 
        //     for (int y = 0; y <= Gillespie_Datatype[x].size(); y++){
        //         if(y == Gillespie_Datatype[x].size()){
        //             cout << Gillespie_Datatype[x][y] << endl;  
        //         }
        //         else{
        //         cout << Gillespie_Datatype[x][y] << ",";
        //         }  
        //     }
        // }
        // cout << endl;

    }

    for (int x = 0; x <= x_max; x++){ // Iterates through rows & columns 
        for (int y = 0; y <= y_max; y++){  
            if(y < y_max) {
                myfile << vec_time[0][x][y] << ",";
            }  
            else {
                myfile << vec_time[0][x][y];
            }
        }
        myfile << "\n";
    }
    myfile << "\n";

    // for (int x = 0; x < Gillespie_Datatype.size(); x++){ // Iterates through rows & columns 
    //     vec_time[0][Gillespie_Datatype[x][1]][Gillespie_Datatype[x][2]] = vec_time[0][Gillespie_Datatype[x][1]][Gillespie_Datatype[x][2]] + 1; 
    // }

    // for (int x = 0; x <= x_max; x++){ // Iterates through rows & columns 
    //     for (int y = 0; y <= y_max; y++){  
    //         if(y < y_max) {
    //             myfile << vec_time[0][x][y] << ",";
    //         }  
    //         else {
    //             myfile << vec_time[0][x][y];
    //         }
    //     }
    //     myfile << "\n";
    // }
    // myfile << "\n";

    // for (int x = 0; x < positions_key.size(); x++){ // Iterates through rows & columns 
    //     for (int y = 0; y < positions_key[x].size(); y++){
    //         if(y == positions_key[x].size()-1){
    //             cout << positions_key[x][y] << endl;  
    //         }
    //         else{
    //         cout << x << " = " << positions_key[x][y] << ",";
    //         }  
    //     }
    //  }
    // cout << endl;

    // for(int i = 0; i < 10; i++){
    //     int index_temp = Assign_Positions(scale);
    //     cout << index_temp << " = " << positions_key[index_temp][0] << "," << positions_key[index_temp][1] << endl;

    // }


    // for (int x = 0; x <= x_max; x++){ // Iterates through rows & columns 
    //     for (int y = 0; y <= y_max; y++){  
    //         if(y < y_max) {
    //             myfile << vec_time[1][x][y] << ",";
    //         }  
    //         else {
    //             myfile << vec_time[1][x][y];
    //         }
    //     }
    //     myfile << "\n";
    // }
    // myfile << "\n";

    myfile.close();
    return(0);
}

int main(void) {
    cout << "Begin" << endl;

    Gillespie(0);

    cout << "End" << endl;
}
