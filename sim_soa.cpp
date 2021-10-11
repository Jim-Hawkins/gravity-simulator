//
// Created by mariwogr on 10/10/21.
//

#include <string>
#include <iostream>
#include <fstream>

#include <random>
#include <cmath>

#include "sim_soa.hpp"

using namespace std;

double pos_generator(int random_seed, float size_enclosure) {

    // Motor
    std::mt19937_64 gen64(random_seed);
    // Distribución uniforme
    std::uniform_real_distribution<> dis(0, size_enclosure);

    // Posición
    double pos = dis(gen64);

    return pos;
}

double weight_generator(int random_seed) {

    // Motor
    std::mt19937_64 gen64(random_seed);
    // Normal distribution
    std::normal_distribution<> d{pow(10, 21), pow(10, 15)};

    // Weight
    double weight = d(gen64);

    return weight;
}

/* *
 * This function will check the parameters at the beginning
 *
 * @param int id                 whether it's the first or last file
 * @param parameters system_data data of the systema (size_enclosure, etc.)
 * @param set objects            structure containing the information of the objects
 * @return 0 on success
 */

int parser(int argc, char* argv[]){

    if (argc != 5 ){
        return -1;
    }
    //checking if the number of objs is smaller than cero
    if ( ( (int) *argv[1] ) < 0){
        return -2;
    }

    //checking if the number of iterations is smaller than cero
    if (((int)*argv[2])<0){
        return -2;
    }

    //checking the seed if it's a positive number
    if (((int)*argv[3])<0){
        return -2;
    }

    //checking if size_enclosure is positive
    if (((float)*argv[4])<0){
        return -2;
    }

    //checking if time_step is a real number positive
    if (((float)*argv[5])<0){
        return -2;
    }

    return 0;

}

double * gravitational_force_calc(set objects, int i, int j) {
    double G = 6.674 * pow(10, -11);
    double fuerza[3];

    double powSqX  = pow((objects.x[i] - objects.x[j]), 2);
    double powSqY  = pow((objects.y[i] - objects.y[j]), 2);
    double powSqZ  = pow((objects.z[i] - objects.z[j]), 2);

    fuerza[0] = (G * objects.m[i] * objects.m[j] * (objects.x[i] - objects.x[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));
    fuerza[1] = (G * objects.m[i] * objects.m[j] * (objects.y[i] - objects.y[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));
    fuerza[2] = (G * objects.m[i] * objects.m[j] * (objects.z[i] - objects.z[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));

    return fuerza;
}

double accel_calc(double m, double F) {
    return (1/(m))*F;
}

int gravitational_force(int i, int j) {

    int grav_force;
    int fuerza[3] = {0,0,0};
    for(int i = 0; i < num_objects; i++){
        for(int j = 0; j < num_objects; j++){
            if (i != j){
                fuerza[0] += gravitational_force_calc(objects, i, j)[0];
                fuerza[1] += gravitational_force_calc(objects, i, j)[1];
                fuerza[2] += gravitational_force_calc(objects, i, j)[2];
            }
        }
        /*
        for (int j = 0; j < i; j++) {
            fuerza -= calculitos();
        }*/
        accel;
        vel;
        pos;
    }
}

int collision_objects(set object1, set object2){

    double mt;
    mt = object1.m[0] + object2.m[0];

}


0->           1 ->          + 0 con 1      + 1 con 2       + 2 con 3       - 3 con 0
| \          /|             + 0 con 2      + 1 con 3       - 2 con 0       - 3 con 1
V             V             + 0 con 3      - 1 con 0       - 2 con 1       - 3 con 2


2             3




int print_error_args(int argc, char* argv[]){
    /*This function will print in the standard output the parameters when the function was called
      and it will show the errors while doing it.*/
    cerr << argv[0] << " invoked with " << argc << " parameters." << endl;
    cerr << "Arguments:" << endl;

    /*If the argc is not 5, it will show a ? character when a parameter is not in the function call*/
    if (1 <= argc) { cerr << "  num_objects: " << argv[1] << endl; }
    else { cerr << "  num_objects: ?" << endl; }

    if (2 <= argc) { cerr << "  num_iterations: " << argv[2] << endl; }
    else { cerr << "  num_iterations: ?" << endl; }

    if (3 <= argc) { cerr << "  random_seed: " << argv[2] << endl; }
    else { cerr << "  random_seed: ?" << endl; }

    if (4 <= argc) { cerr << "  size_enclosure: " << argv[2] << endl; }
    else { cerr << "  size_enclosure: ?" << endl; }

    if (5 <= argc) { cerr << "  time_step: " << argv[2] << endl; }
    else { cerr << "  time_step: ?" << endl; }

    return 0;
}

int write_config(int id, parameters params){

    /*This function will write the parameters from the main program in the init_config file or in the final_config file*/
    ofstream out_file;

    /*If the id is 0 it will write the content in the init_config file*/
    if (id == 0){ out_file.open("init_config.txt"); }
    /*If the id is different than 0 the content will be written in the final_config file*/
    else { out_file.open("final_config.txt"); }
    out_file << "hasta luego";
    out_file.close();
    return 0;
}

int main(int argc, char* argv[]){
    /*The array of parameters argv passes through a parser to check all the arguments are correct*/
    int retcode = parser(argc, argv);

    /*The result of the parser will be equal to -1 or -2 if there are errors with the arguments*/
    if(retcode < 0){
        /*If there is an error, the program will call the function print_error_args to print them
         * through the standard output*/
        print_error_args(argc, argv);

        /*The main function will return retcode, which can be equal to -1 if there are not enough
         * arguments and -2 if the arguments are not correct*/
        return retcode;
    }

    parameters system_data{ (int) *argv[1], (int) *argv[2],
                            (int) *argv[3], (float) *argv[4],
                            (float) *argv[5]};

    set objects;

    for(int i = 0; i < system_data.num_objects; i++){
        objects.x[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.y[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.z[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.m[i] = weight_generator(system_data.random_seed);
    }

    write_config(0, system_data);


    return 0;
}
