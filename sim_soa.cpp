//
// Created by mariwogr on 10/10/21.
//

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <cmath>

#include "sim_soa.hpp"

using namespace std;

double pos_generator(int random_seed, float size_enclosure) {

    // Motor
    std::mt19937_64 gen64(random_seed);

    // Uniform Distribution
    std::uniform_real_distribution<> dis(0, size_enclosure);

    // Position
    double position = dis(gen64);
    return position;
}

double mass_generator(int random_seed) {

    // Motor
    std::mt19937_64 gen64(random_seed);
    // Normal distribution
    std::normal_distribution<> d{pow(10, 21), pow(10, 15)};

    // Mass
    double mass = d(gen64);

    return mass;
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

/*
 * This function will calculate all the gravitational force components of i and j in the structure set objects
 *
 * @param: set objects           structure of objects, with all the components of each point in the simulator
 * @param: int i                 position of the first point
 * @param: int j                 position of the second point
 * @return force                 resulting force vector
 */
double * gravitational_force_calc(set objects, int i, int j) {
    double G = 6.674 * pow(10, -11);
    double force[3];

    double powSqX  = pow((objects.x[i] - objects.x[j]), 2);
    double powSqY  = pow((objects.y[i] - objects.y[j]), 2);
    double powSqZ  = pow((objects.z[i] - objects.z[j]), 2);

    force[0] = (G * objects.m[i] * objects.m[j] * (objects.x[i] - objects.x[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));
    force[1] = (G * objects.m[i] * objects.m[j] * (objects.y[i] - objects.y[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));
    force[2] = (G * objects.m[i] * objects.m[j] * (objects.z[i] - objects.z[j]))/(pow(sqrt(powSqX + powSqY + powSqZ),3));

    return force;
}

/*
 * This function will return a component of the acceleration of a point with the mass m given and the component of the sum of forces F of the point
 *
 * @param: double m             mass of the point
 * @param: double F             one component of the sum of forces F
 *
 * @return (1/(m))*F            the point acceleration point
 */
double accel_calc(double m, double F) {
    return (1/(m))*F;
}


/*
 * This function updates the speed vector v and the position of every point in the set objects of points
 *
 * @param: int num_objects          the total of points in the set of objects
 * @param: set objects              structure of objects, with all the components of each point in the simulator
 * @param: float time_step          time step to obtain the speed and position of the point
 *
 * @return 0                        if the function was executed correctly
 */
int gravitational_force(int num_objects, set objects, float time_step) {

    // The 3 components of the gravitational force will be set to 0
    double force[3] = {0,0,0};
    double accel[3] = {0,0,0};

    // The execution will pass through two nested loops to obtain the sum of gravitational forces of every point with the other points
    for(int i = 0; i < num_objects; i++) {
        for (int j = 0; j < num_objects; j++) {

            // First it checks that the two points are active (not collided)
            if (objects.active[i] && objects.active[j]) {

                // If the two points are not the same, it will sum the force of every component to the total force
                if (i != j) {
                    force[0] += gravitational_force_calc(objects, i, j)[0];
                    force[1] += gravitational_force_calc(objects, i, j)[1];
                    force[2] += gravitational_force_calc(objects, i, j)[2];
                }

        // Update the acceleration
        accel[0] = accel_calc(objects.m[i], force[0]);
        accel[1] = accel_calc(objects.m[i], force[1]);
        accel[2] = accel_calc(objects.m[i], force[2]);

        // Update the speed
        objects.vx[i] = objects.vx[i] + accel[0] * time_step;
        objects.vy[i] = objects.vy[i] + accel[1] * time_step;
        objects.vz[i] = objects.vz[i] + accel[2] * time_step;

        // Update the position
        objects.x[i] = objects.x[i] + objects.vx[i] * time_step;
        objects.y[i] = objects.y[i] + objects.vy[i] * time_step;
        objects.z[i] = objects.z[i] + objects.vz[i] * time_step;
            }
        }
    }
    return 0;
}

int check_bounce(set objects, int obj){

    if(objects.x[obj]<=0){
        objects.x[obj] = 0;
        objects.vx[obj] = -1 * objects.vx[obj];
    }

    if(objects.y[obj]<=0){
        objects.y[obj] = 0;
        objects.vy[obj] = -1 * objects.vy[obj];
    }

    if(objects.z[obj]<=0){
        objects.z[obj] = 0;
        objects.vz[obj] = -1 * objects.vz[obj];
    }



    return 0;

}

int collision_objects(set objects,int i, int j){

    if (objects.active[i] && objects.active[j]){

        objects.m[i] = objects.m[i] + objects.m[j];
        objects.vx[i] = objects.vx[i] + objects.vx[j];
        objects.vy[i] = objects.vy[i] + objects.vy[j];
        objects.vz[i] = objects.vz[i] + objects.vz[j];

        objects.active[j] = false;
    }
}


/*                  DO NOT TOUCH THIS, DEATH PENALTY
0->           1 ->          + 0 con 1      + 1 con 2       + 2 con 3       - 3 con 0
| \          /|             + 0 con 2      + 1 con 3       - 2 con 0       - 3 con 1
V             V             + 0 con 3      - 1 con 0       - 2 con 1       - 3 con 2


2             3

*/


/* *
 * This function will write the errors in the parameter in error case
 *
 * @param int argc                 This is the number of arguments
 * @param char* argv is a pointer to array of chars (strings) which it has the values
 *
 * @return 0 on success
 */
int print_error_args(int argc, char* argv[]) {
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

/* *
 * This function will write the parameters from the main program in the init_config file or in the final_config file
 *
 * @param int id                 whether it's the first or last file
 * @param parameters system_data data of the systema (size_enclosure, etc.)
 * @param set objects            structure containing the information of the objects
 * @return 0 on success
 */
int write_config(int id, parameters system_data, set objects){
    ofstream out_file;
    char res[50];

    /*If the id is 0 it will write the content in the init_config file*/
    if (id == 0){ out_file.open("init_config.txt"); }
    /*If the id is different than 0 the content will be written in the final_config file*/
    else { out_file.open("final_config.txt"); }

    sprintf(res, "%.3f", system_data.size_enclosure);
    out_file << res;
    sprintf(res, "%.3f", system_data.time_step);
    out_file << res;
    sprintf(res, "%d.000", system_data.num_objects);
    out_file << res << endl;

    for(int i = 0; i < system_data.num_objects; i++){
        sprintf(res,
                "%.3f %.3f %.3f %.3f %.3f %.3f %.3f",
                objects.x[i], objects.y[i], objects.z[i], objects.vx[i], objects.vy[i], objects.vz[i], objects.m[i]);
        out_file << res << endl;
    }
    out_file.close();
    return 0;
}

int main(int argc, char* argv[]) {
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

    /* Store simulation arguments in a structure */
    parameters system_data{ (int) *argv[1], (int) *argv[2],
                            (int) *argv[3], (float) *argv[4],
                            (float) *argv[5]};

    /* Declare the structure that holds objects' information */
    set objects;

    /* Initialize x, y, z and m attributes of each object */
    for(int i = 0; i < system_data.num_objects; i++){
        objects.x[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.y[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.z[i] = pos_generator(system_data.random_seed, system_data.size_enclosure);
        objects.m[i] = mass_generator(system_data.random_seed);
        objects.active[i] = true;
    }
    /* Write initial configuration to a file*/
    write_config(0, system_data, objects);

    /* Body of the simulation */
    for(int i = 0; i < system_data.num_iterations; i++){
        gravitational_force(system_data.num_objects, objects, system_data.time_step);
        for(int obj = 0; obj < system_data.num_objects; obj++){
            check_bounce(objects, obj, system_data.size_enclosure);
            for(int c = 0; c < system_data.num_objects; c++){
                if (objects){}
                check_collision(objects, obj, c);
            }
        }
    }

    /* Write final configuration to a file*/
    write_config(1, system_data, objects);
    return 0;
}
