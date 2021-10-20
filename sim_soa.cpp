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

/* *
 * This function will check the parameters at the beginning
 *
 * @param int argc                 its the number about how many parameters are in the execution
 * @param char argv             it's an array of chars, inside it, we have the arguments
 * @return 0 on success
 */

int parser(int argc, char* argv[]){

    if (argc != 6 ){
        return -1;
    }
    //checking if the number of objs is smaller than zero
    if ( ( (int) *argv[1] ) < 0){
        return -2;
    }

    //checking if the number of iterations is smaller than zero
    if ( ( (int) *argv[2] ) < 0){
        return -2;
    }

    //checking the seed if it's a positive number
    if ( ( (int) *argv[3] ) < 0){
        return -2;
    }

    //checking if size_enclosure is positive
    if ( ( (double) *argv[4] ) < 0){
        return -2;
    }

    //checking if time_step is a real number positive
    if ( ( (double) *argv[5] ) < 0){
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
void gravitational_force_calc(set objects, int i, int j, double *force) {
    double G = 6.674 * 1E-11;

    double powSqX  = pow((objects.x[i] - objects.x[j]), 2);
    double powSqY  = pow((objects.y[i] - objects.y[j]), 2);
    double powSqZ  = pow((objects.z[i] - objects.z[j]), 2);
    double norm = std::sqrt(powSqX + powSqY + powSqZ);

    // It will return the three components of the gravitational force between i and j
    force[0] -= (G * objects.m[i] * objects.m[j] * (objects.x[i] - objects.x[j]))/(norm * norm * norm);
    force[1] -= (G * objects.m[i] * objects.m[j] * (objects.y[i] - objects.y[j]))/(norm * norm * norm);
    force[2] -= (G * objects.m[i] * objects.m[j] * (objects.z[i] - objects.z[j]))/(norm * norm * norm);
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
    return (1/m)*F;
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
int gravitational_force(int num_objects, set objects, double time_step, double *force, double *accel) {
    // The execution will pass through two nested loops to obtain the sum of gravitational forces of every point with the other points
    for(int i = 0; i < num_objects; i++) {
        if (!objects.active[i]) { continue; }
        for (int j = 0; j < num_objects; j++) {

            // First it checks that the two points are active (not collided). If the two points are not the same,
            // it will sum the force of every component to the total force
            if (objects.active[j] && i != j) {
                gravitational_force_calc(objects, i, j, &force[3 * i]);
            }
        }
    }
    // Once we have a screenshot of the system in force array, update each active object
    for (int i = 0; i < num_objects; i++) {
        if(objects.active[i]) {
            // Updates the acceleration
            accel[0] = accel_calc(objects.m[i], force[i * 3]);
            accel[1] = accel_calc(objects.m[i], force[(i * 3) + 1]);
            accel[2] = accel_calc(objects.m[i], force[(i * 3) + 2]);

            // Updates the speed
            objects.vx[i] = objects.vx[i] + accel[0] * time_step;
            objects.vy[i] = objects.vy[i] + accel[1] * time_step;
            objects.vz[i] = objects.vz[i] + accel[2] * time_step;

            // Updates the position
            objects.x[i] = objects.x[i] + objects.vx[i] * time_step;
            objects.y[i] = objects.y[i] + objects.vy[i] * time_step;
            objects.z[i] = objects.z[i] + objects.vz[i] * time_step;
        }
    }
    return 0;
}

/*
 * This function checks if the object bounce with a wall and change the values if it's necessary
 *
 * @param: set objects              structure of objects, with all the components of each point in the simulator
 * @param: float size         It's the size of the wall
 * @param: int obj          it's the object which we are going to check
 *
 * @return 0                        if the function was executed correctly
 */
int check_bounce(set objects, int obj, double size){

    //check if the object bounce with a wall

    if(objects.x[obj] <= 0){
        objects.x[obj] = 0;
        objects.vx[obj] = -1 * objects.vx[obj];
    }

    if(objects.y[obj] <= 0){
        objects.y[obj] = 0;
        objects.vy[obj] = -1 * objects.vy[obj];
    }

    if(objects.z[obj] <= 0){
        objects.z[obj] = 0;
        objects.vz[obj] = -1 * objects.vz[obj];
    }

    if(objects.x[obj] >= size){
        objects.x[obj] = size;
        objects.vx[obj] = -1 * objects.vx[obj];
    }

    if(objects.y[obj] >= size){
        objects.y[obj] = size;
        objects.vy[obj] = -1 * objects.vy[obj];
    }

    if(objects.z[obj] >= size){
        objects.z[obj] = size;
        objects.vz[obj] = -1 * objects.vz[obj];
    }

    return 0;
}

/*
 * This function will calculate the euclidean distance between two points i and j in the set of objects: objects (structure
 * of arrays)
 *
 * @param: set objects              structure of objects, with all the components of each point in the simulator
 * @param: int i                    i-position of one point in the set of objects
 * @param: int j                    j-position of one point in the set of objects
 *
 * return 0                         if the execution was correctly executed
 */
int check_collision(set objects, int i, int j){
    double distance = std::sqrt(pow((objects.x[i] - objects.x[j]), 2)\
                            + pow((objects.y[i] - objects.y[j]), 2)\
                            + pow((objects.z[i] - objects.z[j]), 2));

    if(distance < 1.0){
        collision_objects(objects, i, j);
    }
    return 0;
}


/*
* This function will update the objects and their collisions
*
* @param: set objects          array of objects with their properties
* @param: int i                array position of the first object
* @param: int j                array position of the second object
*/
int collision_objects(set objects, int i, int j){

    objects.m[i] = objects.m[i] + objects.m[j];
    objects.vx[i] = objects.vx[i] + objects.vx[j];
    objects.vy[i] = objects.vy[i] + objects.vy[j];
    objects.vz[i] = objects.vz[i] + objects.vz[j];

    objects.active[j] = false;

    return 0;
}

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
    cerr << argv[0] << " invoked with " << argc - 1  << " parameters." << endl;
    cerr << "Arguments:" << endl;

    /* If argc is not 5, it will show a ? character when a parameter is not in the function call */
    if (1 < argc) { cerr << "  num_objects: " << argv[1] << endl; }
    else { cerr << "  num_objects: ?" << endl; }

    if (2 <= argc) { cerr << "  num_iterations: " << argv[2] << endl; }
    else { cerr << "  num_iterations: ?" << endl; }

    if (3 <= argc) { cerr << "  random_seed: " << argv[3] << endl; }
    else { cerr << "  random_seed: ?" << endl; }

    if (4 <= argc) { cerr << "  size_enclosure: " << argv[4] << endl; }
    else { cerr << "  size_enclosure: ?" << endl; }

    if (5 <= argc) { cerr << "  time_step: " << argv[5] << endl; }
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
    char res[100];

    /*If the id is 0 it will write the content in the init_config file*/
    if (id == 0){
        out_file.open("../init_config.txt");
    }
    /*If the id is different from 0 the content will be written in the final_config file*/
    else { out_file.open("../final_config.txt"); }

    sprintf(res, "%.3f ", system_data.size_enclosure);
    out_file << res;
    sprintf(res, "%.3f ", system_data.time_step);
    out_file << res;
    sprintf(res, "%d", system_data.num_objects);
    out_file << res << endl;

    for(int i = 0; i < system_data.num_objects; i++){
        if(objects.active[i]) {
            sprintf(res,
                    "%.3f %.3f %.3f %.3f %.3f %.3f %.3f",
                    objects.x[i], objects.y[i], objects.z[i], objects.vx[i], objects.vy[i], objects.vz[i],
                    objects.m[i]);
            out_file << res << endl;
        }

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
    parameters system_data{ stoi(argv[1]), stoi(argv[2]),
                            stoi(argv[3]), stod(argv[4]),
                            stod(argv[5])};

    /* Declare the structure that holds objects' information */
    set objects{
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (double *) malloc(sizeof(double) * system_data.num_objects),
        (bool *) malloc(sizeof(bool) * system_data.num_objects)
    };

    cout << "Creating simulation:" << endl;
    cout << "  num_objects: " << system_data.num_objects << endl;
    cout << "  num_iterations: " << system_data.num_iterations << endl;
    cout << "  random_seed: " << system_data.random_seed << endl;
    cout << "  size_enclosure: " << system_data.size_enclosure << endl;
    cout << "  time_step: " << system_data.time_step << endl;


	/* Create mersenne-twister generator and create a uniform and a normal distribution */
    mt19937_64 gen64(system_data.random_seed);
    uniform_real_distribution<> position_unif_dist(0, system_data.size_enclosure);
    normal_distribution<> mass_norm_dist{1E21, 1E15};

    double *force = (double *) malloc(sizeof(double) * system_data.num_objects * 3);
    double accel[3] = {0,0,0};

    /* Initialize x, y, z and m attributes of each object */
    for(int i = 0; i < system_data.num_objects; i++){
        objects.x[i] = position_unif_dist(gen64);
        objects.y[i] = position_unif_dist(gen64);
        objects.z[i] = position_unif_dist(gen64);
        objects.m[i] = mass_norm_dist(gen64);
        objects.active[i] = true;
    }
    /* Write initial configuration to a file*/
    write_config(0, system_data, objects);

    /* Initial collision checking */
    for(int i = 0; i < system_data.num_objects; i++){
        if( !objects.active[i] ){ continue; }
        for(int j = 0; j < system_data.num_objects; j++){
            if( i != j && objects.active[j])
                check_collision(objects, i, j);
        }
    }

    /* Body of the simulation */
    for(int i = 0; i < system_data.num_iterations; i++){
        for(int foo=0; foo < system_data.num_objects * 3; foo++){force[foo] = 0;}
        gravitational_force(system_data.num_objects, objects, system_data.time_step, force, accel);

        for(int a = 0; a < system_data.num_objects; a++){

            if (objects.active[a]){
                check_bounce(objects, a, system_data.size_enclosure);
            } else {
                continue;
            }
            for(int b = 0; b < system_data.num_objects; b++){
                if ( a != b && objects.active[b]){
                    check_collision(objects, a, b);
                }
            }
        }
    }

    /* Write final configuration to a file*/
    write_config(1, system_data, objects);
    free(force);
    free(objects.x);
    free(objects.y);
    free(objects.z);
    free(objects.vx);
    free(objects.vy);
    free(objects.vz);
    free(objects.m);
    free(objects.active);
    return 0;
}
