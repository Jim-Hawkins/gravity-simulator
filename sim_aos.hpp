//
// Created by mariwogr on 16/10/21.
//

#ifndef GRAVITY_SIMULATOR_SIM_AOS_HPP
#define GRAVITY_SIMULATOR_SIM_AOS_HPP


struct parameters {
    int num_objects;
    int num_iterations;
    int random_seed;
    double size_enclosure;
    double time_step;
};

struct set {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double m;
    bool active;
};

double accel_calc(double m, double F);
int check_bounce(set objects, int obj, double size);
int check_collision(set objects, int i, int j);
int collision_objects(set objects,int i, int j);
int gravitational_force(int num_objects, set objects, double time_step);
void gravitational_force_calc(set objects, int i, int j, double *force);
int parser(int argc, char* argv[]);
int print_error_args(int argc, char* argv[]);
int write_config(int id, parameters system_data, set params);

#endif //GRAVITY_SIMULATOR_SIM_AOS_HPP
