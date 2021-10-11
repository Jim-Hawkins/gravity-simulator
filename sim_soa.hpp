//
// Created by mariwogr on 10/10/21.
//
#define MAX_OBJECTS 50000

#ifndef GRAVITY_SIMULATOR_SIM_SOA_HPP
#define GRAVITY_SIMULATOR_SIM_SOA_HPP

struct parameters {
    int num_objects;
    int num_iterations;
    int random_seed;
    float size_enclosure;
    float time_step;
};

struct set {
    double x[MAX_OBJECTS];
    double y[MAX_OBJECTS];
    double z[MAX_OBJECTS];
    double vx[MAX_OBJECTS];
    double vy[MAX_OBJECTS];
    double vz[MAX_OBJECTS];
    double m[MAX_OBJECTS];
    bool active[MAX_OBJECTS]{false};
};

double accel_calc(double m, double F);
int gravitational_force(int num_objects, set objects, float time_step);
double * gravitational_force_calc(set objects, int i, int j);
int parser(int argc, char* argv[]);
int print_error_args(int argc, char* argv[]);
double weight_generator(int random_seed);
int write_config(int id, parameters system_data, set params);

#endif //GRAVITY_SIMULATOR_SIM_SOA_HPP
