//
// Created by mariwogr on 10/10/21.
//

#ifndef GRAVITY_SIMULATOR_SIM_SOA_HPP
#define GRAVITY_SIMULATOR_SIM_SOA_HPP

struct parameters {
    int num_objects;
    int num_iterations;
    int random_seed;
    float size_enclosure;
    float time_step
};

struct set {
    double x[MAX_OBJECTS];
    double y[MAX_OBJECTS];
    double z[MAX_OBJECTS];
    double vx[MAX_OBJECTS];
    double vy[MAX_OBJECTS];
    double vz[MAX_OBJECTS];
    double m[MAX_OBJECTS];
};

int write_config(int id, parameters params);
int print_error_args(int argc, char* argv[]);

#endif //GRAVITY_SIMULATOR_SIM_SOA_HPP
