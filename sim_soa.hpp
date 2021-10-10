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
    double x[];
    double y[];
    double z[];
    double vx[];
    double vy[];
    double vz[];
    double ax[];
    double ay[];
    double az[];
    double m[];
};

int write_config(int id, parameters params);
int print_error_args(int argc, char* argv[]);

#endif //GRAVITY_SIMULATOR_SIM_SOA_HPP
