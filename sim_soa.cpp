//
// Created by mariwogr on 10/10/21.
//

#include <string>
#include <iostream>
#include <fstream>
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

double peso_generator(int random_seed) {

    // Motor
    std::mt19937_64 gen64(random_seed);
    // Distribución normal
    std::normal_distribution<> d{pow(10, 21), pow(10, 15)};

    // Peso
    double peso = d(gen64);

    return peso;
}

int parser(int argc, char* argv[]){

    if (argc != 5 ){
        return -1;
    }
    //comprobamos que numero de objetos sea mayor que cero o cero
    if ( ( (int) *argv[1] ) < 0){
        return -2;
    }


}

int gravitational_force_calc(struct parameters i, struct parameters j) {



}

int gravitational_force(int i, int j) {

    int grav_force;

    //for 

}

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
    int retcode = parser(argv);

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

    write_config(0, system_data);


    return 0;
}
