/***********************
 * Code Author: Adam Bryant
 * 
 * 07072020
 * 
 * Driver for determining survival statistics on a BARW system, optimized for
 * parallel programming.
 **********************/
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <random>
#include <thread>
#include "BARW_Nb_pp.h"

using namespace std;


int main(){
    /* Initalization */
    // Quantities of Interest
    double Active_branches = 0;
    double Annihilations = 0;
    double Total_branches = 0;
    double Non_parent_branches = 0;
    double Total_death = 0;

    /* Allocating memory for thread objects */
    vector<VB*> Trees;
    vector<rotor*> Rotors;
    vector<set<branch*>*> Sets;
    rotor* rotptr;
    branch* latt;
    set<branch*>* setptr;
    for(int i = 0; i < NUMBER_OF_THREADS; i++){
        gen genr;
        genr.pos_arr[0] = 0.0;
        genr.pos_arr[1] = 0.0;
        genr.pos_arr[2] = 0.0;
        genr.color = RedCellNumber/LATTICE_SIZE;

        rotptr = new rotor;
        setptr = new set<branch*>;
        Sets.push_back(setptr);
        Rotors.push_back(rotptr);

        latt = new branch;
        setptr->insert(latt);
        latt->lattice_size = LATTICE_SIZE;
        latt->lattice_size_actual = LATTICE_SIZE*1.0;
        latt->dir_arr[0] = 0;
        latt->dir_arr[1] = 0;
        latt->dir_arr[2] = 1;
        latt->surface.push_back(genr);
        VB* tree = new VB;
        (*tree).push_back(latt);
        Trees.push_back(tree);
        SetUpWalkers(*tree);
    }

    // Thread ID's
    int ids[NUMBER_OF_THREADS];
    thread_obj *obj;
    void* retval;
    void* ptr;



    //Bias "s" initialization
    volatile double bias = 0.000;
    redprob = bias + 0.5;
    int Nstep = 5;
    int number_of_Nsteps = 21;

    //For loop parameters
    double INIT_RATE = 0.00;
    rate = INIT_RATE;
    int number_of_bsteps = 21;
    double rate_step = 0.003;

    // Growth parameters
    lambda = 0.005;
    beta = 2;
    radius = (LATTICE_SIZE*1.0)/(2*M_PI);

    // Constructing a variable filename with literals
    stringstream ss;
    ss << LATTICE_SIZE;
    string LATTICE_SIZEstr = ss.str(); ss.str("");
    ss << rate;
    string ratestr = ss.str();ss.str("");
    ss << number_of_bsteps;
    string number_of_bstepsstr = ss.str();ss.str("");
    ss << INIT_RATE + (rate_step * (number_of_bsteps - 1));
    string rate_stepstr = ss.str();ss.str("");
    ss << LATTICE_SIZE + ((number_of_Nsteps-1) * Nstep);
    string Nstepsstr = ss.str();ss.str("");
    ss << MAX_GENERATIONS;
    string gen_str = ss.str();ss.str("");
    ss << NUMBER_OF_STRUCTURES;
    string struct_str = ss.str();ss.str("");
    //Debug initial lattice

    // File Opening and initialization
    ofstream clout;
    string filename = "BARW_Nb_Ni"+LATTICE_SIZEstr+"_Nf"+Nstepsstr+
        "_bi"+ratestr+"_bf"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    clout.open(filename.c_str());


    if(clout.fail()){
        cout << "File failed to open.";
        return 1;
    }

    int i;
    thread_obj* thread_obj_array[NUMBER_OF_THREADS];
    std::thread thread_array[NUMBER_OF_THREADS];
    for(i=0;i<NUMBER_OF_THREADS;i++){
        thread_obj_array[i] = new thread_obj;
        obj = thread_obj_array[i];
        obj->tree = Trees[i];
        obj->active = Sets[i];
        obj->rot  = *Rotors[i];
        obj->ID = i;
    }

    for(int n=0;n<number_of_Nsteps;n++){
        for(int y=0;y<number_of_bsteps;y++){
            printf("N = %d\tb = %f\n",LATTICE_SIZE,rate); fflush(stdout);
            
            // Change angle of rotation distribution
            d2ddist = std::uniform_real_distribution<double> (0.0,(1 - cos(2*M_PI/double(LATTICE_SIZE)))/2);
            radius = (1.0*LATTICE_SIZE)/(2*M_PI);

            for(i=0;i<NUMBER_OF_THREADS;i++){

                // Spin off threads
                obj = thread_obj_array[i];
                thread_array[i] = std::thread(ThreadSpinner,obj);
            }

            for(i = 0; i < NUMBER_OF_THREADS; i++){
                thread_array[i].join();
            }

            // Aggregate data for averaging
            for(i = 0;i < NUMBER_OF_THREADS; i++){
                obj = thread_obj_array[i];
                Active_branches += obj->Active_branches;
                Annihilations += obj->Annihilations;
                Total_branches += obj->Total_branches;
                Non_parent_branches += obj->Non_parent_branches;
                Total_death += obj->Total_death;
                ClearObject(obj);
            }

            // The count variables are actually probabilities now
            Active_branches /= NUMBER_OF_STRUCTURES;
            Annihilations /= NUMBER_OF_STRUCTURES;
            Total_branches /= NUMBER_OF_STRUCTURES;
            Non_parent_branches /= NUMBER_OF_STRUCTURES;
            Total_death /= NUMBER_OF_STRUCTURES;
            // Write to file
            clout << LATTICE_SIZE << " " << rate << " " <<
                Active_branches << " " << Annihilations << " " <<
                Total_branches << " " << Non_parent_branches << " " <<
                Total_death <<  endl;
            
            // Reinitialize all the quantities for next loop
            Active_branches = 0;
            Annihilations = 0;
            Total_branches = 0;
            Non_parent_branches = 0;
            Total_death = 0;
            
            rate += rate_step;
        }
        LATTICE_SIZE += Nstep;
        rate = INIT_RATE;
    }

    for(i = 0; i < NUMBER_OF_THREADS; i++){
        obj = thread_obj_array[i];
        delete Trees[i];
        delete Rotors[i];
        delete Sets[i];
        delete obj;
    }

    printf("Done.\n");
    clout.close();
    return 0;
}
