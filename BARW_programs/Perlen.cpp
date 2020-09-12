/***********************
 * Code Author: Adam Bryant
 * 
 * 06122020
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
#include "Perlen.h"

using namespace std;


int main(){
    /* Initalization */
    // Quantities of Interest
    int Total_attempts = 0;

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
    }

    // Thread ID's
    int ids[NUMBER_OF_THREADS];
    thread_obj *obj;
    void* retval;
    void* ptr;



    //Bias "s" initialization
    volatile double bias = 0.0;
    redprob = bias + 0.5;
    double sstep = 0.04;
    int number_of_ssteps = 1;

    //For loop parameters
    double INIT_RATE = 0.0;
    rate = INIT_RATE;
    int number_of_ysteps = 1;
    double rate_step = 0.003;

    // Growth parameters
    lambda = 0.005;
    beta = 2;
    radius = (LATTICE_SIZE*1.0)/(2*M_PI);

    // Constructing a variable filename with literals
    stringstream ss;
    ss << LATTICE_SIZE;
    string LATTICE_SIZEstr = ss.str(); ss.str("");

    // File Opening and initialization
    ofstream clout;
    string filename = "Perlen_N"+LATTICE_SIZEstr+".dat";
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
        obj->Dirs = new vector<vector<double>>;
        obj->Perlens = new vector<double>;
    }

    for(int s=0;s<number_of_ssteps;s++){
        for(int y=0;y<number_of_ysteps;y++){
            printf("b = %f\ts = %f\n",rate,bias); fflush(stdout);
            
            d2ddist = std::uniform_real_distribution<double> (0.0,(1 - cos(2 * (M_PI/double(LATTICE_SIZE))))/2);
            for(i=0;i<NUMBER_OF_THREADS;i++){

                // Spin of NUMBER of threads
                obj = thread_obj_array[i];
                thread_array[i] = std::thread(ThreadSpinner,obj);
            }

            for(i = 0; i < NUMBER_OF_THREADS; i++){
                thread_array[i].join();
            }

            for(i = 0;i < NUMBER_OF_THREADS; i++){
                obj = thread_obj_array[i];
                ClearObject(obj);
            }

            for(i = 0; i < MAX_GENERATIONS/3; i++){
                clout << obj->Perlens->at(i) << '\n';
            }
        }
    }        
    printf("Done...");


    for(i = 0; i < NUMBER_OF_THREADS; i++){
        obj = thread_obj_array[i];
        delete Trees[i];
        delete Rotors[i];
        delete Sets[i];
        delete obj;
    }

    clout.close();
    printf("Done.\n"); fflush(stdout);
    return 0;
}
