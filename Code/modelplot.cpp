/***********************
 * Code Author: Adam Bryant
 * 
 * 06122020
 * 
 * Driver for determining survival statistics on a BARW system, reworked for
 * parallel processes. This one runs only one full instance of a structure and
 * outputs a .vtk 4 scalar file for position and relative color of strains within the branching structure
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
#include "modelplot.h"

using namespace std;


int main(){
    /* Initalization */
    // Quantities of Interest       --      CHANGE THESE FOR DIFFERENT STRUCTURES
    double bias;
    bias = 0.01;
    rate = 0.03;
    LATTICE_SIZE = 125;
    MAX_GENERATIONS = 1000;
    Selective_color = true;



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
    redprob = bias + 0.5;

    //For loop parameters

    // Growth parameters
    lambda = 0.005;
    beta = 2;
    radius = (LATTICE_SIZE*1.0)/(2*M_PI);


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

    d2ddist = std::uniform_real_distribution<double> (0.0,(1 - cos(2*M_PI/double(LATTICE_SIZE)))/2);

    printf("b = %f\ts = %f\n",rate,bias); fflush(stdout);
    for(i=0;i<NUMBER_OF_THREADS;i++){

        // Spin off NUMBER of threads
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

    for(i = 0; i < NUMBER_OF_THREADS; i++){
        obj = thread_obj_array[i];
        delete Trees[i];
        delete Rotors[i];
        delete Sets[i];
        delete obj;
    }
    return 0;
}
