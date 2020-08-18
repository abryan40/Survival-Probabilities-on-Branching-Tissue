/***********************
 * Code Author: Adam Bryant
 * 
 * 06122020
 * 
 * Driver for determining survival statistics on a BARW system, reworked for
 * parallel processes. The "nd" indicates this version throws out the statistics
 * for processes that resulted in total branch death.
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
#include "BARW_nd_pp.h"

using namespace std;


int main(){
    /* Initalization */
    // Quantities of Interest
    double Average_time = 0;
    double Red_fixation_prob = 0;
    double Average_red_percent = 0;
    double Fixation_count = 0;
    double Fixation_time = 0;
    double Extinction_count = 0;
    double Extinction_time = 0;
    double Steady_state_count = 0;
    double Steady_state_time = 0;
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
        SetUpWalkers(*tree);
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
    int number_of_ssteps = 3;

    //For loop parameters
    double INIT_RATE = 0.057;
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
    ss << bias;
    string biasstr = ss.str();ss.str("");
    ss << rate;
    string ratestr = ss.str();ss.str("");
    ss << rate + (rate_step * (number_of_ysteps - 1));
    string rate_stepstr = ss.str();ss.str("");
    ss << sstep * (number_of_ssteps - 1);
    string sstepstr = ss.str();ss.str("");
    ss << MAX_GENERATIONS;
    string gen_str = ss.str();ss.str("");
    ss << NUMBER_OF_STRUCTURES;
    string struct_str = ss.str();ss.str("");
    //Debug initial lattice

    // File Opening and initialization
    ofstream clout;
    string type = "";	
    if(BARW){ type = "barw"; }
    else if(BRANCH){ type = "branch"; }
    else{ type = "stat"; }
    string filename = type+"_sb_nd_N"+LATTICE_SIZEstr+"_si"+biasstr+"_sf"+sstepstr+
        "_bi"+ratestr+"_bf"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    clout.open(filename.c_str());


    if(clout.fail()){
        cerr << "File failed to open.";
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

    for(int s=0;s<number_of_ssteps;s++){
        for(int y=0;y<number_of_ysteps;y++){
            printf("b = %f\ts = %f\n",rate,bias); fflush(stdout);
            for(i=0;i<NUMBER_OF_THREADS;i++){

                // Spin off threads
                obj = thread_obj_array[i];
                thread_array[i] = std::thread(ThreadSpinner,obj);
            }

            for(i = 0; i < NUMBER_OF_THREADS; i++){
                thread_array[i].join();
            }

            Red_fixation_prob = 0;
            Average_time = 0;
            Total_attempts = 0;

            // Aggregate statistics from the threads
            for(i = 0;i < NUMBER_OF_THREADS; i++){
                obj = thread_obj_array[i];
                Red_fixation_prob += obj->RED_SURV_PROB;
                Average_red_percent += obj->AVERAGE_RED_PERCENT;
                Fixation_count += obj->Fixation_count;
                Fixation_time += obj->Fixation_time;
                Extinction_count += obj->Extinction_count;
                Extinction_time += obj->Extinction_time;
                Steady_state_count += obj->Steady_state_count;
                Steady_state_time += obj->Steady_state_time;
                Total_attempts += obj->attempts;
                ClearObject(obj);
            }
            Red_fixation_prob   /= Total_attempts;
            Average_red_percent /= NUMBER_OF_THREADS;

            // Times are averaged over times they occured, or set to MAX+1 gens.
            if(Fixation_count)
                Fixation_time   /= Fixation_count;
            else
                Fixation_time = MAX_GENERATIONS+1;
            if(Extinction_count)
                Extinction_time /= Extinction_count;
            else
                Extinction_time = MAX_GENERATIONS+1;
            if(Steady_state_count)
                Steady_state_time/= Steady_state_count;
            else
                Steady_state_time = MAX_GENERATIONS+1;

            // The count variables are actually probabilities now
            Fixation_count      /= Total_attempts;
            Extinction_count    /= Total_attempts;
            Steady_state_count  /= Total_attempts;

            // Write to file
            clout << bias << " " << rate << " " <<
                Red_fixation_prob << " " << Average_red_percent << " " <<
                Fixation_count << " " << Fixation_time << " " <<
                Extinction_count << " " << Extinction_time << " " <<
                Steady_state_count << " " << Steady_state_time << " " <<
                Total_attempts <<  endl;

            // Reinitialize all the quantities for next loop
            Average_time = 0;
            Red_fixation_prob = 0;
            Average_red_percent = 0;
            Fixation_count = 0;
            Fixation_time = 0;
            Extinction_count = 0;
            Extinction_time = 0;
            Steady_state_count = 0;
            Steady_state_time = 0;
            Total_attempts = 0;


            rate += rate_step;
        }
        bias += sstep;
        redprob = 0.5 + bias;
        rate = INIT_RATE;
    }


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
