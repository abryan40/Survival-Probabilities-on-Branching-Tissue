/***********************
 * Code Author: Adam Bryant
 * 
 * 07072020
 * 
 * Driver for determining survival statistics on a BARW system, reworked for
 * parallel processes. Collects statistics relating to the state of the structure
 * at every generation timestep and averages for plots of characteristics such
 * as survival over time.
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
#include <algorithm>
#include <thread>
#include "BARW_t_pp.h"

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
    double Total_death_count = 0;
    double Total_death_time = 0;
    float SurProbs[MAX_GENERATIONS]; fill(SurProbs,SurProbs+MAX_GENERATIONS,0);
    float FixProbs[MAX_GENERATIONS]; fill(FixProbs,FixProbs+MAX_GENERATIONS,0);
    float ExtProbs[MAX_GENERATIONS]; fill(ExtProbs,ExtProbs+MAX_GENERATIONS,0);
    float StsProbs[MAX_GENERATIONS]; fill(StsProbs,StsProbs+MAX_GENERATIONS,0);
    float DieProbs[MAX_GENERATIONS]; fill(DieProbs,DieProbs+MAX_GENERATIONS,0);

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
    double sstep = 0.002;
    int number_of_ssteps = 1;

    //For loop parameters
    double INIT_RATE = 0;
    rate = INIT_RATE;
    int number_of_ysteps = 7;
    double rate_step = 0.01;

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
    ss << number_of_ysteps;
    string number_of_ystepsstr = ss.str();ss.str("");
    ss << rate_step;
    string rate_stepstr = ss.str();ss.str("");
    ss << MAX_GENERATIONS;
    string gen_str = ss.str();ss.str("");
    ss << NUMBER_OF_STRUCTURES;
    string struct_str = ss.str();ss.str("");
    //Debug initial lattice

    // File Opening and initialization
    ofstream Pout,Fout,Eout,Sout,Dout;
    string type = "";
    if(BARW){ type = "barw"; }
    else if(BRANCH){ type = "branch"; }
    else{ type = "stat"; }
    string filename = "Psurv_"+type+""+"_bt_N"+LATTICE_SIZEstr+"_s"+biasstr+"_b"
        +ratestr+"_ty"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    Pout.open(filename.c_str());
    filename = "Pfix_"+type+""+"_bt_N"+LATTICE_SIZEstr+"_s"+biasstr+"_y"
        +ratestr+"_ty"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    Fout.open(filename.c_str());
    filename = "Pext_"+type+""+"_bt_N"+LATTICE_SIZEstr+"_s"+biasstr+"_y"
        +ratestr+"_ty"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    Eout.open(filename.c_str());
    filename = "Psts_"+type+""+"_bt_N"+LATTICE_SIZEstr+"_s"+biasstr+"_y"
        +ratestr+"_ty"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    Sout.open(filename.c_str());
    filename = "Pdie_"+type+""+"_bt_N"+LATTICE_SIZEstr+"_s"+biasstr+"_y"
        +ratestr+"_ty"+rate_stepstr+"_gens"+gen_str+"_str"+struct_str+".dat";
    Dout.open(filename.c_str());


    if(Pout.fail()){
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

    d2ddist = std::uniform_real_distribution<double>(0.0,(1-cos(2*M_PI/double(LATTICE_SIZE)))/2);

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

            // Aggregate data from arrays
            for(i = 0;i < NUMBER_OF_THREADS; i++){
                obj = thread_obj_array[i];
                for(int j=0;j<MAX_GENERATIONS;j++){
                    SurProbs[j] += obj->SurProbs[j];
                    FixProbs[j] += obj->FixProbs[j];
                    ExtProbs[j] += obj->ExtProbs[j];
                    StsProbs[j] += obj->StsProbs[j];
                    DieProbs[j] += obj->DieProbs[j];
                }
                ClearObject(obj);
            }

            // Average
            for(i = 0; i < MAX_GENERATIONS; i++){
                SurProbs[i] /= NUMBER_OF_STRUCTURES;
                FixProbs[i] /= NUMBER_OF_STRUCTURES;
                ExtProbs[i] /= NUMBER_OF_STRUCTURES;
                StsProbs[i] /= NUMBER_OF_STRUCTURES;
                DieProbs[i] /= NUMBER_OF_STRUCTURES;
            }
            
            // Write to file
            Pout << rate;
            for(i = 0; i < MAX_GENERATIONS; i++){
                Pout << " " << SurProbs[i]; SurProbs[i] = 0;}
            Pout << endl;
            Fout << rate;
            for(i = 0; i < MAX_GENERATIONS; i++){
                Fout << " " << FixProbs[i]; FixProbs[i] = 0;}
            Fout << endl;
            Eout << rate;
            for(i = 0; i < MAX_GENERATIONS; i++){
                Eout << " " << ExtProbs[i]; ExtProbs[i] = 0;}
            Eout << endl;
            Sout << rate;
            for(i = 0; i < MAX_GENERATIONS; i++){
                Sout << " " << StsProbs[i]; StsProbs[i] = 0;}
            Sout << endl;
            Dout << rate;
            for(i = 0; i < MAX_GENERATIONS; i++){
                Dout << " " << DieProbs[i]; DieProbs[i] = 0;}
            Dout << endl;

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
            Total_death_count = 0;
            Total_death_time = 0;


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

    Pout.close();
    Fout.close();
    Eout.close();
    Sout.close();
    Dout.close();
    printf("Done.\n"); fflush(stdout);
    return 0;
}
