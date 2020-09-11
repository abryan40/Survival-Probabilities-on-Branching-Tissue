/* Code author: Adam Bryant
 * summer 2020
 * 
 * The header that holds every important function for the BARW cell lattice model.
 * Modified from BARW_nd.h for collection of probabilities of various states as a function
 * of time.
 */


#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstdio>
#include <random>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <set>
#include <pthread.h>
#include <thread>
#include <algorithm>
// Self-defined 
#include "matrix.h"

using namespace std;

/* Globals */
int LATTICE_SIZE = 50;      // N_0
int RedCellNumber = 1;      // x_0
double lambda = 0.005;
double beta = 2;
volatile double redprob;
volatile double rate;       // b
double length_unit = 1;
short annihilations;
double radius;
int NUMBER_OF_THREADS = 25;
int NUMBER_OF_STRUCTURES = 2000;
int WORKLOAD = NUMBER_OF_STRUCTURES/NUMBER_OF_THREADS;
int MAX_GENERATIONS = 1000;
/* Control variables */
bool BARW = true;       // If false, does not update position (no termination)
bool BRANCH = true;

typedef vector< vector<double > > M;

// Random number initialization using mt19937
std::random_device rd;
std::mt19937 mt(time(NULL));
std::uniform_real_distribution<double> dist(0.0,1.0);
std::uniform_real_distribution<double> realdist(-1.0,1.0);
std::uniform_real_distribution<double>  branchdist(0.0,1.0);
std::uniform_real_distribution<double> alphadist    ((1-cos((100/3)*(M_PI/180)))/2,
        ((1-cos((200/3)*(M_PI/180))))/2);
std::uniform_real_distribution<double> pidist(0.0,2*M_PI);
std::uniform_real_distribution<double> d2ddist(0.0,(1 - cos(2 * (M_PI/double(100))))/2);

// Populates the "backbone" vector<gen> surface of the branch structure.
// The pos_arr[3] is the center 
struct gen{
    float pos_arr[3];
    float color;        // Number 0.0 - 1.0, red percentage
};

struct branch{
    vector<gen> surface;
    short walker1;
    short walker2;
    short prevwalker1;
    short prevwalker2;
    unsigned short growgen;
    unsigned short branchGeneration;
    bool termination;
    bool fixation;
    bool growing;
    int generation;
    char fixcolor;
    int lattice_size;
    double lattice_size_actual;
    double radius;
    double pos_arr[3] = {0.,0.,0.};
    double dir_arr[3] = {0.,0.,1.};
    double pos_wall[6] = {0.,0.,0.,0.,0.,0.}; //mx,Mx,my,My,mz,Mz
    branch* parent;
    branch* sister;
    branch(){
        fixation = false;
        fixcolor = '0';
        lattice_size = LATTICE_SIZE;
        lattice_size_actual = LATTICE_SIZE*1.0;
        radius = lattice_size_actual/(2*M_PI);
        termination = false;
        growing = false;
        growgen = 1;
        generation = 1;
        branchGeneration = 1;
        parent = NULL;
        sister = NULL;
    };
};

typedef vector<branch*> VB;

struct rotor{
    double s,c;
    vector<double> a_vec,z_vec,new_dir,a2_vec,v_vec;
    M sksym,ID3,R,v_x2,Inv_R;
    rotor(){
        a_vec.resize(3);
        z_vec.resize(3);
        z_vec[0] = 0;
        z_vec[1] = 0;
        z_vec[2] = 1;
        new_dir.resize(3);
        a2_vec.resize(3);
        v_vec.resize(3);
        sksym.resize(3);
        sksym[0].resize(3);
        sksym[1].resize(3);
        sksym[2].resize(3);
        ID3 = I(3);
        v_x2.resize(3);
        v_x2[0].resize(3);
        v_x2[1].resize(3);
        v_x2[2].resize(3);
        R.resize(3);
        R[0].resize(3);
        R[1].resize(3);
        R[2].resize(3);
        Inv_R.resize(3);
        Inv_R[0].resize(3);
        Inv_R[1].resize(3);
        Inv_R[2].resize(3);
    }
};

// A struct passed to threads that holds:
/*
 * rotor: Collection of presized matrices for rotation
 * tree : The branching structure
 * AVERAGE_TIME: Variable for average fixation time for red cells
 * RED_SURV_PROB: Variable for P_fix
 * ID: Integer ID
 */
struct thread_obj{
    rotor rot;
    VB* tree;
    set<branch*>* active;
    set<branch*>::iterator sit;
    short SurProbs[1000];
    bool Fixation;    
    int Fixation_count = 0;
    int Fixation_time = 0;
    short FixProbs[1000];
    bool Extinction = false;
    int Extinction_count = 0;
    int Extinction_time = 0;
    short ExtProbs[1000];
    //int Extinction_time[WORKLOAD];
    bool Steady_State = false;
    int Steady_state_count = 0;
    int Steady_state_time = 0;
    short StsProbs[1000];
    //int Steady_state_time[WORKLOAD];
    bool Total_death = false;
    int Total_death_count = 0;
    int Total_death_time = 0;
    short DieProbs[1000];
    //int Total_death_times[WORKLOAD];
    int NUMBER_OF_GENERATIONS;
    double AVERAGE_TIME = 0;
    double RED_SURV_PROB = 0;
    double AVERAGE_RED_PERCENT;
    int ID;
};

// Lattice competition functions
bool CellEvolution(thread_obj*,int);
void Competition(thread_obj*);
void Competition2(thread_obj*);
void CheckForFixation(VB&);
void SetUpWalkers(VB&);
void Status(const VB&);

// CHANGE TO POINTERS :/
template <typename T>
double avg(vector<T>);
template <typename T>
double avg(T*,int);

// Termination functions
bool AreAnyRed(VB&);
bool AreTheyAllFixated(thread_obj*);
bool AreHalfFixated(VB&);
bool AreTheyAllTerminated(thread_obj*);
double RedPercent(thread_obj*);
bool minicomp();

// Functions that change population numbers
void growth(branch&,VB&);
void BranchProcess(thread_obj*,branch*,int);
void Reinitialize(thread_obj*);

// Positional functions
double PowerLaw(int);
double Length(int);
void InitializePos_Dir(branch&,branch&,branch*,rotor&);
void UpdatePos(thread_obj*);
void RotateDir(thread_obj*);
void DirMut(branch&,rotor&,double,double,double);
void RecordPosData(const VB&,ofstream&,int);
void CollisionDetection(VB&,set<branch*>*);

/*******************************
 * 
 *  ~ Function Definitions ~
 * 
 ******************************/
void* ThreadSpinner(void* t){
    thread_obj* obj = (thread_obj*) t;
    for(int i = 0; i < WORKLOAD; i++){
        if(!CellEvolution(obj,i)){
            i--;
        }
        
        //Status(obj->tree);
        Reinitialize(obj);
    }
    return (void*) obj;
}

void ClearObject(thread_obj* o){
    o->Fixation = false;
    o->Fixation_count = 0;
    o->Fixation_time = 0;
    o->Extinction = false;
    o->Extinction_count = 0;
    o->Extinction_time = 0;
    o->Steady_State = false;
    o->Steady_state_count = 0;
    o->Steady_state_time = 0;
    o->Total_death = false;
    o->Total_death_count = 0;
    o->Total_death_time = 0;
    o->NUMBER_OF_GENERATIONS = 0;
    o->AVERAGE_TIME = 0;
    o->RED_SURV_PROB = 0;
    o->AVERAGE_RED_PERCENT = 0;
    fill(o->SurProbs,o->SurProbs+MAX_GENERATIONS,0);
    fill(o->FixProbs,o->FixProbs+MAX_GENERATIONS,0);
    fill(o->ExtProbs,o->ExtProbs+MAX_GENERATIONS,0);
    fill(o->StsProbs,o->StsProbs+MAX_GENERATIONS,0);
    fill(o->DieProbs,o->DieProbs+MAX_GENERATIONS,0);
}
//The main program within the main program
bool CellEvolution(thread_obj* obj,int index){
    int i;
    branch* b;
    int ss_time = 0;
    while(obj->NUMBER_OF_GENERATIONS < MAX_GENERATIONS){
        obj->NUMBER_OF_GENERATIONS++;

        // All branches are fixated, check for homogeneity for statistics.
        if(!obj->Steady_State && AreTheyAllFixated(obj)){
            obj->Steady_State = true;
            ss_time = obj->NUMBER_OF_GENERATIONS;
            obj->Steady_state_time += ss_time;
            obj->Steady_state_count++;

            char fixcolor;
            for(i = 0; i < obj->tree->size(); i++){
                if(obj->tree->at(i) != NULL && !obj->tree->at(i)->termination){
                    fixcolor = obj->tree->at(i)->fixcolor; break;
                }
            }
            bool homogeneous = true;
            for(i = 0; i < obj->tree->size(); i++){
                if( obj->tree->at(i) != NULL && !obj->tree->at(i)->termination &&
                     obj->tree->at(i)->fixcolor != fixcolor){
                    homogeneous = false; break;
                }
            }
            if(homogeneous){
                if(fixcolor == 'r'){
                    obj->Fixation = true;
                    obj->Fixation_time += ss_time;
                    obj->Fixation_count++;
                    //printf("Fixating after %d gens.\n",obj->NUMBER_OF_GENERATIONS);
                    if(!BARW){ return true; }
                }else if(fixcolor == 'g'){
                    obj->Extinction = true;
                    obj->Extinction_time += ss_time;
                    obj->Extinction_count++;
                    if(!BARW){ return false; }
                }
            }else{
                // Returns if steady state is reached not entire green. Does away with total death.
                if(!BARW){ return true; }
            }

            // Total death is a wholly different entity, and is the only thing in
            // BARW that can end the structure before generation 1000.
        }else if(BARW && AreTheyAllTerminated(obj)){
            obj->Total_death = true;
            if(obj->Extinction){
                obj->Extinction_time -= ss_time; obj->Extinction_count--;
            }
            if(obj->Fixation){
                obj->Fixation_time -= ss_time; obj->Fixation_count--;
            }
            if(obj->Steady_State){
                obj->Steady_state_time -= ss_time; obj->Steady_state_count--;
            }

            // If not counting total death, this rectifies
            if(ss_time){
                if(obj->Fixation){
                    for(int k=ss_time-1;k < obj->NUMBER_OF_GENERATIONS-1;k++){
                        obj->StsProbs[k]--;
                        obj->FixProbs[k]--;
                    }
                    for(int k=0;k < obj->NUMBER_OF_GENERATIONS-1;k++)
                        obj->SurProbs[k]--;
                }else if(obj->Extinction){
                    for(int k=ss_time-1;k < obj->NUMBER_OF_GENERATIONS-1;k++){
                        obj->ExtProbs[k]--;
                        obj->StsProbs[k]--;
                    }
                    for(int k=0;k < ss_time-1;k++)
                        obj->SurProbs[k]--;
                }else{
                    for(int k=ss_time-1;k < obj->NUMBER_OF_GENERATIONS-1;k++){
                        obj->StsProbs[k]--;
                    }
                    for(int k=0;k < obj->NUMBER_OF_GENERATIONS-1;k++)
                        obj->SurProbs[k]--;
                }
            }else{
                for(int k=0;k < obj->NUMBER_OF_GENERATIONS-1;k++)
                    obj->SurProbs[k]--;
            }
            obj->Total_death_time += obj->NUMBER_OF_GENERATIONS;
            obj->Total_death_count++;
            /*
            for(i = obj->NUMBER_OF_GENERATIONS-1; i < MAX_GENERATIONS; i++){
                obj->DieProbs[i] += 1;
            }
            */
            return false;
        }

        // Main computation functions

        // Cell competition
        if(!obj->Steady_State && obj->NUMBER_OF_GENERATIONS%2){
            Competition(obj);
        }else{
            Competition2(obj);
        }

        if(BRANCH && !(obj->active->empty())){
            for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
                b = *(obj->sit);
                if(!b->growing){
                    if(dist(mt) < (rate)){
                        b->growing = true;
                        growth(*b,*obj->tree);
                    }
                }else{
                    growth(*b,*obj->tree);
                }
            }
        }


        if(!obj->Steady_State){
            CheckForFixation(*(obj->tree));
        }


        UpdatePos(obj);


        if(BARW){
            CollisionDetection(*(obj->tree),obj->active);
            RotateDir(obj);
        }

        // Growth and bifurcation loops
        if(BRANCH && !(obj->active->empty())){
            if(BARW){
                for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
                    b = *(obj->sit);
                    if(b->growing && !(b->termination)){
                        if(b->lattice_size >= 2*LATTICE_SIZE){
                            b->growing = false;
                            BranchProcess(obj,*(obj->sit),0);

                            // This hurts, but to keeps things safe
                            obj->sit = obj->active->begin();
                        }
                    }
                }
            }else{
                for(i = 0; i < obj->tree->size(); i++){
                    if(obj->tree->at(i) != NULL && obj->tree->at(i)->growing &&
                            !obj->tree->at(i)->termination){
                        if(obj->tree->at(i)->lattice_size >= 2*LATTICE_SIZE){
                            obj->tree->at(i)->growing = false;
                            BranchProcess(obj,obj->tree->at(i),i);
                        }
                    }
                }
            }
        }

        // Filling in time dependant quantities
        if(obj->Steady_State){
            obj->StsProbs[obj->NUMBER_OF_GENERATIONS-1] += 1;
            if(obj->Extinction){
                obj->ExtProbs[obj->NUMBER_OF_GENERATIONS-1] += 1;
            }else{
                if(obj->Fixation){
                    obj->FixProbs[obj->NUMBER_OF_GENERATIONS-1] += 1;
                }
                obj->SurProbs[obj->NUMBER_OF_GENERATIONS-1] += 1;
            }
        }else{
            obj->SurProbs[obj->NUMBER_OF_GENERATIONS-1] += 1;
        }
    }
    // At generation 1000, mutant strain survives if not extinct
    //if(obj->Extinction) return false;
    return true;
}

// Stochastic walker growth, stagnation, or decay checks
// Growth also occurs here
void Competition(thread_obj* obj){
    branch* b;
    for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
        b = *(obj->sit);
        b->generation++;
        if(!b->fixation){
            b->prevwalker1 = b->walker1;
            b->prevwalker2 = b->walker2;

            //Last cell checks
            //b->Walker1 - competes for growth
            if(dist(mt)<=redprob){
                b->walker1 --;

                //Wrapping mechanic
                if(b->walker1<=-1){ //Smallest space is 0
                    b->walker1 = b->lattice_size - 1;
                }
            }

            //b->Walker2 - competes against loss
            if(dist(mt)>=redprob){
                b->walker2 --;

                //Wrapping mechanic
                if(b->walker2<=-1){
                    b->walker2 = b->lattice_size - 1;
                }
            }
        }
    }
}


// Competition2() inputs the result of 1 v. 2 into 2, 2 v. 3 into 3, etc.
// The value for 1 v. LAST is put into 1.
void Competition2(thread_obj* obj){
    branch* b;
    for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
        b = *(obj->sit);
        b->generation++;
        if(!b->fixation){

            b->prevwalker1 = b->walker1;
            b->prevwalker2 = b->walker2;

            //Last cell checks
            //b->Walker1 - competes against loss
            if(dist(mt)>=redprob){
                b->walker1 ++;

                //Wrapping mechanic
                if(b->walker1>=(b->lattice_size)){ //Largest space is LATTICE_SIZE-1
                    b->walker1 = 0;
                }
            }

            //b->Walker2 - competes for growth
            if(dist(mt)<=redprob){
                b->walker2 ++;

                //Wrapping mechanic
                if(b->walker2>=(b->lattice_size)){
                    b->walker2 = 0;
                }
            }
        }
    }
}

// This seems complicated, but is better than checking every single cell
// It is also size invariant, only ever keeping track of 4 ints
void CheckForFixation(VB& tree){
    for(int i=0;i<tree.size();i++){
        if(tree[i] != NULL && !tree[i]->fixation){
            if(tree[i]->walker1==tree[i]->walker2){
                if((tree[i]->prevwalker1==0)&&
                        (tree[i]->prevwalker2==tree[i]->lattice_size-1)){
                    tree[i]->fixcolor = 'r';
                    tree[i]->fixation = true;
                }else if((tree[i]->prevwalker1==tree[i]->lattice_size-1)&&
                        (tree[i]->prevwalker2==0)){
                    tree[i]->fixcolor = 'g';
                    tree[i]->fixation = true;
                }else if(tree[i]->prevwalker1>
                        tree[i]->prevwalker2){
                    tree[i]->fixcolor = 'r';
                    tree[i]->fixation = true;
                }else if(tree[i]->prevwalker2>
                        tree[i]->prevwalker1){
                    tree[i]->fixcolor = 'g';
                    tree[i]->fixation = true;
                }else{
                    cout << "\nUnknown fixation condition";
                    cout << "L = " << LATTICE_SIZE << endl;
                    Status(tree);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    return;
}

// Total system fixation. Often occurs very early.
bool AreTheyAllFixated(thread_obj* obj){
    branch* b;
    for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
        b = *(obj->sit);
        if(!b->fixation){
            return false;
        }
    }
    return true;
}

// Can check for majority fixation. Not used.
bool AreHalfFixated(VB& tree){
    int sum=0;
    for(int i=0;i<tree.size();i++){
        if(tree[i]->fixation){
            sum++;
        }
    }  
    if(tree.size()/2 < sum){
        return true;
    }else{
        return false;
    }
}

// Check for total termination, a self-annihilated system.
bool AreTheyAllTerminated(thread_obj* obj){
    obj->sit = obj->active->begin();
    branch* b;
    while(obj->sit != obj->active->end()){
        b = *(obj->sit);
        if(!(b->termination)){
            return false;
        }
        // If a branch is terminated, remove it from active
        obj->active->erase(obj->sit);
        if(obj->active->size() < 1){
            return true;
        }else{
            obj->sit = obj->active->begin();
        }
    }
    return true;
}

double RedPercent(thread_obj* obj){
    if(obj->Total_death) return 0;
    if(obj->Extinction) return 0;
    if(obj->Fixation) return 1;
    double sum = 0;
    for(int i=0;i<obj->tree->size();i++){
        if( obj->tree->at(i) != NULL && !obj->tree->at(i)->termination &&
            obj->tree->at(i)->surface.size() > 1){
            if(obj->tree->at(i)->fixcolor == 'r')
                sum += 1;
            else if(obj->tree->at(i)->fixcolor != 'g')
                sum += obj->tree->at(i)->surface[obj->tree->at(i)->surface.size() - 1].color;
        }
    }
    if(obj->active->size() > 0)
        return sum/obj->active->size();
    else
        return 0;
}

// Defines walker numbers based on globals
void SetUpWalkers(VB& tree){
    for(int i=0;i<tree.size();i++){
        tree[i]->walker1 = (LATTICE_SIZE/2) - ((RedCellNumber/2) + (RedCellNumber%2));
        tree[i]->walker2 = tree[i]->walker1 + RedCellNumber;
        tree[i]->prevwalker1 = tree[i]->walker1;
        tree[i]->prevwalker2 = tree[i]->walker2;
    }
}


// Growth according to Length() function, determines number of
// new cells to add.
void growth(branch& branch,VB& tree){

    //First, determine how many new cells to add
    branch.lattice_size_actual = Length(branch.growgen+1);
    int newcells = floor(branch.lattice_size_actual - branch.lattice_size);

    //Truncates spilloff at the end
    if(branch.lattice_size_actual > (2.0*LATTICE_SIZE)){
        newcells = 2*LATTICE_SIZE - branch.lattice_size;
    }
    //Now we have the number of new cells and the actual increase in size
    //This loop choses a number of random slots for new cells
    bool gw1,gw2;
    if(branch.walker1 > branch.walker2){
        gw1 = true;
        gw2 = false;
    }else if(branch.walker1 < branch.walker2){
        gw1 = false;
        gw2 = true;
    }else{
        gw1 = false; gw2 = false;
    }


    for(int j=0;j<newcells;j++){
        if(branch.fixation){
            branch.lattice_size++;
            //Non-fixed branches
        }else{
            int newspot = floor(dist(mt)*branch.lattice_size);

            //Not on either walker
            if((newspot!=branch.walker1)&&(newspot!=branch.walker2)){
                if(gw1){
                    //Red
                    if(newspot > branch.walker1){
                        branch.lattice_size++;
                        //Green
                    }else if((newspot < branch.walker1)&&
                            (newspot > branch.walker2)){
                        branch.walker1++;
                        branch.lattice_size++;
                        //Red
                    }else if((newspot < branch.walker1)&&
                            (newspot < branch.walker2)){
                        branch.walker1++;
                        branch.walker2++;
                        branch.lattice_size++;
                    }
                }else if(gw2){
                    //Green
                    if(newspot > branch.walker2){
                        branch.lattice_size++;
                        //Green
                    }else if((newspot < branch.walker2)&&
                            (newspot > branch.walker1)){
                        branch.walker2++;
                        branch.lattice_size++;
                        //Red
                    }else if((newspot < branch.walker1)&&
                            (newspot < branch.walker2)){
                        branch.walker1++;
                        branch.walker2++;
                        branch.lattice_size++;
                    }
                }
                //New spot is on walker 1
            }else if(newspot == branch.walker1){
                int newcell = minicomp();
                if(gw1){
                    if(newcell){
                        branch.lattice_size++;
                    }else{
                        branch.walker1++;
                        branch.lattice_size++;    
                    }
                }else if(gw2){
                    if(newcell){
                        branch.walker2++;
                        branch.lattice_size++;
                    }else{
                        branch.walker1++;
                        branch.walker2++;
                        branch.lattice_size++;
                    }
                }
            }else if(newspot == branch.walker2){
                int newcell = minicomp();
                if(gw1){
                    if(newcell){
                        branch.walker1++;
                        branch.walker2++;
                        branch.lattice_size++;
                    }else{
                        branch.walker1++;
                        branch.lattice_size++;
                    }
                }else if(gw2){
                    if(newcell){
                        branch.walker2++;
                        branch.lattice_size++;
                    }else{
                        branch.lattice_size++;
                    }
                }
            }
            newspot = 0;
        }
    }
    branch.growgen++;
    branch.radius = branch.lattice_size/(2*M_PI);
}

// Adds new branches in place of the previous tip
void BranchProcess(thread_obj* obj,branch* b,int i){
    int initsize = obj->tree->size();

    int rand = floor(branchdist(mt)*b->lattice_size);
    b->walker1 += rand;
    b->walker2 += rand;
    if(b->walker1 >= b->lattice_size)   b->walker1 -= b->lattice_size;   
    if(b->walker2 >= b->lattice_size)   b->walker2 -= b->lattice_size;

    branch* newcells1 = new branch;
    branch* newcells2 = new branch;
    obj->active->insert(newcells1);
    obj->active->insert(newcells2);
    newcells1->branchGeneration = b->branchGeneration + 1;
    newcells2->branchGeneration = b->branchGeneration + 1;
    newcells1->sister = newcells2;
    newcells2->sister = newcells1;
    newcells1->parent = b;
    newcells2->parent = b;

    InitializePos_Dir(*newcells1,*newcells2,b,obj->rot);
    obj->tree->push_back(newcells1);
    obj->tree->push_back(newcells2);


    if((b->walker1 > b->lattice_size)||
            (b->walker2 > b->lattice_size)){
        if(!b->fixation){
            cout << "\n\nSomehow, the walkers are exceeding lattice size...";
            Status(*(obj->tree));
            exit(EXIT_FAILURE);
        }
    }

    //Already fixated branch
    if(b->fixation){
        obj->tree->at(initsize)->walker1 = 0;
        obj->tree->at(initsize+1)->walker2 = 0;
        obj->tree->at(initsize)->walker2 = 0;
        obj->tree->at(initsize+1)->walker1 = 0;
        obj->tree->at(initsize+1)->fixation = true;
        obj->tree->at(initsize+1)->fixcolor = b->fixcolor;
        obj->tree->at(initsize)->fixation = true;
        obj->tree->at(initsize)->fixcolor = b->fixcolor;
    }

    //Most common occurances

    // ## 1 ##
    // Splitting in the red between two walkers
    else if((b->walker1 < LATTICE_SIZE)&&
            (LATTICE_SIZE < b->walker2)&&
            (b->walker1 != 0)){
        obj->tree->at(initsize)->walker1 = b->walker1;
        obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;
        obj->tree->at(initsize)->walker2 = 0;
        obj->tree->at(initsize+1)->walker1 = 0;
    }

    // ## 2 ##
    //Splitting in the green between two walkers
    else if((b->walker1 > LATTICE_SIZE)&&
            (LATTICE_SIZE > b->walker2)&&
            (b->walker2 != 0)){
        obj->tree->at(initsize)->walker2 = b->walker2;
        obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
        obj->tree->at(initsize)->walker1 = 0;
        obj->tree->at(initsize+1)->walker2 = 0;
    }

    // ## 3 ##
    //Right branch fixated green/re
    else if((b->walker1 < LATTICE_SIZE)&&
            (b->walker2 < LATTICE_SIZE)&&
            (b->walker1 > 0)&&
            (b->walker1 < b->walker2)){
        obj->tree->at(initsize+1)->fixation = true;
        obj->tree->at(initsize+1)->fixcolor = 'g';
        obj->tree->at(initsize)->walker1 = b->walker1;
        obj->tree->at(initsize)->walker2 = b->walker2;
    }else if((b->walker1 < LATTICE_SIZE)&&
            (b->walker2 < LATTICE_SIZE)&&
            (b->walker2 > 0)&&
            (b->walker1 > b->walker2)){
        obj->tree->at(initsize+1)->fixation = true;
        obj->tree->at(initsize+1)->fixcolor = 'r';
        obj->tree->at(initsize)->walker1 = b->walker1;
        obj->tree->at(initsize)->walker2 = b->walker2; 
    }

    // ## 4 ##
    //Left branch fixated red/green
    else if((b->walker1 > LATTICE_SIZE)&&
            (b->walker2 > LATTICE_SIZE)&&
            (b->walker1 < b->walker2)){
        obj->tree->at(initsize)->fixation = true;
        obj->tree->at(initsize)->fixcolor = 'g';
        obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
        obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;
    }else if((b->walker1 > LATTICE_SIZE)&&
            (b->walker2 > LATTICE_SIZE)&&
            (b->walker1 > b->walker2)){
        obj->tree->at(initsize)->fixation = true;
        obj->tree->at(initsize)->fixcolor = 'r';
        obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
        obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;
    }

    // ## 5 ##
    //Walker1 = 0 block     3,4,11.
    else if(b->walker1 == 0){
        if(b->walker2 < LATTICE_SIZE){
            //Right branch fixates green
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'g';
            obj->tree->at(initsize)->walker1 = b->walker1;
            obj->tree->at(initsize)->walker2 = b->walker2;        
        }else if(b->walker2 > LATTICE_SIZE){
            //Left branch fixates red
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'r';
            obj->tree->at(initsize+1)->walker1 = b->walker1;
            obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;    
        }else if(b->walker2 == LATTICE_SIZE){
            //Left fixates red AND Right fixates green
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'r';
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'g';    
        }
    }

    // ## 6 ##
    //Walker2 = 0 block     5,6,12.
    else if(b->walker2 == 0){
        if((b->walker1 < LATTICE_SIZE)&&
                (b->walker1 != 0)){
            //Right branch fixates red
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'r';
            obj->tree->at(initsize)->walker1 = b->walker1;
            obj->tree->at(initsize)->walker2 = b->walker2;    
        }else if(b->walker1 > LATTICE_SIZE){
            //Left branch fixates green
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'g';
            obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
            obj->tree->at(initsize+1)->walker2 = b->walker2;    
        }else if(b->walker1 == LATTICE_SIZE){
            //Left fixates green AND Right fixates red
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'g';
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'r';    
        }
    }

    // ## 7 ##
    //Walker1 = LATTICE_SIZE    7,8.
    else if(b->walker1 == LATTICE_SIZE){
        if((b->walker2 < LATTICE_SIZE)&&
                (b->walker2 != 0)){
            //Right branch fixates red
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'r';
            obj->tree->at(initsize)->walker1 = b->walker1;
            obj->tree->at(initsize)->walker2 = b->walker2;    
        }else if(b->walker2 > LATTICE_SIZE){
            //Left branch fixates green
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'g';
            obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
            obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;    
        }
    }

    //Walker2 = LATTICE_SIZE    9,10.
    else if(b->walker2 == LATTICE_SIZE){
        if((b->walker1 < LATTICE_SIZE)&&
                (b->walker1 != 0)){
            //Right branch fixates green
            obj->tree->at(initsize+1)->fixation = true;
            obj->tree->at(initsize+1)->fixcolor = 'g';
            obj->tree->at(initsize)->walker1 = b->walker1;
            obj->tree->at(initsize)->walker2 = b->walker2;

        }else if(b->walker1 > LATTICE_SIZE){
            //Left branch fixates green
            obj->tree->at(initsize)->fixation = true;
            obj->tree->at(initsize)->fixcolor = 'r';
            obj->tree->at(initsize+1)->walker1 = b->walker1 - LATTICE_SIZE;
            obj->tree->at(initsize+1)->walker2 = b->walker2 - LATTICE_SIZE;    
        }
    }
    b->termination = true;
    b->growing = false;

    // Parent branch is terminated, no longer active
    obj->active->erase(b);
    // Cuts down memory useage. These structures get crazy big.
    if(!BARW){ b->surface.clear(); delete b; obj->tree->at(i) = NULL; }
}

// Debug status function
void Status(const VB& tree){
    for(int i=0;i<tree.size();i++){
        cout << "\n\tBranch " << i;
        if(tree[i]->fixation)
            cout << " - FIX";
        if(tree[i]->termination)
            cout << " - TERM";
        cout << "\nBranch ID = \t" << tree[i];
        cout << "\nGeneration # =\t" << tree[i]->surface.size();
        cout << "\nwalker1 = \t" << tree[i]->walker1;
        cout << "\nwalker2 = \t" << tree[i]->walker2;
        cout << "\nprevwalker1 =\t" << tree[i]->prevwalker1;
        cout << "\nprevwalker2 =\t" << tree[i]->prevwalker2;
        cout << "\nfixation bool =\t" << tree[i]->fixation;
        cout << "\nfixcolor =\t" << tree[i]->fixcolor;
        cout << "\ngrowth flag = \t" << tree[i]->growing;
        cout << "\ngrwowth gen = \t" << tree[i]->growgen;
        cout << "\nlattice size =\t" << tree[i]->lattice_size;
        cout << "\nNlattice size =\t" << tree[i]->lattice_size_actual;
        cout << "\nx =\t" << tree[i]->pos_arr[0];
        cout << "\ny =\t" << tree[i]->pos_arr[1];
        cout << "\nz =\t" << tree[i]->pos_arr[2];
        cout << "\nxdir =\t" << tree[i]->dir_arr[0];
        cout << "\nydir =\t" << tree[i]->dir_arr[1];
        cout << "\nzdir =\t" << tree[i]->dir_arr[2];
        if(i != 0){
            for(int j=0;j<tree.size();j++){
                if(tree[j] == tree[i]->parent){
                    cout << "\nparent = \t# " << j;
                    break;
                }
            }
        }
        if(i != 0){
            for(int j=0;j<tree.size();j++){
                if(tree[j] == tree[i]->sister){
                    cout << "\nsister = \t# " << j;
                    break;
                }
            }
        }
        cout << "\nbranch gen =\t" << tree[i]->branchGeneration;
        cout << endl;
    }
}


template <typename T>
double avg(vector<T> numbers){
    double sum;
    if(numbers.size()){
        for(int i=0;i<numbers.size();i++){
            sum += numbers[i];
        }
        return ((double)sum)/numbers.size();
    }else{
        return 0;
    }
}

template <typename T>
double avg(T* numbers,int size){
    double sum = 0;
    if(size){
        for(int i=0; i<size; i++){
            sum += numbers[i];
        }
        return ((double)sum)/size;
    }else{
        return 0;
    }
}

// Sets up a single new branch at the origin (for new systems)
void Reinitialize(thread_obj* obj){
    obj->active->clear();
    for(int i=0;i<obj->tree->size();i++){
        if(obj->tree->at(i) != NULL){
            delete obj->tree->at(i);
        }
    }
    obj->tree->clear();
    branch* latt = new branch;
    latt->lattice_size = LATTICE_SIZE;
    latt->lattice_size_actual = LATTICE_SIZE*1.0;
    obj->tree->push_back(latt);
    obj->active->insert(latt);

    gen genr;
    genr.pos_arr[0] = 0.0;
    genr.pos_arr[1] = 0.0;
    genr.pos_arr[2] = 0.0;
    latt->dir_arr[0] = 0;
    latt->dir_arr[1] = 0;
    latt->dir_arr[2] = 1;
    genr.color = RedCellNumber/LATTICE_SIZE;
    latt->surface.push_back(genr);

    SetUpWalkers(*(obj->tree));
    obj->NUMBER_OF_GENERATIONS = 0;
    obj->Fixation = false;
    obj->Extinction = false;
    obj->Steady_State = false;
    obj->Total_death = false;
}

// Returns true if any branch is in a steady state
bool AreAnyRed(VB& tree){
    for(int i=0;i<tree.size();i++){
        if(tree[i]->fixcolor == 'r'){
            return true;
        }
    }
    return false;
}

double PowerLaw(int t){
    return (lambda * pow(t,beta));
}

// Single cell competition function for growth on an interface
bool minicomp(){
    if(dist(mt)<=redprob){
        return 1;
    }else{
        return 0;
    }
}


// Extension of the power law that defines growth rate * 2rPi
double Length(int growgen){
    return (2*M_PI*radius)*(1+PowerLaw(growgen));
}



/************************************
 * 
 *  ~ Position Functions ~
 * 
 ***********************************/

// 
void InitializePos_Dir(branch& lattice_1,branch& lattice_2,branch* b,rotor& rot){
    double x_1,x_2,x_3;
    vector<double> v;
    x_1 = pidist(mt);
    x_2 = pidist(mt);
    x_3 = alphadist(mt);

    lattice_1.pos_arr[0] = b->pos_arr[0];
    lattice_1.pos_arr[1] = b->pos_arr[1];
    lattice_1.pos_arr[2] = b->pos_arr[2];
    lattice_2.pos_arr[0] = b->pos_arr[0];
    lattice_2.pos_arr[1] = b->pos_arr[1];
    lattice_2.pos_arr[2] = b->pos_arr[2];
    lattice_1.dir_arr[0] = b->dir_arr[0];
    lattice_1.dir_arr[1] = b->dir_arr[1];
    lattice_1.dir_arr[2] = b->dir_arr[2];
    lattice_2.dir_arr[0] = b->dir_arr[0];
    lattice_2.dir_arr[1] = b->dir_arr[1];
    lattice_2.dir_arr[2] = b->dir_arr[2];

    DirMut(lattice_1,rot,x_1,x_2,x_3);
    // Pay attention to the -x_2
    DirMut(lattice_2,rot,x_1,-x_2,x_3);
}


// Fast matrix rotation (from an off axis initial angle)
void DirMut(branch& lat,rotor& rot,double x_1,double x_2,double x_3){
    double theta;
    // Establish a,z and cross v
    // =========================
    rot.a_vec[0] = lat.dir_arr[0];
    rot.a_vec[1] = lat.dir_arr[1];
    rot.a_vec[2] = lat.dir_arr[2];

    norm(rot.a_vec);
    rot.v_vec = cross3(rot.a_vec,rot.z_vec);
    if((mag(rot.v_vec) < 0.000000001)){
        if(rot.a_vec[2] > 0){
            rot.v_vec[0] = 0;
            rot.v_vec[1] = 0;
            rot.v_vec[2] = 1;
        }else{
            rot.v_vec[0] = 0;
            rot.v_vec[1] = 0;
            rot.v_vec[2] = -1;
        }
    }
    norm(rot.v_vec);

    // Build components of R
    // =====================
    rot.sksym = vec_to_SkewSym(rot.v_vec);
    rot.s = mag(rot.v_vec);
    rot.c = dot(rot.a_vec,rot.z_vec);
    theta = acos(rot.c);
    rot.v_x2 = matmul(rot.sksym,rot.sksym);
    s_x_M(sin(theta),rot.sksym);
    s_x_M(2*pow(sin(theta/2),2),rot.v_x2);

    // Build R
    // =======
    rot.R = matadd(rot.ID3,rot.sksym);
    rot.R = matadd(rot.R,rot.v_x2);

    // Apply angular rotations
    // ======================
    Uniform_Sphere(rot.z_vec,x_1,x_2,x_3);
    rot.new_dir[0] = rot.z_vec[0];
    rot.new_dir[1] = rot.z_vec[1];
    rot.new_dir[2] = rot.z_vec[2];

    // Rotate back to "original direction" + rotations
    // ===============================================
    rot.Inv_R = inv3(rot.R);
    rot.a2_vec = matmul(rot.Inv_R,rot.new_dir);
    norm(rot.a2_vec);

    // Update new branch's direciton
    // =============================
    lat.dir_arr[0] = rot.a2_vec[0];
    lat.dir_arr[1] = rot.a2_vec[1];
    lat.dir_arr[2] = rot.a2_vec[2];

    rot.z_vec[0] = 0;
    rot.z_vec[1] = 0;
    rot.z_vec[2] = 1;
}

// Adds position for new tips based on previous position, moving (1) unit in "dir(ection)""
void UpdatePos(thread_obj* obj){
    branch* b;
    for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
        b = *(obj->sit);
        for(int j=0;j<3;j++)
            b->pos_arr[j] += length_unit * b->dir_arr[j];
        gen tip;
        tip.pos_arr[0] = b->pos_arr[0];
        tip.pos_arr[1] = b->pos_arr[1];
        tip.pos_arr[2] = b->pos_arr[2];         
        if(!b->fixation){
            if(b->walker2 > b->walker1){
                tip.color =  1.0*abs(b->walker2 - b->walker1)/b->lattice_size;
            }else if(b->walker1 > b->walker2){
                tip.color = 1 - (1.0*abs(b->walker2 - b->walker1)/b->lattice_size);
            }
        }else{
            b->fixcolor == 'r' ? tip.color = 1.0 : tip.color = 0.0;
        }
        if(tip.pos_arr[0] < b->pos_wall[0]) b->pos_wall[0] = tip.pos_arr[0];
        if(tip.pos_arr[0] > b->pos_wall[1]) b->pos_wall[1] = tip.pos_arr[0];
        if(tip.pos_arr[1] < b->pos_wall[2]) b->pos_wall[2] = tip.pos_arr[1];
        if(tip.pos_arr[1] > b->pos_wall[3]) b->pos_wall[3] = tip.pos_arr[1];
        if(tip.pos_arr[2] < b->pos_wall[4]) b->pos_wall[4] = tip.pos_arr[2];
        if(tip.pos_arr[2] > b->pos_wall[5]) b->pos_wall[5] = tip.pos_arr[2];

        b->surface.push_back(tip);    // Adds new tip cell to the tree
    }
    return;
}

// See paper for angle distribution reasons
void RotateDir(thread_obj* obj){
    double x_1,x_2,x_3;
    // Angle of rotation = pi/10
    for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
        x_1 = pidist(mt);
        x_2 = pidist(mt);
        x_3 = d2ddist(mt);
        DirMut(*(*obj->sit),obj->rot,x_1,x_2,x_3);
    }
}

// Writes position data to files for model-making
void RecordPosData(const VB& tree,ofstream& clout,int NUMBER_OF_GENERATIONS){
    for(int i=0;i<tree.size();i++){
        for(int j=0;j<tree[i]->surface.size();j++){
            clout   << tree[i]->surface[j].pos_arr[0] << "\t"
                << tree[i]->surface[j].pos_arr[1] << "\t"
                << tree[i]->surface[j].pos_arr[2] << "\t"
                << tree[i]->surface[j].color << "\t";
            clout << "<! B:" << i << " C: " << j << " >\n";
        }
    }
}


/* The code for collision detection between each tip.
 * 
 * For each tip, every generation within the dimensions of the branch itself,
 * like a cube with verticies at (0 - r,0 - r,0 - r) and (10 + r,10 + r,10 + r)
 * is checked for collision. 
 * 
 * 
 */
void CollisionDetection(VB& tree,set<branch*>* active){
    double self_ann_limit = LATTICE_SIZE/2.5;
    for(unsigned int i=0;i<tree.size();i++){
        if(!tree[i]->termination){
            for(unsigned int j=0;j<tree.size();j++){
                if(j != i){
                    // Large, massively confusing flag
                    // Means "if the generation < 21 (for 100), don't check for
                    // collision on sister or parent branches."
                    if(!((tree[i]->generation < 2*self_ann_limit) &&
                                ((tree[j] == tree[i]->parent) || (tree[j] == tree[i]->sister)))){   
                        if((tree[i]->surface.back().pos_arr[0] > tree[j]->pos_wall[0] - radius) &&
                                (tree[i]->surface.back().pos_arr[0] < tree[j]->pos_wall[1] + radius) &&
                                (tree[i]->surface.back().pos_arr[1] > tree[j]->pos_wall[2] - radius) &&
                                (tree[i]->surface.back().pos_arr[1] < tree[j]->pos_wall[3] + radius) &&
                                (tree[i]->surface.back().pos_arr[2] > tree[j]->pos_wall[4] - radius) &&
                                (tree[i]->surface.back().pos_arr[2] < tree[j]->pos_wall[5] + radius)){
                            for(int k=tree[j]->surface.size()-1;k>=0;k--){
                                if(!tree[i]->termination){
                                    tree[i]->termination = 
                                        sqrt(pow(abs(tree[i]->surface.back().pos_arr[0] - tree[j]->surface[k].pos_arr[0]),2) +
                                                pow(abs(tree[i]->surface.back().pos_arr[1] - tree[j]->surface[k].pos_arr[1]),2) +
                                                pow(abs(tree[i]->surface.back().pos_arr[2] - tree[j]->surface[k].pos_arr[2]),2))
                                        < tree[i]->radius ? 1 : 0;
                                    if(tree[i]->termination){

                                        annihilations++;
                                        tree[i]->growing = false;
                                        active->erase(tree[i]);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }else{
                    // Self-annihilation
                    if(tree[i]->surface.size() > self_ann_limit){
                        for(int k=0;k<tree[i]->surface.size()-self_ann_limit;k++){
                            tree[i]->termination =
                                sqrt(pow(abs(tree[i]->surface.back().pos_arr[0] - tree[i]->surface[k].pos_arr[0]),2) +
                                        pow(abs(tree[i]->surface.back().pos_arr[1] - tree[i]->surface[k].pos_arr[1]),2) +
                                        pow(abs(tree[i]->surface.back().pos_arr[2] - tree[i]->surface[k].pos_arr[2]),2))
                                < tree[i]->radius ? 1 : 0;
                            if(tree[i]->termination){
                                //printf("Branch %i collision with Branch %i, cell %i\n",i,j,k);
                                annihilations++;
                                tree[i]->growing = false;
                                active->erase(tree[i]);
                                break;

                            }
                        }
                    }
                }
                if(tree[i]->termination) break;
            }
        }
    }
    //cout << "Exiting Collision\n";
    return;
}
