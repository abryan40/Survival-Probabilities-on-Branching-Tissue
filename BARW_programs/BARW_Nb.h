/* Code author: Adam Bryant
 * summer 202-
 * 
 * The header that holds every important function for the BARW cell lattice model.
 * Modified for collection only of BARW statistics, forgoing intercellular competition
 * for statistics relating to the structure itself.
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
double radius;
int NUMBER_OF_THREADS = 50;
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
std::uniform_real_distribution<double> alphadist    ((1-cos((100/3)*(M_PI/180.0)))/2,
        ((1-cos((200/3)*(M_PI/180.0))))/2);
std::uniform_real_distribution<double> pidist(0.0,2.0*M_PI);
std::uniform_real_distribution<double> d2ddist(0.0,(1.0 - cos(2.0*M_PI/double(LATTICE_SIZE)))/2.0);
//Ignore everything intelliSense tries to tell you

// These are only used when making 3D models, not for data collection
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
    bool is_parent;
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
        is_parent = false;
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
    volatile long Active_branches;
    volatile long Annihilations;
    volatile long Total_branches;
    volatile long Non_parent_branches;
    volatile int Total_death;
    int NUMBER_OF_GENERATIONS;
    int ID;
};

// Lattice competition functions
bool CellEvolution(thread_obj*,int&);
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
void CollisionDetection(VB&,set<branch*>*,int&);

/*******************************
 * 
 *  ~ Function Definitions ~
 * 
 ******************************/
void* ThreadSpinner(void* t){
    thread_obj* obj = (thread_obj*) t;
    int np = 0;
    int junk = 0;
    int annihilations = 0;
    obj->Non_parent_branches = 0;
    obj->Active_branches = 0;
    obj->Total_branches = 0;
    obj->Annihilations = 0;
    obj->Total_death = 0;
    for(int i = 0; i < WORKLOAD; i++){
        if(!CellEvolution(obj,annihilations)){
            obj->Total_death++;
        }

        obj->Active_branches += obj->active->size();
        obj->Annihilations += annihilations;
        obj->Total_branches += obj->tree->size();
        for(int j = 0; j < obj->tree->size(); j++){
            if(obj->tree->at(j) != NULL && !(obj->tree->at(j)->is_parent)){
                obj->Non_parent_branches++;
            }
        }
        annihilations = 0;
        //Status(obj->tree);
        Reinitialize(obj);
    }
    return (void*) obj;
}

void ClearObject(thread_obj* o){
    o->Active_branches = 0;
    o->Annihilations = 0;
    o->Total_branches = 0;
    o->Non_parent_branches = 0;  
    o->Total_death = 0;
    o->NUMBER_OF_GENERATIONS = 0;
    Reinitialize(o);  
}
//The main program within the main program
bool CellEvolution(thread_obj* obj,int& annihilations){
    branch* b;
    while(obj->NUMBER_OF_GENERATIONS < MAX_GENERATIONS){
        obj->NUMBER_OF_GENERATIONS++;

        // Exit conditions
        if(AreTheyAllTerminated(obj)){
            return false;
        }

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

        UpdatePos(obj);        
        CollisionDetection(*(obj->tree),obj->active,annihilations);
        RotateDir(obj);

        // Growth and bifurcation loops
        for(obj->sit=obj->active->begin();obj->sit!=obj->active->end();obj->sit++){
            b = *(obj->sit);
            if(b->growing && !(b->termination)){
                if(b->lattice_size >= 2*LATTICE_SIZE){
                    b->growing = false;
                    BranchProcess(obj,*(obj->sit),0);

                    if(obj->active->empty()) break;
                    obj->sit = obj->active->begin();
                }
            }
        }
    }
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
    if(branch.lattice_size_actual > (2*LATTICE_SIZE)){
        newcells = 2*LATTICE_SIZE - branch.lattice_size;
    }
    branch.lattice_size += newcells;
    branch.growgen++;
    branch.radius = double(branch.lattice_size)/(2*M_PI);
    // Adds the cells but we don't care about them
}

// Adds new branches in place of the previous tip
void BranchProcess(thread_obj* obj,branch* b,int i){
    branch* newcells1 = new branch;
    branch* newcells2 = new branch;
    obj->active->insert(newcells1);
    obj->active->insert(newcells2);

    InitializePos_Dir(*newcells1,*newcells2,b,obj->rot);
    obj->tree->push_back(newcells1);
    obj->tree->push_back(newcells2);
    newcells1->branchGeneration = b->branchGeneration + 1;
    newcells2->branchGeneration = b->branchGeneration + 1;
    newcells1->sister = newcells2;
    newcells2->sister = newcells1;
    newcells1->parent = b;
    newcells2->parent = b;
    b->is_parent = true;
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
    latt->lattice_size_actual = double(LATTICE_SIZE);
    latt->radius = double(LATTICE_SIZE)/(2*M_PI);
    obj->tree->push_back(latt);
    obj->active->insert(latt);

    gen genr;
    genr.pos_arr[0] = 0.0;
    genr.pos_arr[1] = 0.0;
    genr.pos_arr[2] = 0.0;
    latt->dir_arr[0] = 0;
    latt->dir_arr[1] = 0;
    latt->dir_arr[2] = 1;
    genr.color = float(RedCellNumber)/LATTICE_SIZE;
    latt->surface.push_back(genr);

    SetUpWalkers(*(obj->tree));
    obj->NUMBER_OF_GENERATIONS = 0;
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
    if(abs(lat.dir_arr[2] - 1) < 0.0000001){
        rot.v_vec[0] = 0;
        rot.v_vec[1] = 0;
        rot.v_vec[2] = 0;
    }else{
        rot.v_vec = cross3(rot.a_vec,rot.z_vec);
        norm(rot.v_vec);
    }
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
void CollisionDetection(VB& tree,set<branch*>* active,int& ann){

    int self_ann_limit = LATTICE_SIZE/2.5;
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
                                        
                                        ann++;
                                        k = 0;
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
                                
                                ann++;
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
