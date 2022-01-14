/* Code author: Adam Bryant
 * summer 2019
 * 
 * The header that holds every important function for the BARW cell lattice model.
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
int LATTICE_SIZE = 20;
int RedCellNumber = 1;
double lambda = 0.005;
double beta = 2;
volatile double redprob;
volatile double rate;
double length_unit = 1;
short annihilations;
double radius;
int NUMBER_OF_THREADS = 1;
int NUMBER_OF_STRUCTURES = 100;
int WORKLOAD = NUMBER_OF_STRUCTURES/NUMBER_OF_THREADS;
int MAX_GENERATIONS = 3000;
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
std::uniform_real_distribution<double> d2ddist(0.0,(1 - cos((2*M_PI)/double(LATTICE_SIZE)))/2);


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
    vector<vector<double>>* Dirs;
    vector<double>* Perlens;
    int ID;
    int attempts;
    int NUMBER_OF_GENERATIONS;
};

// Lattice competition functions
bool CellEvolution(thread_obj*,int);
void Competition(thread_obj*);
void Competition2(thread_obj*);
void CheckForFixation(VB&);
void SetUpWalkers(VB&);
void Status(const VB&);

// CHANGE TO POINTERS :/
double avg(const vector<double>&);
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
 i 
 *  ~ Function Definitions ~
 * 
 ******************************/
void* ThreadSpinner(void* t){
    thread_obj* obj = (thread_obj*) t;
    short attempts = 0;
    int i,j,k;
    int mult;
    vector<double> ToBeAvg;
    ToBeAvg.resize(0);
    ofstream fout;
    fout.open("pos.vtk");

    obj->Perlens->clear();
    obj->Perlens->resize(MAX_GENERATIONS/3);

    fill(obj->Perlens->begin(),obj->Perlens->end(),double(0));

    obj->Dirs->clear();
    obj->Dirs->resize(MAX_GENERATIONS);
    vector<double> d = {0,0,0};
    fill(obj->Dirs->begin(),obj->Dirs->end(),d);

    for(i = 0; i < WORKLOAD; i++){
        if(CellEvolution(obj,i)){ ; }
        printf("%d\n",i);
        if(!i){
            for(j = 0; j < MAX_GENERATIONS; j++){
                fout << obj->Dirs->at(j)[0] << " "
                << obj->Dirs->at(j)[1] << " "
                << obj->Dirs->at(j)[2] << "\n";
            }
            fout.close();
        }
        for(j = 1; j <= MAX_GENERATIONS/3; j++){
            mult = j;
            for(k = 0; k < MAX_GENERATIONS - mult; k++){
                ToBeAvg.push_back(dot(obj->Dirs->at(k),obj->Dirs->at(k+mult)));
            }
            obj->Perlens->at(j-1) += avg(ToBeAvg);
            ToBeAvg.clear();
        }

        Reinitialize(obj);
    }

    double div = double(WORKLOAD);
    for(i = 0; i < MAX_GENERATIONS/3; i++)
        obj->Perlens->at(i) /= div;

    return (void*) obj;
}

void ClearObject(thread_obj* o){
    
}
//The main program within the main program
bool CellEvolution(thread_obj* obj,int index){
    int i;
    branch* b;
    int ss_time = 0;
    while(obj->NUMBER_OF_GENERATIONS < MAX_GENERATIONS){
        
        obj->Dirs->at(obj->NUMBER_OF_GENERATIONS)[0] = obj->tree->at(0)->dir_arr[0];
        obj->Dirs->at(obj->NUMBER_OF_GENERATIONS)[1] = obj->tree->at(0)->dir_arr[1];
        obj->Dirs->at(obj->NUMBER_OF_GENERATIONS)[2] = obj->tree->at(0)->dir_arr[2];

        
        obj->NUMBER_OF_GENERATIONS++;
        // Exit conditions


        RotateDir(obj);
    }
    return true;
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

    obj->NUMBER_OF_GENERATIONS = 0;
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

double avg(const vector<double>& a){
    double sum = 0;
    for(unsigned int i = 0; i < a.size(); i++){
        sum += a[i];
    }
    return sum/double(a.size());
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
            }
        }
    }
    //cout << "Exiting Collision\n";
    return;
}
