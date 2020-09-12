#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <iterator>

using namespace std;

typedef vector<vector <double> > M;

vector<double> arr_to_vec(double* a,int size){
    vector<double> vec(size);
    for(unsigned int i=0;i<size;i++) vec[i] = a[i];
    return vec;
}



double dot(const double* a,const double* b,int size){    
    double sum = 0;
    for(unsigned int i=0;i<size;i++)
        sum += a[i] * b[i];
    return sum;       
}

float dot(const float* a,const float* b,int size){    
    float sum = 0;
    for(unsigned int i=0;i<size;i++)
        sum += a[i] * b[i];
    return sum;       
}

double dot(const vector<double>& a,const vector<double>& b){
    if(a.size() != b.size())
        return -1;
    
    double sum = 0;
    for(unsigned int i=0;i<a.size();i++)
        sum += a[i] * b[i];
    return sum; 
}

/*
double* cross3(const double* a,const double* b){
    if(sizeof(a) != sizeof(b))
        return NULL;
    double c[3];
    c[0] = (a[1]*b[2]) - (a[2]*b[1]);
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);
    c[2] = (a[0]*b[1]) - (a[1]*b[0]);
    return c;
}
*/



vector<double> cross3(const vector<double> a,const vector<double> b){
    vector<double> c(3);
    if(a.size() != b.size())
        return c;
    c[0] = (a[1]*b[2]) - (a[2]*b[1]);
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);
    c[2] = (a[0]*b[1]) - (a[1]*b[0]);
    return c;
}

/*
M matmul(const double** a,const double** b){
    if(sizeof(a) != sizeof(b[0])){
        M mat;
        return mat;
    }
    unsigned int w,h;
    w = sizeof(a)/sizeof(a[0][0]);
    h = sizeof(b[0])/sizeof(b[0][0]);
    M mat(w,vector<double>(h));
    double sum;
    for(unsigned int i=0;i<w;i++){
        for(unsigned int j=0;j<h;j++){
            sum = 0;
            for(unsigned int k=0;k<h;k++){
                sum += a[i][k] * b[k][j];
            }
            mat[i][j] = sum;
        }
    }
    return mat;
}
*/

M matmul(const M& a,const M& b){
    if(a.size() != b[0].size()){
        M mat;
        return mat;
    }
    unsigned int w,h;
    w = a.size();
    h = b[0].size();
    M mat(w,vector<double>(h));
    double sum;
    for(unsigned int i=0;i<w;i++){
        for(unsigned int j=0;j<h;j++){
            sum = 0;
            for(unsigned int k=0;k<h;k++){
                sum += a[i][k] * b[k][j];
            }
            mat[i][j] = sum;
        }
    }
    return mat;
}

vector<double> matmul(const M& a,const vector<double>& b){
    if(a.size() != b.size()){
        vector<double> mat;
        return mat;
    }
    unsigned int w,h;
    w = a.size();
    h = b.size();
    vector<double> v(w);
    double sum;
    for(unsigned int i=0;i<w;i++){
        sum = 0;
        for(unsigned int j=0;j<h;j++){
            sum += a[i][j] * b[j];
        }
        v[i] = sum;
    }
    return v;
}

// matmul for v^T * v = M
M matmul(const M& a,const vector<double>& b,bool vtv){
    unsigned int w,h;
    w = a.size();
    h = b.size();    
    if(w != h){
        M mat;
        return mat;
    }
    M mat(w,vector<double>(h));

    for(unsigned int i=0;i<w;i++){
        for(unsigned int j=0;j<h;j++){
            mat[i][j] = a[i][0] * b[j];
        }
    }
    return mat;
}

/*
M matadd(const double** a,const double** b){
    unsigned int w,h;
    w = sizeof(a)/sizeof(a[0][0]);
    h = sizeof(a[0])/sizeof(a[0][0]);
    M mat(w,vector<double>(h));
    if((sizeof(a) != sizeof(b)) || (sizeof(a[0]) != sizeof(b[0])))
        return mat;
    
    for(unsigned int i=0;i<w;i++)
        for(unsigned int j=0;j<h;j++)
            mat[i][j] = a[i][j] + b[i][j];
    return mat;
}
*/

M matadd(const M& a,const M& b){
    unsigned int w,h;
    w = a.size();
    h = a[0].size();
    M mat(w,vector<double>(h));
    if((a.size() != b.size()) || (a[0].size() != b[0].size()))
        return mat;
    
    for(unsigned int i=0;i<w;i++)
        for(unsigned int j=0;j<h;j++)
            mat[i][j] = a[i][j] + b[i][j];
    return mat;
}

vector<double> vecadd(const vector<double>& a,const vector<double>& b){
    vector<double> c(a.size());
    if(a.size() != b.size())
        return c;
    for(unsigned int i=0;i<a.size();i++)
        c[i] = a[i] + b[i];
    return c;
}

/*
double** I(int n){
    if(n <= 0)
        return NULL;
    double ID[n][n];
    for(unsigned int i=0;i<n;i++)
        for(unsigned int j=0;j<n;j++)
            i == j ? ID[i][j] = 1 : ID[i][j] = 0;
    
    return ID;
}
*/

M I(int n){
    M ID(n,vector<double>(n));
    if(n <= 0)
        return ID;
    for(unsigned int i=0;i<n;i++)
        for(unsigned int j=0;j<n;j++)
            i == j ? ID[i][j] = 1 : ID[i][j] = 0;
    return ID;
}


double mag(double* a,int size){
    double sum = 0;
    if(size > 0){
        for(unsigned int i=0;i<size;i++){
            sum += pow(a[i],2);
        }
        return sqrt(sum);
    }else{
        return 0;
    }
}


double mag(const vector<double>& a){
    double sum = 0;
    if(a.size() > 0){
        for(unsigned int i=0;i<a.size();i++)
            sum += pow(a[i],2);
        return sqrt(sum);
    }else{
        return 0;
    }
}


void norm(double* a,int size){
    short i=0;
    double m = mag(a,size);
    while(m != 1.0 && i < 10){
        for(unsigned int i=0;i<size;i++)
            a[i] /= m;
        m = mag(a,size);
        i++;
    }
}


void norm(vector<double>& a){
    short i=0;
    double m = mag(a);
    while(m != 1.0 && i < 10){
        for(unsigned int i=0;i<a.size();i++)
            a[i] /= m;
        m = mag(a);
        i++;
    }
}

void XMatRotation(double phi,double arr[3]){
    arr[1] +=   (cos(phi)*arr[1]) - (sin(phi)*arr[2]);
    arr[2] +=   (sin(phi)*arr[1]) + (cos(phi)*arr[2]);
    norm(arr,3);
}

void XMatRotation(double phi,vector<double>& arr){
    arr[1] +=   (cos(phi)*arr[1]) - (sin(phi)*arr[2]);
    arr[2] +=   (sin(phi)*arr[1]) + (cos(phi)*arr[2]);
    norm(arr);
}

void YMatRotation(double theta,double arr[3]){
    arr[0] +=   (cos(theta)*arr[0]) + (sin(theta)*arr[2]);
    arr[2] +=  -(sin(theta)*arr[0]) + (cos(theta)*arr[2]);
    norm(arr,3);
}

void YMatRotation(double theta,vector<double>& arr){
    arr[0] +=   (cos(theta)*arr[0]) + (sin(theta)*arr[2]);
    arr[2] +=  -(sin(theta)*arr[0]) + (cos(theta)*arr[2]);
    norm(arr);
}

void ZMatRotation(double psi,double arr[3]){
    arr[0] +=   (cos(psi)*arr[0]) - (sin(psi)*arr[1]);
    arr[1] +=   (sin(psi)*arr[0]) + (cos(psi)*arr[1]);
    norm(arr,3);
}

void ZMatRotation(double psi,vector<double>& arr){
    arr[0] +=   (cos(psi)*arr[0]) - (sin(psi)*arr[1]);
    arr[1] +=   (sin(psi)*arr[0]) + (cos(psi)*arr[1]);
    norm(arr);
}

void MatRotation(double arr[3],double a,double b,double c){
    double dir0,dir1,dir2;
    dir0 = arr[0];
    dir1 = arr[1];
    dir2 = arr[2];
    arr[0] += ((cos(a)*cos(b))*dir0) +
    ((cos(a)*sin(b)*sin(c) - sin(a)*cos(c))*dir1) +
    ((cos(a)*sin(b)*cos(c) + sin(a)*sin(c))*dir2);
    arr[1] += ((sin(a)*cos(b))*dir0) +
    ((sin(a)*sin(b)*sin(c) + cos(a)*cos(c))*dir1) +
    ((sin(a)*sin(b)*cos(c) - cos(a)*sin(c))*dir2);
    arr[2] += ((-sin(b))*dir0) +
    ((cos(b)*sin(c))*dir1) +
    ((cos(b)*cos(c))*dir2);

    norm(arr,3);
    return;
}

void MatRotation(vector<double>& arr,double a,double b,double c){
    double dir0,dir1,dir2;
    dir0 = arr[0];
    dir1 = arr[1];
    dir2 = arr[2];
    arr[0] += ((cos(a)*cos(b))*dir0) +
    ((cos(a)*sin(b)*sin(c) - sin(a)*cos(c))*dir1) +
    ((cos(a)*sin(b)*cos(c) + sin(a)*sin(c))*dir2);
    arr[1] += ((sin(a)*cos(b))*dir0) +
    ((sin(a)*sin(b)*sin(c) + cos(a)*cos(c))*dir1) +
    ((sin(a)*sin(b)*cos(c) - cos(a)*sin(c))*dir2);
    arr[2] += ((-sin(b))*dir0) +
    ((cos(b)*sin(c))*dir1) +
    ((cos(b)*cos(c))*dir2);

    norm(arr);
    return;
}

void MatRotation_Sphere(double arr[3],double a,double b,double c){
    double dir0,dir1,dir2;
    dir0 = arr[0];
    dir1 = arr[1];
    dir2 = arr[2];
    arr[0] = ((cos(a)*cos(b))*dir0) +
    ((cos(a)*sin(b)*sin(c) - sin(a)*cos(c))*dir1) +
    ((cos(a)*sin(b)*cos(c) + sin(a)*sin(c))*dir2);
    arr[1] = ((sin(a)*cos(b))*dir0) +
    ((sin(a)*sin(b)*sin(c) + cos(a)*cos(c))*dir1) +
    ((sin(a)*sin(b)*cos(c) - cos(a)*sin(c))*dir2);
    arr[2] = ((-sin(b))*dir0) +
    ((cos(b)*sin(c))*dir1) +
    ((cos(b)*cos(c))*dir2);

    norm(arr,3);
    return;
}

void MatRotation_Sphere(vector<double>& arr,double a,double b,double c){
    double dir0,dir1,dir2;
    dir0 = arr[0];
    dir1 = arr[1];
    dir2 = arr[2];
    arr[0] = ((cos(a)*cos(b))*dir0) +
    ((cos(a)*sin(b)*sin(c) - sin(a)*cos(c))*dir1) +
    ((cos(a)*sin(b)*cos(c) + sin(a)*sin(c))*dir2);
    arr[1] = ((sin(a)*cos(b))*dir0) +
    ((sin(a)*sin(b)*sin(c) + cos(a)*cos(c))*dir1) +
    ((sin(a)*sin(b)*cos(c) - cos(a)*sin(c))*dir2);
    arr[2] = ((-sin(b))*dir0) +
    ((cos(b)*sin(c))*dir1) +
    ((cos(b)*cos(c))*dir2);

    norm(arr);
    return;
}

M T(const M& a){
    unsigned int w,h;
    w = a.size();
    h = a[0].size();
    M mat(w,vector<double>(h));
    for(unsigned int i=0;i<w;i++)
        for(unsigned int j=0;j<h;j++)
            mat[i][j] = a[j][i];
    return mat;
}

M T(const vector<double>& a){
    unsigned int w,h;
    w = a.size();
    M mat(w,vector<double>(1));
    for(unsigned int i=0;i<w;i++)
        mat[i][0] = a[i];
    return mat;
}

M vec_to_SkewSym(const vector<double>& a){
    M mat(3,vector<double>(3));
    mat[0][0] = 0;
    mat[0][1] = -a[2];
    mat[0][2] = a[1];
    mat[1][0] = a[2];
    mat[1][1] = 0;
    mat[1][2] = -a[0];
    mat[2][0] = -a[1];
    mat[2][1] = a[0];
    mat[2][2] = 0;
    return mat;
}

double det3(const double a[3][3]){
    double d =  (a[0][0]*((a[1][1]*a[2][2]) - (a[1][2]*a[2][1]))) -
                (a[0][1]*((a[1][0]*a[2][2]) - (a[1][2]*a[2][0]))) +
                (a[0][2]*((a[1][0]*a[2][1]) - (a[1][1]*a[2][0])));
    return d;
}

double det3(const M& a){
    if(a.size() != 3 || a[0].size() != 3)
        return 0;
    double d =  (a[0][0]*((a[1][1]*a[2][2]) - (a[1][2]*a[2][1]))) -
                (a[0][1]*((a[1][0]*a[2][2]) - (a[1][2]*a[2][0]))) +
                (a[0][2]*((a[1][0]*a[2][1]) - (a[1][1]*a[2][0])));
    return d;
}


double det2(const double a[2][2]){
    return (a[0][0]*a[1][1]) - (a[0][1]*a[1][0]);
}


double det2(const M& a){
    if(a.size() != 2 || a[0].size() != 2)
        return 0;
    return (a[0][0]*a[1][1]) - (a[0][1]*a[1][0]);
}

/*
void s_x_M(double s,double** a){
    for(unsigned int i=0;i<sizeof(a)/sizeof(a[0][0]);i++)
        for(unsigned int j=0;j<sizeof(a[0])/sizeof(a[0][0]);j++)
            a[i][j] *= s;
}
*/

void s_x_M(double s,M& a){
    for(unsigned int i=0;i<a.size();i++)
        for(unsigned int j=0;j<a[0].size();j++)
            a[i][j] *= s;
}

void s_x_v(double s,vector<double>& a){
    for(unsigned int i=0;i<a.size();i++)
        a[i] *= s;
}


M inv3(const M& a){
    double det = det3(a);
    if(det == 0){
        cerr << "No inverse\n";
        exit(EXIT_FAILURE);
    }
    M mat(3,vector<double>(3));
    M t = T(a);
    // Create the Adjugate Matrix
    vector<M> minors;
    M dets(3,vector<double>(3));
    for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            M mini(2,vector<double>(2));
            vector<double> nums;
            for(unsigned int q=0;q<3;q++){
                for(unsigned int r=0;r<3;r++){
                    if((q != i) && (r != j)){
                        nums.push_back(t[q][r]);
                    }
                }
            }
            mini[0][0] = nums[0];
            mini[0][1] = nums[1];
            mini[1][0] = nums[2];
            mini[1][1] = nums[3];
            dets[i][j] = det2(mini);
        }
    }

    // Cofactors
    dets[0][1] *= -1;
    dets[1][0] *= -1;
    dets[1][2] *= -1;
    dets[2][1] *= -1;

    s_x_M(1/det,dets);
    return dets;
}

void Uniform_Sphere(vector<double>& dir,double x_1,double x_2,double x_3){
    double theta    = x_1;
    double phi      = x_2;
    double z        = x_3;

    vector<double> v(3);
    v[0] = cos(phi)*sqrt(z);
    v[1] = sin(phi)*sqrt(z);
    v[2] = sqrt(1 - z);
    
    M ID3 = I(3);
    M t = T(v);
    M vvt = matmul(T(v),v,1);

    s_x_M(2,vvt);
    s_x_M(-1,ID3);
    M H = matadd(vvt,ID3);

    M R = I(3);
    R[0][0] = cos(theta);
    R[0][1] = sin(theta);
    R[1][0] = -sin(theta);
    R[1][1] = R[0][0];

    M rotation_mat = matmul(H,R);
    dir = matmul(rotation_mat,dir);

    return;
}

double Edistance(float* a,float* b){
    return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2));
}
