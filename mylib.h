#ifndef MYLIB_H
#define MYLIB_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
//#include <cmath>
#include <vector>
#include <deque>
#include <ctime>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <float.h>
#include <math.h>
//time counting
#include <stdio.h>
#include <time.h>

#define OVERLAP 0
#define NO_OVERLAP 1

#define TRUE 1
#define FALSE 0

#define PI 3.1415926535
#define XI 0.1

#define MAX_SPHERES 200
#define NSWEEP 100 //5000
#define NSWEEP_SW 50 //5000
#define NSWEEP_SW_IN 10000
#define NSWEEP_LATTICE 3000


#define NSWEEP_SCY_IN 2000
#define NSWEEP_SCY 500

#define IRATIO 100

#define REALIZ 100 //500
#define NUM_SHELL 1
#define EPSLON_SHELL 0.001
#define PHI_F 0.49 // hard spheres freezing packing fraction 
#define PHI_P 0.64 // hard spheres maximum packing fraction 

#define NBIN_GR 60
#define DR_GR 0.1

#define NBIN_GR_ANGLE 70
#define DR_GR_ANGLE 0.1

#define NBIN_MINDIST 70
#define DR_MINDIST 0.1

#define NBIN_HP 300
#define DR_HP 0.01
 
#define EPS_CGM 1e-10
#define ITER_MAX 1000
#define MAX_DECIMATED 1


using namespace std;

/* Macro definitions for integer arguments only */


#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(x,y)        ((x) > (y) ? (x) : (y))
#define DEG_TO_RAD(x)   ((double) (x) * M_PI / 180.0)
#define RAD_TO_DEG(x)   ((double) (x) * 180.0 / M_PI)
#define MIN_RIJ(x, L) ( (fabs(x)<=L/2)?x:(x-((x >0)?L:-L) ) )


#define MakeVector(x, y, z, v)          (v)[0]=(x),(v)[1]=(y),(v)[2]=(z)

#define VecCopy(a,b)    ((b)[0] = (a)[0], (b)[1] = (a)[1],\
                         (b)[2] = (a)[2])

#define VecScale(s,a,b)  ((b)[0] = (s) * (a)[0], (b)[1] = (s) * (a)[1],\
                         (b)[2] = (s) * (a)[2])
                         
#define VecAdd(a,b,c)    ((c)[0] = (a)[0] + (b)[0],\
                         (c)[1] = (a)[1] + (b)[1],\
                         (c)[2] = (a)[2] + (b)[2])

#define VecSub(a,b,c)    ((c)[0] = (a)[0] - (b)[0],\
                         (c)[1] = (a)[1] - (b)[1],\
                         (c)[2] = (a)[2] - (b)[2])

#define VecDot(a,b)      (double) ((a)[0] * (b)[0] +\
                         (a)[1] * (b)[1] +\
                         (a)[2] * (b)[2])

#define VecNegate(a, b)  ((b)[0] = - (a)[0], (b)[1] = - (a)[1],\
                         (b)[2] = - (a)[2])

#define VecAddS(s,a,b,c) ((c)[0] = (s) * (a)[0] +(b)[0],\
                         (c)[1] = (s) * (a)[1] +(b)[1],\
                         (c)[2] = (s) * (a)[2] +(b)[2])

#define VecCross(a,b,c)  ((c)[0] = (a)[1] * (b)[2] -(a)[2] * (b)[1],\
                         (c)[1] = (a)[2] * (b)[0] -(a)[0] * (b)[2],\
                         (c)[2] = (a)[0] * (b)[1] -(a)[1] * (b)[0])

#define VecLength(s)     sqrt(VecDot((s), (s)))

#define VecDist(a, b)     sqrt(VecDot((a), (b)))

#define PW2(x) (x*x)

#define pivot_index() (begin+(end-begin)/2)
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))


//physical structures creation functions

void penetrable_spheres(double shell, double diam, long int Nsphere, double Lato);

void impenetrable_spheres (double shellin, double shellfin, int numshell, double *percVect ,double diam, long int Nsphere, double Lato, double * conductance, double Tstar, double lambda);

void place_impenetrable_sphere(double shell,double diam, long int Nsphere, double Lato,  double **pos, int * NpartInSquare, int ** PartInSquarePos);

void place_impenetrable_lattice(double shell,double diam, long int Nsphere, double Lato, double **pos, int * NpartInSquare, int ** PartInSquarePos);

void MC_impenetrable_sphere(double shell,double diam, long int Nsphere, double Lato,  double **pos, int * NpartInSquare, int ** PartInSquarePos, int **connect, double Tstar, double lambda, int sweep, double &passo);

void MC_impenetrable_pure(double shell,double diam, long int Nsphere, double Lato,  double **pos, int * NpartInSquare, int ** PartInSquarePos, int **connect);

//segregation functions
void segregation (double shellin, double shellfin, int numshell, double *percVect, double diamCond, double diamIns, long int NsphereCond,long int NsphereIns, double Lato);

void place_segregation(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns);

int place_segregation_lattice(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns);

void MC_segregation(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns, int **connect);

double adaptMC(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int * NpartInSquareCond, int ** PartInSquarePosCond, std::ofstream& nameFile);

int **createAdj(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int **connect);

void createAdj2(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int **connect, int **AdjList);


//spherocylinder function 
void impenetrable_spherocylinder( int numshell, double *percVect ,double diam, double length, long int Nsphere, double Lato, double *condVect, double Tstar, double lambda, string suffix);

void place_impenetrable_spherocylinder(double shell,double diam, double length, long int Nsphere, double Lato, double **pos, double **direct, int * NpartInSquare, int ** PartInSquarePos);

void MC_spherocyl_pure(double diam, long int Nsphere, double Lato, double length, double **pos, double **direct, int * NpartInSquare, int ** PartInSquarePos, int nsweep, double &passo, double &ds);

void MC_spherocyl_SW(double diam, long int Nsphere, double Lato, double length, double **pos, double **direct, int * NpartInSquare, int ** PartInSquarePos, int Nsweep_cyl, double &pas, double &ds, double Tstar, double lambda );

double dist2_rods(double normr12, double *r12, double *w1, double *w2,double lh1, double lh2);

void flip ( double *wold,  double dotstep);

int **  createAdjcylinder(double shell, double diamCond, double length, int NsphereCond, double Lato, double **posCond, double **direct, int **connect);


//time functions
void getTime (clock_t start);
clock_t start_time();

//memory allocation functions
double * allocate_double_vector(int n);
double ** allocate_double_matrix(int row, int col);
int ** allocate_int_matrix(int row, int col);
int * allocate_int_vector(int n);
void DestroyDoubleMatrix(double ** Mat,int nrow);
void DestroyIntMatrix(int ** Mat, int nrow);
void DestroyDoubleVector(double * vect);
void DestroyIntVector(int * vet);

//random number generation functions
double genrand(); 
void sgenrand(unsigned long);/* initializing random generator with a NONZERO seed */

//cluster functions
int find_max_el(int *cluster, int Nsphere);
void findClusterRec(int * cluster,  int Nsphere, int** AdjList, int col);
void growClusterRec(int * cluster,  int Nsphere, int** AdjList, int col, int x, int count);
void adjust_group(int * cluster, int Nsphere);
int Ngroup_search(int *vet, int n);
int bool_search(int *vet, int n, int num);
void ClusterRec(int *clusterRec, int Nsphere, int **AdjList, std::ofstream& nameFile);

void ClusterHK(int *clusterHK, int Nsphere, int **AdjList, std::ofstream& nameFile);
int count_zeros(int *vet, int x, int m, int **AdjList);

int percolation (int Nsphere, int ** connect, int *cluster, double *count_perc, int index, int realiz);
void percolation_wrap (int Nsphere, int ** connect, int *clusterHK, int *count_perc, int realiz, double** posCond, double Lato, double diam, double shell);
 
//print functions
void PrintAdjList(int Nsphere, int **AdjList);
void print_filePos_impenetrable_beforeMC(long int Nsphere,  double **pos);
void print_fileDir_impenetrable_beforeMC(long int Nsphere,  double **pos);
void print_filePos_impenetrable_afterMC(long int Nsphere,  double **pos);
void print_fileDir_impenetrable_afterMC(long int Nsphere,  double **pos);
void print_int_vect(int * vect, int n , std::ofstream& nameFile);


//observable functions
void update_hist_gr(double ** pos, int N, double *hist, double Lato, double deltar, int Nbin);
void normal_write_gr(int Nsphere, double *hist, double Lato, double deltar, int Nbin, int realiz);
void update_hist_Hp(double **pos, int ** AdjList, int Nsphere, double Lato, double *hist, double deltar, int Nbin);
void normal_write_Hp(int Nsphere, double *hist, double deltar, int Nbin, int realiz);
void update_hist_gr_angle(double ** direct, double **pos, int Nsphere, double *hist, double Lato, double deltar, int Nbin);
void normal_write_gr_angle(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz);
void normal_write_cc(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz);
void update_hist_gr_mindist(double ** direct, double **pos, int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin);
void normal_write_gr_mindist(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz, string suffix);
double check_order(double ** direct,  int Nsphere, double length);

//Conductivity functions
void sort(int* array, int begin, int end, int* index);
void sortDeque(deque <int> &array, int begin, int end, deque <double> &index);
int ricercaSequenziale(int *array, int x, int n);
int ricercaSequenzialeDeque(deque <int> array, int x, int n);
void swapDouble(double * a, double * b, double tt);

void createCondStruct(double diamCond, int Nsphere, double Lato, double **pos, int *cluster, int **AdjList, int **connect, deque <deque <int> > &Adj, deque < deque <double> > &CondList, int * NumNeigh, int span_index, int *countPart);

void createCondStructCyl(double diamCond, double length, int Nsphere, double Lato, double **pos, double **direct, int *cluster, int **AdjList, int **connect, deque <deque <int> > &Adj, deque < deque <double> > &CondList, int * NumNeigh, int span_index, int *countPart);

double decimation(double diamCond, int Nsphere, double Lato, deque <deque <int> > &Adj, deque < deque <double> > &CondList,  int * NumNeigh, int *NunNeighIndex, int span_index, int countPart, int max_decimated, double min_cond);

void calcConductivity(double diamCond, int Nsphere, double Lato, double **pos, int *cluster, int **AdjList, int **connect, int span_index, double *conductance, double toll, double min_cond);

void calcConductivityCyl(double diamCond, double length, int Nsphere, double Lato, double **pos, double **direct, int *cluster, int **AdjList, int **connect, int span_index, double *conductance, double toll, double min_cond);

double createCGMstruct(double diamCond, int Nsphere, double Lato, deque <deque <int> > &Adj, deque < deque <double> > &CondList, double toll);

void kersh(int nequ,int nterm, vector <int> ia , vector <int> ja,  vector <double> sysmat, double *prec);
void lsolve(int nequ,int nterm, vector <int> ia , vector <int> ja, double *lmat, double *vec, double *pvec);
void lsolve2(int n,int nt, vector <int> ia , vector <int> ja, double *prec, double *v, double *w);
void MatVecMult(int n,int nterm, vector <int> ia , vector <int> ja,  vector <double> sysmat, double *invec, double *outvec);
#endif

