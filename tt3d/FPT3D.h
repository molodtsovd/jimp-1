#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <vector>
#include "dumb_logger.h"
#include "point3d.h"
extern "C" 
{
#include "fdtimes.h"
}


using namespace std;
#define TIME_3D_MSG 	0
#define	TT3D_INFINITY  	0.500e+19
#ifndef EPS
#define EPS		1.e-9
#endif
#ifndef M_SQRT3
#define M_SQRT3     	1.732050807568877076
#endif
#ifndef PI
#define PI		3.141592653589793
#endif


class TT3D
{
	// logger  
	dumb::logger jiout;
	
	// model
	float h;
	int nx, ny, nz, Nnode, Nyz, Nele, Neyz, nez;
	point3d R0, R1;
	
	double* m;//[Nele]
	//char* fixed;//[Nele]

	// memory used by FP
	float* hs;//[Nnode];

	// data
	int Nt, Nsrc, Nrcv;
	vector<double> Tobs; //[Nt]
	vector<double> Tcal; //[Nt]
	vector<double> Wd; //[Nt]
	vector<int> srcID; //[Nt]
	vector<int> rcvID; //[Nt]
	vector<point3d> src; //[Nsrc]
	vector<point3d> rcv; //[Nrcv]

	// values of functionals
	double Phi_d;


public:
	// constructor: define log file path
	TT3D(point3d R0_, point3d R1_, double h_, char* ft_filename, char* bath_filename);
	void setVref(double* ref);
	double FP(double* vel, double& RMS);
	void getTcal(float* t, vector<point3d>& rcvID, vector<double>& T);
	int update_hs(double* vel);
	
	// TEMP
	void saveTcal(char* fname, double err);
	void saveSD(char* fname);
	void saveT(char* fname);
	void sliceNnode(char* fname, float* v, bool vel);
	void saveNele(char* fname, double* v, bool vel);
	void sliceNele(char* fname, bool vel);

};	
double GaussPseudoRand(double mean, double sd);

extern "C"
{
    TT3D* TT3D_(point3d R0_, point3d R1_, double h_, char* ft_filename, char* bath_filename)
    {
      return new TT3D(R0_, R1_, h_, ft_filename, bath_filename);
    }
    int FP_(TT3D* a, double* m, double& RMS)
    {
      return a->FP(m, RMS);
    }
    void saveSD_(TT3D* a, char* fname)
    {
      return a->saveSD(fname);
    }
    void saveTcal_(TT3D* a, char* fname, double err)
    {
      return a->saveTcal(fname, err);
    }
    void setVref_(TT3D* a, double* v)
    {
      return a->setVref(v);
    }
}

