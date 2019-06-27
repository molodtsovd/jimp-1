#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <vector>
#include <math.h>
#include "dumb_logger.h"
#include "point3d.h"


using namespace std;

#define PI2		1.5707963267948966
#define BG_DENSITY	2.67 // Bouguer density
#define G 		6.674 // gravitational constant: g[mGal] = G * density[g/cm3] * distance[km]

class grav
{
	// logger  
	dumb::logger jiout;
	
	// model
	float h;
	int nx, ny, nz, Nnode, Nyz, Nele, Neyz, nez;
	int Nele_ext, Neyz_ext, nez_ext, nez_air; // vertically extended mesh
	point3d R0, R1;
	
	double* m;//anomalous density [Nele]
	double* u;//potential (for Poisson solver) [Nele_ext]
	char* mask;//[Nele]
	float* A;//[Nd * Nele]

	// data
	int Nd;
	vector<double> g_cal; //[Nd]
	vector<double> g_mask; //[Nd]
	vector<double> g_obs; //[Nd]
	vector<double> Wd; //[Nd]
	vector<point3d> rcv; //[Nd]
	double g_obs_av;
	// for isostasy
	double* elev; // 2D: ny-1
	double* elev_cal; // 2D: ny-1
	double* elev_obs; // TMP for synthetic 'observed' topo, 2D: ny-1
	double rho_a_iso, H0_iso, Wd_iso;

	// values of functionals
	double Phi_d, Phi_iso;


public:
	// constructor: define log file path
	grav(point3d R0_, point3d R1_, double h_, int nez_air, char* d_filename, char* bath_filename, 
	     char* iso_filename, char* elev_obs_filename);
	
	double FP_SOR(double* den, double tol, double omega);
	int SOR(double* den, double tol, double omega);
	
	void getA_2D(); // char* mask
	double FP_2D(double* den, double& RMSg, double& RMSe);
	
	void FP_mask(double* den, char* mask_);
	double FP(double* den, bool fixed=true);
	
	// TEMP
	void saveGD(char* fname);
	void saveElev(char* fname);
};	

inline double quadPlane(double dx, double dy, double dz, double sgn);


extern "C"
{
    grav* grav_(point3d R0_, point3d R1_, double h_, int nez_air_, char* d_filename, char* bath_filename,
		char* iso_filename, char* elev_obs_filename)
    {
      return new grav(R0_, R1_, h_, nez_air_, d_filename, bath_filename, iso_filename, elev_obs_filename);
    }
//     void FP_mask_(grav* a, double* m, char* mask_)
//     {
//       return a->FP_mask(m, mask_);
//     }
    double FP_SOR_(grav* a, double* m, double tol, double omega)
    {
      return a->FP_SOR(m, tol, omega);
    }
    double FP_2D_(grav* a, double* m, double& RMSg, double& RMSe)
    {
      return a->FP_2D(m, RMSg, RMSe);
    }
    void saveGD_(grav* a, char* fname)
    {
      return a->saveGD(fname);
    }
    void saveElev_(grav* a, char* fname)
    {
      return a->saveElev(fname);
    }
}

