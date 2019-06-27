#include <stdio.h>
#include <omp.h>
#include <time.h>

extern void compro_1_(int, char*, int, double*, double*, double*, double*, 
int, int, double, double, double*, double*, double*);

#define NUMPHA_TD 	21
#define MAXENDMEM_TD 	7
#define MAX_NAMLEN 	10
#define MAX_NOX 	15
#define MAX_NPHASES 	10

#define MANTLE		1

struct material
{
  //char name[128];
  char name; // 0=water/air, 1=crust, 2=mantle
  double Q, k;
  int Nphases; // 0 if fixed
  int Noxides; // 0 if fixed
  char namphatxt[MAX_NPHASES * MAX_NAMLEN];
  double phacomtxt[MAX_NPHASES * (MAX_NOX+4)];
  double endmemrat[MAX_NPHASES * MAXENDMEM_TD];
  double molref[MAX_NPHASES], volref[MAX_NPHASES];
  // NOTE:
  // volref/molref vs phacomtxt[:][1:2]
  // wt% (phacomtxt[:][0]) is unused
  // mol (phacomtxt[:][3]) is recomputed by compro_1
  // oxides are unused
};

// NOTE: h[m], vel[km/s], den[g/cm^3], T[K], Q[W/m^3], k[W/(m*K)], serp[weight %]
int getvel(struct material* layers, int* layerID, double* T, double* vel, double* den,
	   double h, int nex, int ney, int nez, double* serp,
	   // for heat solver
	   int updateT, double Ttop, double Tbot, char* mask, double tol, double omega)
{
  int i=0, j=0, k=0;
  int Neyz = ney*nez;

  double dP = 2670.*9.8*1.e-5*h; // [bar]

  // serpentinized mantle of Carlson & Miller, 2003
  double vm = 8.2; // vp of peridotite mantle 
  double a = vm*0.13/0.31 * 0.01; // 0.13 - weight fraction of water in serpentinite 
  double dm = 3.3; // density of peridotite mantle 
  double b = (dm - 2.485)/2.485 * 0.01; // 2.485 - density of serpentinite
  
// built-in heat solver
  if(updateT)
      getT_SOR_nl(layers, layerID, mask, T, nex, ney, nez, h, Ttop, Tbot, tol, omega);
//
  
#pragma omp parallel for private(i,j,k)
  for(k = 0; k < nez; ++k)
  {
    double P = 1.+dP*k;
    for(i = 0; i < nex; ++i)
      for(j = 0; j < ney; ++j)
      {
	int ijk = i*Neyz + j*nez + k;
	struct material* m = &(layers[layerID[ijk]]);

	if(m->Nphases > 0)
	{
	  double vssys;
	  
	  compro_1_(m->Nphases, m->namphatxt, m->Noxides, m->molref, m->volref, m->phacomtxt, 
	  m->endmemrat, NUMPHA_TD, MAXENDMEM_TD, P, T[ijk], &vel[ijk], &vssys, &den[ijk]);

	  den[ijk] *= 1.e-3;
	}
	else if(m->name == MANTLE)
	{
	  // serpentinized mantle of Carlson & Miller, 2003
	  vel[ijk] = vm - a * serp[ijk];
	  den[ijk] = dm / (1. + b * serp[ijk]);
	}
      }
  }
  
  return 0;
}


int getT_SOR_nl(struct material* layers, int* layerID,char* mask, double* T, 
		     int nex, int ney, int nez, double h, 
		     double Ttop, double Tbot, double tol, double omega)
{
    int i=0, j=0, k=0, iter=0;
    int Neyz = ney*nez;
    int Nele = Neyz*nex;
    
    // set Dirichlet BC at top and bottom, linearly interpolate in between
    int nez1 = nez-1;
    for(k = 0; k < nez; ++k)
    {
	double T0 = (k*Tbot + (nez1-k)*Ttop) / nez1;
	for(i = 0; i < nex; ++i)
	    for(j = 0; j < ney; ++j)
		T[i*Neyz + j*nez + k] = T0;
    }
    // set Ttop at masked area
    for(i = 0; i < Nele; ++i)
	if(mask[i])
	    T[i] = Ttop;
	
    double* W = (double*)malloc(Nele * sizeof(double));
    double* b = (double*)malloc(Nele * sizeof(double));
    double h2 = h*h; 
    
    int iterLim = 1;
    for(iter = 0; iter < iterLim; ++iter)
    {
// 	printf("\niteration: %d\n", iter);
	for(i = 0; i < Nele; ++i)
	{
	    struct material* m = &(layers[layerID[i]]);
	    W[i] = m->k;
	    b[i] = m->Q * h2;
	    
	    // k(T) of Levy, Jaupart & Mareschal, 2010
	    if(m->name == MANTLE) // mantle
		W[i] = W[i] * sqrt(298./T[i]) + 0.368e-9*T[i]*T[i]*T[i];
	    else // crust and mask (for continuity)
		W[i] = 2.26 + (355.576*W[i] - 618.241)/T[i] - 0.30247*W[i];
	    // Ben
//		W[i] = W[i] * pow(298./T[i], 0.33) + 0.368e-9*T[i]*T[i]*T[i];
	}  
	
	if(getT_SOR(W, b, mask, T, nex, ney, nez, tol, omega))
	  return 1;
    }
    
    free((char*)W);
    free((char*)b);    

    return 0;
}

// TODO: red-black ordering for parallelization
int getT_SOR(double* W, double* b, char* mask, double* x, 
		     int nex, int ney, int nez, double tol, double omega)
{
    int i=0, j=0, k=0;
    int Neyz = ney*nez;
    int Nele = Neyz*nex;
    
    // double omega = 1.97; 
    double omega1 = 1. - omega;
    int max_iter = Nele;
    double rchange = 0.;

//     printf("\nSOR:\nomega = %g, max_iter = %d, tol = %g\n", omega, max_iter, tol);
    time_t t0, t1;
    time(&t0);
    
//     double x0 = (Ttop+Tbot)/2.;
//     for(i = 0; i < Nele; ++i)
//       x[i] = x0;
//     // set Dirichlet BC at top and bottom
//     for(i = 0; i < nex; ++i)
// 	for(j = 0; j < ney; ++j)
// 	{
// 	    x[i*Neyz + j*nez] = Ttop;
// 	    x[i*Neyz + j*nez + nez-1] = Tbot;
// 	}
//     // set Ttop at masked area
//     for(i = 0; i < Nele; ++i)
// 	if(mask[i])
// 	    x[i] = Ttop;
    
    // SOR
    int iter;
    for(iter = 0; iter < max_iter; ++iter)
    {
        double x2 = 0., dx2 = 0.;
	
#pragma omp parallel for private(i, j, k) /*reduction(+ : x2)*/ reduction(max : dx2)
	for(j = 0; j < ney; ++j)
		for(i = 0; i < nex; ++i)
			for(k = 1; k < nez-1; ++k)
			{
				int ijk = i*Neyz + j*nez + k;
				
				if(mask[ijk])
					continue;
				
				double rhs = b[ijk];
				double lhs = 0.;
				
				// x
				if(i != 0)
				{
					rhs += W[ijk - Neyz] * x[ijk - Neyz];
					lhs += W[ijk - Neyz];
				}
				if(i != nex-1)
				{
					rhs += W[ijk] * x[ijk + Neyz];
					lhs += W[ijk];
				}

				// y
				if(j != 0)
				{
					rhs += W[ijk - nez] * x[ijk - nez];
					lhs += W[ijk - nez];
				}
				if(j != ney-1)
				{
					rhs += W[ijk] * x[ijk + nez];
					lhs += W[ijk];
				}

				// z
				rhs += W[ijk - 1] * x[ijk - 1] + W[ijk] * x[ijk + 1];
				lhs += W[ijk - 1] + W[ijk];

				double x_old = x[ijk];
				x[ijk] = omega * rhs / lhs + omega1 * x[ijk];

				// calculate relative change wrt the previous iteration
				// use L_infinity norm or L2 norm
				//x2 += x_old * x_old;
				double b = x_old * x_old;
				x_old -= x[ijk];
				//dx2 += x_old * x_old;

				double a = x_old * x_old / b;
				if(a > dx2)
				  dx2 = a;
			}

	rchange = sqrt(dx2);//sqrt(dx2 / x2);
        if(rchange < tol)
        {
	  time(&t1);
//	  printf("\nSOR solver converged after %d iterations\nElapsed time: %g s\n", iter,
//		 difftime(t1, t0));
          return 0;
        }
    }

    printf("\nSOR solver failed! Iteration limit is reached with rchange %g\n", rchange);

    return 1;
}
