#include "grav.h"

grav::grav(point3d R0_, point3d R1_, double h_, int nez_air_, char* d_filename, char* bath_filename, 
	   char* iso_filename, char* elev_filename)
{
	jiout.init("grav-log.txt");
	
	// model dimensions
	h = h_;
	R0 = R0_;
	R1 = R1_;
	point3d R = R1-R0;
	
	nx = ceil(R.x / h) + 1;
	ny = ceil(R.y / h) + 1;
	nz = ceil(-R.z / h) + 1;
	
	printf("\nnx=%d, ny=%d, nz=%d, nez_air=%d\n",nx,ny,nz,nez_air_);

	// 
	Nnode = nx*ny*nz;
	Nyz = ny*nz;
	Nele = (nx-1)*(ny-1)*(nz-1);
	Neyz = (ny-1)*(nz-1);
	nez = nz-1;

	// allocate memory for inversion arrays
	m = new double[Nele];
	mask = new char[Nele];
	for(int i = 0; i < Nele; ++i)
	  mask[i] = 0;

// vertically extended mesh for potential
	nez_air = nez_air_;
	nez_ext = nez + nez_air;
	Neyz_ext = (ny-1)*nez_ext;
	Nele_ext = (nx-1)*Neyz_ext;
	u = new double[Nele_ext];

	
// load geometry & observed data
	FILE* f = fopen(d_filename, "r");
	char buf[1024];

	// X,Y,Z,GOBS,ERR,GCAL
	double X,Y,Z,GOBS,ERR,GCAL;

	while(fscanf(f, "%s", buf) != EOF)
	{
		if(buf[0] == '/') 
			continue; // ignore comments
		
		if(!sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf", &X, &Y, &Z, &GOBS, &ERR, &GCAL))	
		{
			jiout << sprintf(jiout.buf, "Error while parsing ASCII string \"%s\"", buf);
		}
		point3d p = {X - R0.x, Y - R0.y, R0.z - Z};
		rcv.push_back(p);

		g_obs.push_back(GOBS);
		g_cal.push_back(0.);
		g_mask.push_back(0.);
		double var_d = ERR*ERR;
		Wd.push_back(1./(var_d > 0. ? var_d : 1.));
	}
	Nd = g_obs.size();

	fclose(f);


// constant component of observed data
	g_obs_av = 0.;
	for(int n = 0; n < Nd; ++n)
	  g_obs_av += g_obs[n];
	g_obs_av /= Nd;
	
	
// compute 2D FP matrix A
	getA_2D();
	
// // define starting density model
// 
// NOTE: bathymetry is currently 2D

	// load bathymetry
	f = fopen(bath_filename, "r");

	vector<double> bathx_buf;
	vector<double> bathz_buf;
	elev = new double[ny-1];
	elev_cal = new double[ny-1];
	double x, z;
	while(fscanf(f, "%s", buf) != EOF)
	{
		if(!sscanf(buf, "%lf,%lf", &x, &z))	
		{
			jiout << sprintf(jiout.buf, "Error while parsing ASCII string \"%s\"", buf);
		}
		bathx_buf.push_back(x);
		bathz_buf.push_back(z);
	}
	fclose(f);
	

// load synthetic 'observed' elevation
	f = fopen(elev_filename, "r");

	elev_obs = new double[ny-1];
	int i = 0;
	while(fscanf(f, "%s", buf) != EOF)
	{
		if(buf[0] == '/') 
		  continue;
		if(!sscanf(buf, "%*lf,%*lf,%lf", &elev_obs[i]))	
		{
			jiout << sprintf(jiout.buf, "Error while parsing ASCII string \"%s\"", buf);
		}
		++i;
	}
	fclose(f);
///

	// translate bathymetry to the model grid
	// X in bathymetry file must be ascending
	int j = 0;
	for(int i = 0; i < ny-1; ++i)
	{
		double x = (double(i) + 0.5) * h + R0.y;
		while(bathx_buf[j] < x)
			++j;
		elev[i] = bathz_buf[j-1] + (bathz_buf[j] - bathz_buf[j-1]) / (bathx_buf[j] - bathx_buf[j-1]) * (x - bathx_buf[j-1]);
	}


	// load isostasy parameters
	f = fopen(iso_filename, "r");
	fscanf(f, "Wd=%lf\nH0=%lf\nrho_a=%lf", &Wd_iso, &H0_iso, &rho_a_iso);	
	//printf("Wd=%g\nH0=%g\nrho_a=%g", Wd_iso, H0_iso, rho_a_iso);	
	fclose(f);
	
// 
// 	// build the model
// 	double dw = 1. - BG_DENSITY;
// 	for(int i = 0; i < nx-1; ++i)
// 		for(int j = 0; j < ny-1; ++j)
// 			for(int k = 0; k < nz-1; ++k)
// 			{
// 				if(R0.z - (double(k)+0.5)*h > elev[j])
// 					m[i*Neyz + j*nez + k] = dw;
// 				else
// 					m[i*Neyz + j*nez + k] = 0.;
// 			}
}

//////// 2D

double g_rect(double dx, double dz)
{
	if(dx == 0. && dz == 0.)
	  return 0.;
	
	return dx * log(dx*dx + dz*dz) + 2. * dz * atan(dx / dz); 
}
double g_halflayer(double dx, double dz1, double dz2)
{
	if(dx == 0. && (dz1 == 0. || dz2 == 0.))
	  return 0.;
	
	// NOTE: z axis points upward (cf. standard formula)
	double dx2 = dx*dx;
	return -dx * log((dx2 + dz2*dz2) / (dx2 + dz1*dz1)) - 2. * dz2 * (PI2 + atan(dx / dz2)) + 2. * dz1 * (PI2 + atan(dx / dz1)); 
}
// consider YZ plane
void grav::getA_2D() // char* mask
{
	A = new float[Nd * Nele];
#pragma omp parallel for
	for(int n = 0; n < Nd; ++n)
	{
		double y = rcv[n].y;
		double z = rcv[n].z;
		for(int k = 0; k < nez; ++k)
		{
			double y2 = h;
			double y1 = (ny-2)*h;
			double z1 = -k*h;
			double z2 = -(k+1)*h;
			for(int i = 0; i < nx-1; ++i)
			{
				A[n*Nele + i*Neyz + k] = G * g_halflayer(y-y2, z1-z, z2-z);
				A[n*Nele + i*Neyz + (ny-2)*nez + k] = G * g_halflayer(y1-y, z1-z, z2-z);
			}
		}
		for(int j = 1; j < ny-2; ++j)
			for(int k = 0; k < nez; ++k)
			{
				double y1 = j*h;
				double y2 = (j+1)*h;
				double z1 = -k*h;
				double z2 = -(k+1)*h;
				for(int i = 0; i < nx-1; ++i)
					A[n*Nele + i*Neyz + j*nez + k] = \
					G * (g_rect(y1-y, z1-z) + g_rect(y2-y, z2-z) - \
					g_rect(y1-y, z2-z) - g_rect(y2-y, z1-z));
			}
	}
}

double grav::FP_2D(double* den, double& RMSg, double& RMSe)
{
  double Phi_d = 0.;
// #pragma omp parallel for schedule(dynamic) reduction(+ : Phi_d)
//   for(int n = 0; n < Nd; ++n)
//   {
//     g_cal[n] = 0.;
//     for(int i = 0; i < Nele; ++i)
//       g_cal[n] += A[n*Nele + i] * den[i];
// 
//     double r = g_cal[n] - g_obs[n];
//     Phi_d += r*Wd[n]*r;
//   }

  // substract constant component from g_obs and g_cal
#pragma omp parallel for schedule(dynamic)
  for(int n = 0; n < Nd; ++n)
  {
    g_cal[n] = 0.;
    for(int i = 0; i < Nele; ++i)
      g_cal[n] += A[n*Nele + i] * den[i];
  }
  double g_cal_av = 0.;
  for(int n = 0; n < Nd; ++n)
    g_cal_av += g_cal[n];
  g_cal_av /= Nd;
  for(int n = 0; n < Nd; ++n)
  {
    double r = g_cal[n] - g_cal_av - g_obs[n] + g_obs_av;
    Phi_d += r*Wd[n]*r;
  }
  RMSg = sqrt(Phi_d / Nd);

  
  // isostasy (Lachenbruch & Morgan, 1990)
  // NOTE: tested only in marine case; 2D
  double Phi_d_iso = 0.;
  int Nd_iso = ny-1;
  
  double rho_w = 1.;
  
// //// for local calibration:
//   int iref = ny/2; // index of reference column
//   double Lref = nez*h + elev[iref];
//   int j = round(-elev[iref] / h); // number of water cells
//   double rho_ref = (elev[iref]/h + j) * den[iref*nez + j]; // fraction of the cell below bathymetry
// 
//   for(; j < nez; ++j)
//       rho_ref += den[iref*nez + j];
//   rho_ref *= h/Lref;
//   
//   H0_iso = ((rho_a_iso - rho_ref) * Lref - (rho_a_iso - rho_w) * elev[iref]) / rho_a_iso;     
// ////
  
  double a = h / rho_a_iso;
  double w = rho_a_iso / (rho_a_iso - rho_w);
  for(int i = 0; i < Nd_iso; ++i)
  {
    double L = nez*h + elev[i];
    int j = round(-elev[i] / h); // number of water cells
    double rho = (elev[i]/h + j) * den[i*nez + j]; // fraction of the cell below bathymetry
    for(; j < nez; ++j)
	rho += den[i*nez + j];
    
    elev_cal[i] = L - rho*a - H0_iso;
    if(elev[i] < 0.)
      elev_cal[i] *= w;

    double r = elev_cal[i] - elev_obs[i];
    Phi_d_iso += r*Wd_iso*r;
  }
  RMSe = sqrt(Phi_d_iso / Nd_iso);

  return 0.5 * (Phi_d + Phi_d_iso);
}

////////


double grav::FP_SOR(double* b, double tol, double omega)
{  
    if(SOR(b, tol, omega))
      return NAN;
    
    double Phi_d = 0.;
    for(int n = 0; n < Nd; ++n)
    {
      int i = static_cast<int>(rcv[n].x / h);
      int j = static_cast<int>(rcv[n].y / h);
      int k = static_cast<int>(rcv[n].z / h) + nez_air/2;
      int ijk = i*Neyz_ext + j*nez_ext + k;
      
      g_cal[n] = (u[ijk - 1] - u[ijk + 1]) / (2.*h);
      
      double r = g_cal[n] - g_obs[n];
      Phi_d += r*Wd[n]*r;
    }

  // normalize by number of data
  Phi_d /= Nd;
  
  return Phi_d;
}


// TO DO: red-black ordering for parallelization
int grav::SOR(double* b, double tol, double omega)
{  
    // double omega = 1.97; 
    double omega1 = 1. - omega;
    int max_iter = Nele_ext;
    double rchange = 0.;

    printf("\nSOR:\nomega = %g, max_iter = %d, tol = %g\n", omega, max_iter, tol);
    time_t t0, t1;
    time(&t0);
    
    for(int i = 0; i < Nele_ext; ++i)
      u[i] = 0.;
    // set Dirichlet BC at top and bottom
    double u0 = 0.;//h*nez_air * 2. * 3.141592653589793 * G * (3.-2.67) * 40;
    for(int i = 0; i < nx-1; ++i)
	for(int j = 0; j < ny-1; ++j)
	{
	    u[i*Neyz_ext + j*nez_ext] = u0;
	    u[i*Neyz_ext + j*nez_ext + nez_ext-1] = 0.;
	}
    
    // SOR
    for(int iter = 0; iter < max_iter; ++iter)
    {
        double x2 = 0., dx2 = 0.;

#pragma omp parallel for /*reduction(+ : x2)*/ reduction(max : dx2)
	for(int j = 0; j < ny-1; ++j)
		for(int i = 0; i < nx-1; ++i)
			for(int k = nez_ext-2; k > 0; --k)
			{
				int ijk = i*Neyz_ext + j*nez_ext + k;
				
				double rhs = b[ijk];
				double lhs = 0.;

				// x
				if(i != 0)
				{
					rhs += u[ijk - Neyz_ext];
					lhs += 1.;
				}
				if(i != nx-2)
				{
					rhs += u[ijk + Neyz_ext];
					lhs += 1.;
				}

				// y
				if(j != 0)
				{
					rhs += u[ijk - nez_ext];
					lhs += 1.;
				}
				if(j != ny-2)
				{
					rhs += u[ijk + nez_ext];
					lhs += 1.;
				}

				// z
				rhs += u[ijk - 1];
				lhs += 1.;
				if(k != nez_ext-1)
				{
					rhs += u[ijk + 1];
					lhs += 1.;
				}

				double x_old = u[ijk];
				u[ijk] = omega * rhs / lhs + omega1 * u[ijk];

				// calculate relative change wrt the previous iteration
				// use L_infinity norm or L2 norm
				//x2 += x_old * x_old;
				double b = x_old * x_old;
				x_old -= u[ijk];
				//dx2 += x_old * x_old;

				double a = x_old * x_old / b;
				if(a > dx2)
				  dx2 = a;
			}

	rchange = sqrt(dx2);//sqrt(dx2 / x2);
        if(rchange < tol)
        {
	  time(&t1);
	  printf("\nSOR solver converged after %d iterations\nElapsed time: %g s\n", iter,
		 difftime(t1, t0));
          return 0;
        }
    }

    printf("\nSOR solver failed! Iteration limit is reached with rchange %g\n", rchange);

    return 1;
}



void grav::FP_mask(double* den, char* mask_)
{
  for(int i = 0; i < Nele; ++i)
    mask[i] = mask_[i];
  for(int i = 0; i < Nd; ++i)
    g_mask[i] = 0.;

  FP(den, false);

  for(int i = 0; i < Nd; ++i)
    g_mask[i] = g_cal[i];
}


double grav::FP(double* den, bool fixed)
{
  double Gh3 = G*h*h*h;
  double Gh2 = G*h*h;
  double Gh = G*h;
  
#pragma omp parallel for
  for(int n = 0; n < Nd; ++n)
  {
    g_cal[n] = g_mask[n];
    for(int i = 1; i < nx-2; ++i)
	    for(int j = 1; j < ny-2; ++j)
		    for(int k = 0; k < nz-1; ++k)
		    {
		      if(mask[i*Neyz + j*nez + k] == fixed)
			continue;
		      
		      point3d R = {((double)i+0.5)*h, ((double)j+0.5)*h, -((double)k+0.5)*h};
		      R = rcv[n] - R;
		      double d = r(R);
		      g_cal[n] += Gh3 * R.z / (d*d*d) * den[i*Neyz + j*nez + k];
		    }
    // i = 0, nx-2
    for(int j = 1; j < ny-2; ++j)
	    for(int k = 0; k < nz-1; ++k)
	    {
	      int jk = j*nez + k;
	      if(mask[jk] == fixed && mask[(nx-2)*Neyz + jk] == fixed)
		continue;
	      
	      point3d R = {0.5*h, ((double)j+0.5)*h, -((double)k+0.5)*h};
	      R = rcv[n] - R;
	      double a = Gh2 * R.z / (R.y*R.y + R.z*R.z);
	      if(mask[jk] != fixed)
		g_cal[n] += a * (1. + R.x / r(R)) * den[jk];
	      
	      R.x -= (nx-2)*h;
	      if(mask[(nx-2)*Neyz + jk] != fixed)
		g_cal[n] += a * (1. - R.x / r(R)) * den[(nx-2)*Neyz + jk];
	    }

    // j = 0, ny-2
    for(int i = 1; i < nx-2; ++i)
	    for(int k = 0; k < nz-1; ++k)
	    {
	      int ik = i*Neyz + k;
	      if(mask[ik] == fixed && mask[ik + (ny-2)*nez])
		continue;

	      point3d R = {((double)i+0.5)*h, 0.5*h, -((double)k+0.5)*h};
	      R = rcv[n] - R;
	      double a = Gh2 * R.z / (R.x*R.x + R.z*R.z);
	      if(mask[ik] != fixed)
		g_cal[n] += a * (1. + R.y / r(R)) * den[ik];

	      R.y -= (ny-2)*h;
	      if(mask[ik + (ny-2)*nez] != fixed)
		g_cal[n] += a * (1. - R.y / r(R)) * den[ik + (ny-2)*nez];
	    }

    // corners
    for(int k = 0; k < nz-1; ++k)
    {
      // i = 0, j = 0
      int ijk = k;
      if(mask[ijk] != fixed)
      {  
	point3d R = {0.5*h, 0.5*h, -((double)k+0.5)*h};
	R = rcv[n] - R;
	g_cal[n] += Gh * quadPlane(R.x, R.y, R.z, -1.) * den[ijk];
      }
      // i = nx-2, j = 0
      ijk = (nx-2)*Neyz + k;
      if(mask[ijk] != fixed)
      {  
	point3d R = {((double)(nx-2) + 0.5)*h, 0.5*h, -((double)k+0.5)*h};
	R = rcv[n] - R;
	g_cal[n] += Gh * quadPlane(-R.x, R.y, R.z, -1.) * den[ijk];
      }
      // i = nx-2, j = ny-2
      ijk = (nx-2)*Neyz + (ny-2)*nez + k;
      if(mask[ijk] != fixed)
      {  
	point3d R = {((double)(nx-2) + 0.5)*h, ((double)(ny-2) + 0.5)*h, -((double)k+0.5)*h};
	R = rcv[n] - R;
	g_cal[n] += Gh * quadPlane(R.x, R.y, R.z, 1.) * den[ijk];
      }
      // i = 0, j = ny-2
      ijk = (ny-2)*nez + k;
      if(mask[ijk] != fixed)
      {  
	point3d R = {0.5*h, ((double)(ny-2) + 0.5)*h, -((double)k+0.5)*h};
	R = rcv[n] - R;
	g_cal[n] += Gh * quadPlane(-R.x, R.y, R.z, 1.) * den[ijk];
      }
    }
  }
  
  Phi_d = 0.;
  for(int n = 0; n < Nd; ++n)
  {
    double r = g_cal[n] - g_obs[n];
    Phi_d += r*Wd[n]*r;
  }
  return Phi_d;
}


inline double quadPlane(double dx, double dy, double dz, double sgn)
{
  
  return sgn*PI2 + atan(dy / dz) + sgn*atan(dy * dx / (dz * sqrt(dx*dx + dy*dy + dz*dz))) + atan(dx / dz);
}


double GaussPseudoRand(double mean, double sd)
{
	const double PI = 3.141592653589793;
	const double eps = 1.e-30;

	double u1 = double(rand())/RAND_MAX + eps;
	double u2 = double(rand())/RAND_MAX + eps;

	return sqrt(-2.*log(u1)) * cos(2.*PI*u2) * sd + mean;
}

void grav::saveGD(char* fname)
{
	FILE* f = fopen(fname, "w");
	fprintf(f, "//X,Y,Z,GOBS,ERR,GCAL\n");

	for(int k = 0; k < Nd; ++k)
	{
		float rx = R0.x + rcv[k].x;
		float ry = R0.y + rcv[k].y;
		float rz = R0.z - rcv[k].z;
		fprintf(f, "%f,%f,%f,%f,%f,%f\n", rx,ry,rz,g_obs[k],1./sqrt(Wd[k]),g_cal[k]);
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n Gravity data saved to file %s\n", fname);
}

void grav::saveElev(char* fname)
{
	FILE* f = fopen(fname, "w");
	fprintf(f, "//Y,ELEV_OBS,ELEV_CAL\n");

	float y = R0.y - 0.5*h;
	for(int k = 0; k < ny-1; ++k)
	{
		y += h;
		fprintf(f, "%f,%f,%f\n", y, elev_obs[k], elev_cal[k]);
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n Elevation data saved to file %s\n", fname);
}
