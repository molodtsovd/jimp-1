#include "FPT3D.h"

/*
int main()
{	
	TT3D tomo("test/rapids4_100.ttime", "test/rapids4_100.bath");

	double* m;
	tomo.FP(m);
	
	tomo.saveSD("test/fp.sd");

	return 0;
}
*/
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
double TT3D::FP(double* vel, double& RMS)
{
	int exitcode = 0;

	// update slowness model
	if(update_hs(vel))
		return NAN;

	// make threadprivate copies
	int Nth = 0;
	#pragma omp parallel reduction(+ : Nth)
	++Nth;
	//printf("\n\nNth=%d", Nth);
	// NOTE: hs_buf = hs is threadprivate in time_3d(), 
	// if hs is shared, when time_3d() is run in parallel it corrupts hs
	// anyway this copy will be required in MPI implementation
	float* _hs = new float[Nnode * Nth];
	for(int i = 0; i < Nth; ++i)
	    for(int j = 0; j < Nnode; ++j)
		_hs[i*Nnode + j] = hs[j];

	float* t = new float[Nnode * Nth];
	

	// parallel over sources
	#pragma omp parallel for
	for(int s = 0; s < Nsrc; ++s)
	{
		int thread = omp_get_thread_num();
		
		if(!exitcode)
		{
		//// make threadprivate copies of current source data
		vector<point3d> _rcv;
		vector<double> _Tcal;

		int i = 0;
		while(srcID[i] != s) ++i;
		while(i < Nt && srcID[i] == s)
		{
			_rcv.push_back(rcv[rcvID[i]]);
			_Tcal.push_back(0.);
			++i;
		}
		////

		// calculate eikonal
//		if(time_3d(_hs + thread*Nnode, t + thread*Nnode, nx, ny, nz, src[s].x, src[s].y, src[s].z, 0.001, TIME_3D_MSG))
		if(time_2d(_hs + thread*Nnode, t + thread*Nnode, ny, nz, src[s].y, src[s].z, 0.001, TIME_3D_MSG))
			exitcode = -1;
		
	// NOTE: for time_2d:
		for(int i = 0; i < Nyz; ++i)
		    t[thread*Nnode + Nyz + i] = t[thread*Nnode + i];
		
		// interpolate calculated eikonal at receivers
		getTcal(t + thread*Nnode, _rcv, _Tcal);

		// accumulate Tcal
		i -= _Tcal.size();
		for(int j = 0; j < _Tcal.size(); ++j)
		{
			Tcal[i] = _Tcal[j];
			++i;
		}
		
		}
	}
	delete[] t;
	delete[] _hs;

	// accumulate functional value
	Phi_d = 0.;
	for(int i = 0; i < Nt; ++i)
	{
		double dT = Tcal[i] - Tobs[i];
		Phi_d += dT * Wd[i] * dT;
	}
	RMS = sqrt(Phi_d / Nt);
	Phi_d *= 0.5;
	
	if(exitcode)
	{
		jiout << sprintf(jiout.buf, "Eikonal solver failed!");
		RMS = NAN;
		return NAN;
	}
	
	return Phi_d;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void TT3D::getTcal(float* t, vector<point3d>& _rcv, vector<double>& _Tcal)
{
	//jiout << sprintf(jiout.buf, "Interpolating traveltime field at receivers");

	int i0, i1, j0, j1, k0, k1;
	float t000, t100, t010, t001, t110, t101, t011, t111;
	float xd, yd, zd, _xd, _yd, _zd;

	// get time at receivers via trilinear interpolation
	for(int i = 0; i < _rcv.size(); ++i)
	{
		i0 = static_cast<int>(_rcv[i].x / h);
		j0 = static_cast<int>(_rcv[i].y / h);
		k0 = static_cast<int>(_rcv[i].z / h);
		if(i0 == nx-1) i0 -= 1;
		if(j0 == ny-1) j0 -= 1;
		if(k0 == nz-1) k0 -= 1;
		i1 = i0 + 1;
		j1 = j0 + 1;
		k1 = k0 + 1;


		if(i0 < 0 || j0 < 0 || k0 < 0 || i1 > nx-1 || j1 > ny-1 || k1 > nz-1)
		{
			jiout << sprintf(jiout.buf, "\n Receiver #%d (%g, %g, %g) out of mesh!\n", i, _rcv[i].x, _rcv[i].y, _rcv[i].z);
			continue;
		}

		xd = _rcv[i].x / h - i0;
		yd = _rcv[i].y / h - j0;
		zd = _rcv[i].z / h - k0;
		_xd = 1. - xd;
		_yd = 1. - yd;
		_zd = 1. - zd;

		t000 = t[i0*Nyz + j0*nz + k0];
		t100 = t[i1*Nyz + j0*nz + k0];
		t010 = t[i0*Nyz + j1*nz + k0];
		t001 = t[i0*Nyz + j0*nz + k1];
		t110 = t[i1*Nyz + j1*nz + k0];
		t101 = t[i1*Nyz + j0*nz + k1];
		t011 = t[i0*Nyz + j1*nz + k1];
		t111 = t[i1*Nyz + j1*nz + k1];

		_Tcal[i] = ((t000*_xd + t100*xd)*_yd + (t010*_xd + t110*xd)*yd)*_zd + ((t001*_xd + t101*xd)*_yd + (t011*_xd + t111*xd)*yd)*zd;
	}
}


int TT3D::update_hs(double* vel)
{
	for(int i = 0; i < nx-1; ++i)
		for(int j = 0; j < ny-1; ++j)
			for(int k = 0; k < nz-1; ++k)
			{
				double v = vel[i*Neyz + j*nez + k];
				if(v == -1.) // if fixed
				  v = m[i*Neyz + j*nez + k];
// 				else if(v < v_min || v > v_max)
// 				{
// 					jiout << sprintf(jiout.buf, "\n Invalid velocity at cell (%d, %d, %d) !\n", i, j, k);
// 					return 1;
// 				}
				
				hs[i*Nyz + j*nz + k] = h/v;
			}

	return 0;
}





TT3D::TT3D(point3d R0_, point3d R1_, double h_, char* ft_filename, char* bath_filename)
{
	jiout.init("tt3d-log.txt");
	
	// model dimensions
	h = h_;
	R0 = R0_;
	R1 = R1_;
	point3d R = R1-R0;
	
	nx = ceil(R.x / h) + 1;
	ny = ceil(R.y / h) + 1;
	nz = ceil(-R.z / h) + 1;
	
	printf("\nnx=%d, ny=%d, nz=%d\n",nx,ny,nz);

// 	double sigma_dT = 0.015; // seconds
// 	double var_dT = sigma_dT*sigma_dT;

	// 
	Nnode = nx*ny*nz;
	Nyz = ny*nz;
	Nele = (nx-1)*(ny-1)*(nz-1);
	Neyz = (ny-1)*(nz-1);
	nez = nz-1;

	// allocate memory for inversion arrays
	hs = new float[Nnode];
	m = new double[Nele];



// load geometry & Tobs
	FILE* ft = fopen(ft_filename, "r");
	char buf[1024];

	// NS,XS,ZS,NR,XR,ZR,T,ERR
	double NS,XS,ZS,NR,XR,ZR,T,ERR;

	int sID = -1, rID = -1;
	vector<int> srcIDshort; // [Nsrc] unique src IDs used to detect new srcs
	vector<int> rcvIDshort; // [Nrcv] unique rcv IDs used to detect new rcvs
	//int* rcvIDcheck = new int[Nrcv]; // TO DO: add number of rcv/src to header

	while(fscanf(ft, "%s", buf) != EOF)
	{
		if(buf[0] == '/') 
			continue; // ignore comments
		
		if(!sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &NS, &XS, &ZS, &NR, &XR, &ZR, &T, &ERR))	
		{
			jiout << sprintf(jiout.buf, "Error while parsing ASCII string \"%s\"", buf);
		}
		// add src
		if(!srcIDshort.size() || srcIDshort.back() != int(NS))
		{
			point3d p = {((double)nx - 1.)/2., (XS - R0.y)/h, (R0.z - ZS)/h};
			src.push_back(p);
			srcIDshort.push_back(int(NS));
			++sID;
		}
		// add rcv
		bool new_rcv = true;
		rID = 0;
		while(rID < rcvIDshort.size())
		{
			if(rcvIDshort[rID] == NR)
			{
				new_rcv = false;
				break;
			}
			++rID;
		}
		if(new_rcv)
		{
			point3d p = {h * ((double)nx - 1.)/2., XR - R0.y, R0.z - ZR};
			rcv.push_back(p);
			rcvIDshort.push_back(int(NR));
		}

		Tobs.push_back(T);
		Tcal.push_back(0.);
		Wd.push_back(1./(ERR != 0. ? ERR*ERR : 1.));
		srcID.push_back(sID);
		rcvID.push_back(rID);
	}
	Nsrc = src.size();
	Nrcv = rcv.size();
	Nt = Tobs.size();

	fclose(ft);

// define starting velocity model

// NOTE: bathymetry is currently 2D

	// load bathymetry
	FILE* txt = fopen(bath_filename, "r");

	vector<double> bathx_buf;
	vector<double> bathz_buf;
	double* bath = new double[ny-1];
	double x, z;
	while(fscanf(txt, "%s", buf) != EOF)
	{
		if(!sscanf(buf, "%lf,%lf", &x, &z))	
		{
			jiout << sprintf(jiout.buf, "Error while parsing ASCII string \"%s\"", buf);
		}
		bathx_buf.push_back(x);
		bathz_buf.push_back(z);
	}
	fclose(txt);

	// translate bathymetry to the model grid
	// X in bathymetry file must be ascending
	int j = 0;
	for(int i = 0; i < ny-1; ++i)
	{
		double x = (double(i) + 0.5) * h + R0.y;
		while(bathx_buf[j] < x)
			++j;
		bath[i] = bathz_buf[j-1] + (bathz_buf[j] - bathz_buf[j-1]) / (bathx_buf[j] - bathx_buf[j-1]) * (x - bathx_buf[j-1]);
	}

	// build the model
	double v0 = 1.5;
	double v_max = 9.5;
	double sw = 1/1.5;
	double dvdz = (v_max - v0)/sqrt(40.);//(v_max-v0)/40.*h;//

	for(int i = 0; i < nx-1; ++i)
		for(int j = 0; j < ny-1; ++j)
			for(int k = 0; k < nz-1; ++k)
			{
				if(R0.z - (double(k)+0.5)*h > bath[j])
					m[i*Neyz + j*nez + k] = sw;
				else
				{
					//double var = v_max - v0 - 0.2;
					m[i*Neyz + j*nez + k] = 1. / (dvdz * sqrt(k * h) + v0);//1./(v0 + double(k)*dvdz);//1./(double(rand())/RAND_MAX * var + v0 + 0.1);//
				}
			}
// load velocity
	//double vpvs = 1.7;
	//FILE* v_1 = fopen("v_1_.dat", "r");
	//for(int k = 0; k < nz-1; ++k)
	//{
	//	for(int j = 0; j < ny-1; ++j)
	//	{
	//		double v;
	//		fscanf(v_1, "%lf ", &v);
	//		if(R0.z - (double(k)+0.5)*h <= bath[j])
	//			m0[j*nez + k] = vpvs/v;
	//	}
	//	fscanf(v_1, "/n");
	//}
	//for(int i = 1; i < nx-1; ++i)
	//	for(int j = 0; j < ny-1; ++j)
	//		for(int k = 0; k < nz-1; ++k)
	//			m0[i*Neyz + j*nez + k] = m0[j*nez + k];
	//fclose(v_1);
//

}

void TT3D::setVref(double* v)
{
  for(int i = 0; i < Nele; ++i)
    m[i] = v[i];
}









///// TEMP //////////////////////////////////

void TT3D::saveNele(char* fname, double* m_, bool vel)
{
	FILE* f = fopen(fname, "w");
	fprintf(f, "%d %d %d\n", nx-1, ny-1, nz-1);
	for(int i = 0; i < nx-1; ++i)
		for(int k = 0; k < nz-1; ++k)
		{
			for(int j = 0; j < ny-1; ++j)
			{
				double a = m_[i*Neyz + j*nez + k];
				if(vel)
				{
					a = 1./a;
				}
				fprintf(f, "%g ", a);
			}
			fprintf(f, "\n");
		}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n saveNele: volume saved to file %s\n", fname);
}
void TT3D::sliceNele(char* fname, bool vel)
{
	FILE* f = fopen(fname, "w");
	int i = nx/2-1;
	for(int k = 0; k < nez; ++k)
	{
		for(int j = 0; j < ny-1; ++j)
		{
			double a = m[i*Neyz + j*nez + k];
			if(vel)
			{
				a = 1./a;
			}
			fprintf(f, "%f ", a);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n sliceNele: slice saved to file %s\n", fname);
}
void TT3D::sliceNnode(char* fname, float* v, bool vel)
{
	FILE* f = fopen(fname, "w");
	int i = nx/2;
	for(int k = 0; k < nz; ++k)
	{
		for(int j = 0; j < ny; ++j)
		{
			double a = v[i*Nyz + j*nz + k];
			if(vel)
				a = h/a;
			fprintf(f, "%f ", a);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n sliceNnode: slice saved to file %s\n", fname);
}

void TT3D::saveT(char* fname)
{
	FILE* f = fopen(fname, "w");
	//fprintf(f, "sID	rID	sX	rX	Tobs	Tcal\n");
	for(int k = 0; k < Nt; ++k)
		fprintf(f, "%d	%d	%f	%f	%f	%f\n", srcID[k], rcvID[k], src[srcID[k]].z, rcv[rcvID[k]].z, Tobs[k], Tcal[k]);
	fclose(f);

	jiout << sprintf(jiout.buf, "\n Traveltimes saved to file %s\n", fname);
}

void TT3D::saveSD(char* fname)
{
	FILE* f = fopen(fname, "w");
	fprintf(f, "ID,RGridE,RGridN,RElev,SGridE,SGridN,SElev,WaveID,TimeObs,TimeErr,TimeMod\n");

	for(int k = 0; k < Nt; ++k)
	{
		float rx = R0.x + rcv[rcvID[k]].x;
		float ry = R0.y + rcv[rcvID[k]].y;
		float rz = R0.z - rcv[rcvID[k]].z;
		float sx = R0.x + src[srcID[k]].x*h;
		float sy = R0.y + src[srcID[k]].y*h;
		float sz = R0.z - src[srcID[k]].z*h;
		fprintf(f, "%d,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n", k,ry,rx,rz,sy,sx,sz,0,Tobs[k],1./sqrt(Wd[k]),Tcal[k]);
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n Traveltimes saved to file %s\n", fname);
}

void TT3D::saveTcal(char* fname, double err)
{
	FILE* f = fopen(fname, "w");
// 	fprintf(f, "NS,XS,ZS,NR,XR,ZR,T,ERR\n");

	for(int k = 0; k < Nt; ++k)
	{
	      double T = Tcal[k] + GaussPseudoRand(0., err);
// 	      float sx = R0.x + src[srcID[k]].x*h;
	      float sy = R0.y + src[srcID[k]].y*h;
	      float sz = R0.z - src[srcID[k]].z*h;
// 	      float rx = R0.x + rcv[rcvID[k]].x;
	      float ry = R0.y + rcv[rcvID[k]].y;
	      float rz = R0.z - rcv[rcvID[k]].z;
	      fprintf(f, "%d,%.3f,%.3f,%d,%.3f,%.3f,%.3f,%.3f\n", srcID[k],sy,sz,rcvID[k],ry,rz,T,err);
	}
	fclose(f);

	jiout << sprintf(jiout.buf, "\n Traveltimes saved to file %s\n", fname);
}

double GaussPseudoRand(double mean, double sd)
{
// 	const double PI = 3.141592653589793;
	const double eps = 1.e-30;

	double u1 = double(rand())/RAND_MAX + eps;
	double u2 = double(rand())/RAND_MAX + eps;

	return sqrt(-2.*log(u1)) * cos(2.*PI*u2) * sd + mean;
}


