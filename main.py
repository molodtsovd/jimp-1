import numpy as np
from scipy import interpolate
import ctypes as ct
from copy import deepcopy
from grav import *
from tt3d import *
from compro import *
import ownplt
import time
from mpi4py import MPI


class compot():

    def __init__(self, h=0.5):
        
        # build FD mesh
        
        self.h = h
        self.R0 = R0 = point3d(0, 0, 0)
        self.R1 = R1 = point3d(h, 231, -41)#[km]
        
        self.nex = int(np.ceil((R1.x - R0.x) / h))
        self.ney = int(np.ceil((R1.y - R0.y) / h))
        self.nez = int(np.ceil((R0.z - R1.z) / h))
        nex,ney,nez = self.nex,self.ney,self.nez
        
        self.Yh = np.linspace(R0.y + h/2, R1.y - h/2, ney)
        self.Zh = np.linspace(R0.z - h/2, R1.z + h/2, nez)
        
        self.vp = np.zeros((nez, ney, nex))
        self.den = np.zeros((nez, ney, nex))
        self.serpwFD = np.zeros((nez, ney, nex))
        self.layerID = np.zeros((nez, ney, nex), dtype=int)
        self.T = np.zeros((nez, ney, nex))
        self.Tmask = np.zeros((nez, ney, nex), dtype=bool)
        
        self.Ttop = 300.
        self.Tbot = 1500.
        
        # build FD mesh and load seismic data in TT3D C++ class
        self.seis = TT3D(R0, R1, h, b"test/rapids4_100.ttime", b"test/rapids4_100.bath")
        
        # build FD mesh and load gravity data in grav C++ class
        self.grav = Grav(R0, R1, h, 0, b"test/synth.gd", b"test/rapids4_100.bath", 
                    b"test/isopar.txt", b"test/synth.elev")
        
        # load tomographic vel
        
        vtom = np.genfromtxt("test/tomo/v5.dat", skip_header=83, max_rows=82)
        htom = 0.5
        y = np.linspace(R0.y + htom/2, R1.y - htom/2, int((R1.y-R0.y)/htom))
        z = np.linspace(R1.z + htom/2, R0.z - htom/2, int((R0.z-R1.z)/htom))
        spl = interpolate.RectBivariateSpline(z, y, vtom)
        self.vp[:,:,0] = spl(np.flip(self.Zh, axis=0), self.Yh)
        
        # density(vp) - Prada et al., 2018:
        self.den = 0.357 + 1.114*self.vp - 0.182*self.vp**2 + 0.010*self.vp**3        
        
        # load interfaces of RAPIDS4 model (O'Reilly et al., 2006)
        
        # bathymetry
        y0, z0 = np.genfromtxt("test/rapids4/in.txt", max_rows=2)
        y0 = np.flip(R1.y - y0,0)
        z0 = np.flip(z0,0)
        # basement top
        y1, z1 = np.genfromtxt("test/rapids4/in.txt", skip_header=6, max_rows=2)
        y1 = np.flip(R1.y - y1,0)
        z1 = np.flip(z1,0)
        # Moho
        y2, z2 = np.genfromtxt("test/rapids4/in.txt", skip_header=8, max_rows=2)
        y2 = np.flip(R1.y - y2,0)
        z2 = np.flip(z2,0)
        # interpolate Moho on the inversion mesh
        y = y1 #np.linspace(R0.y, R1.y, num=25)
        spl = interpolate.PchipInterpolator(y2, z2)
        z2 = spl(y)
        y2 = y
        
        # build layered model
        
        self.Nlayers = 6
        Nlayers = self.Nlayers
        self.layers = np.empty(Nlayers, dtype=object)
        fnames = ['bulk_SH025met_1', 'bulk_TM025met_1', 'bulk_TMMFequ_1']
        Nphases = [0, 0, 6, 7, 6, 0]
        Noxides = [0, 0, 8, 8, 7, 0]
        self.H = np.empty(Nlayers-1, dtype=list)
        
        # mantle model
        self.serpw = 33. * np.exp(-((y-97.5)/20.)**2) # weight %
        self.serpDepth = 3.
        
        self.Hfd = np.zeros((Nlayers+1, ney)) #2D
        self.Hfd[0,:] = np.linspace(R0.y + h/2, R1.y - h/2, ney)
        
        self.layers[0] = material(name='water', k = 0.6)
        self.layers[1] = material(name='sediments', Q = 1.e-6, k = 2.5)
        self.layers[Nlayers-1] = material(name='mantle', Q = 0., k = 5.3)
        # NOTE: removing conductivity contrast between mask and earth to stabilize SOR:
        self.layers[0].k = self.layers[1].k
        
        self.H[0] = []
        for i in range(len(y0)):
            self.H[0].append(point3d(0, y0[i], z0[i]))
        self.H[1] = []
        for i in range(len(y1)):
            self.H[1].append(point3d(0, y1[i], z1[i]))
        
        for i in range(2, Nlayers-1):
            
            self.layers[i] = material(Q = 0.5e-6, k = 3.5)
            self.layers[i].load("compro/%s.txt" % fnames[i-2], Nphases[i], Noxides[i])
            
            self.H[i] = []
            for j in range(len(y)):
                self.H[i].append(point3d(0, y[j], (z2[j] - z1[j])/(Nlayers-3)))
        # modify the center...
        for i in range(11, 13):
            self.H[3][i].z += self.H[2][i].z / 2.
            self.H[4][i].z += self.H[2][i].z / 2.
            self.H[2][i].z = 0.
        for i in range(13, 20):
            self.H[4][i].z += self.H[2][i].z + self.H[3][i].z
            self.H[2][i].z = self.H[3][i].z = 0.
        for i in range(20, 23):
            self.H[3][i].z += self.H[2][i].z / 2.
            self.H[4][i].z += self.H[2][i].z / 2.
            self.H[2][i].z = 0.
        
        # NOTE: making upper crustal layer bulk_SH025met_1 really felsic
        self.layers[2].phacomtxt[:,1] = [0.49,39.98,8.10,0.51,9.74,41.18]                
#        self.layers[2].phacomtxt[:,1] = [0.,44.,0.0,0.,13.,43.]                

        # select 'free' model parameters to invert for
        # NOTE: view() only works if slices are contiguous!
        # TODO: replace by vectors of prior PDF
        self.layers_free = self.layers[2:5].view()
        self.H_free = self.H[2:5].view()
        self.serpw_free = self.serpw[9:27].view()
        # backup copies for rejecting a step
        self.layers_free_back = deepcopy(self.layers_free)
        self.H_free_back = deepcopy(self.H_free)
        self.serpw_free_back = deepcopy(self.serpw_free)
        if not(np.shares_memory(self.layers_free, self.layers) 
               and np.shares_memory(self.H_free, self.H) 
               and np.shares_memory(self.serpw_free, self.serpw)):
            from sys import exit
            exit()


    # map the model to FD mesh
    def map2fd(self, newGeometry=True):
        
        if newGeometry == True:                   
            for i in range(self.Nlayers-1):
                spl = interpolate.PchipInterpolator([p.y for p in self.H[i]], 
                                                    [p.z for p in self.H[i]])
                self.Hfd[i + 1, :] = spl(self.Hfd[0])            
            for i in range(self.nex):
                for j in range(self.ney):
                    l = 0
                    z = self.Hfd[1, j]
                    for k in range(self.nez):
                        while (k + 0.5)*self.h > z and l < self.Nlayers-1:
                            l += 1 
                            z += self.Hfd[l + 1, j]
                        self.layerID[k,j,i] = l

        # mantle serpentinization
        spl = interpolate.PchipInterpolator([p.y for p in self.H[-1]], self.serpw)
        serpFDtop = spl(self.Hfd[0])
        z0 = -self.Hfd[1:self.Nlayers, :].sum(axis=0)
        
        for i in range(self.nex):
            for j in range(self.ney):
                self.serpwFD[:,j,i] = serpFDtop[j] * \
                    np.exp(-abs(self.Zh - z0[j])/self.serpDepth)
        
    
    # set the temperature Dirichlet mask and water density after the 
    # first call of map2fd()
    def setDenTmask(self):   
        self.Tmask[self.layerID == 0] = 1
        self.Tmask = (ct.c_byte * (self.nex*self.ney*self.nez))(*self.Tmask.ravel('F'))
        self.den[self.layerID == 0] = 1.

    def fwd(self, updateTemp):
        # calculate velocity and density
        [self.vp, self.den, self.T] = getvel(self.layers, self.layerID, self.T, self.h, 
            self.vp, self.den, self.serpwFD, int(updateTemp), self.Ttop, self.Tbot, self.Tmask)
        # solve seismic forward
#        start = time.time()
        phi_s, self.RMSs = self.seis.FP(self.vp)
#        print("seis took %f ms" % ((time.time() - start) * 1000.))
        # solve gravity forward
        phi_g, self.RMSg, self.RMSe = self.grav.FP_2D(self.den)
        return (self.RMSs + self.RMSg + self.RMSe) #/ (1./88381+1./201+1./462) #phi_s + phi_g
        
        
    def plotSaveData(self): 
#        # plot T
#        ax = ownplt.plotTemp(self.T, self.h)
#        ownplt.plotH(self.H, self.Hfd, ax)
        # plot vp
        ax = ownplt.plotV2d(self.vp[:,:,0], self.h, v0=1.5, v1=8.2)
#        ownplt.plotH(self.H, self.Hfd, ax) 
        # plot den
        ax = ownplt.plotV2d(self.den[:,:,0], self.h, 
                            v0=1., v1=3.3, dv=0.5, prop='Density $[g/cm^3]$')
#        ownplt.plotH(self.H, self.Hfd, ax)              
        # traveltimes
        self.seis.saveSD(b"test/t.sd")        
        ownplt.plotSD("test/t.sd")
        # gravity        
        self.grav.saveGD(b"test/mod_fwd_flt.gd")
        ownplt.plotGD_2D("test/mod_fwd_flt.gd")
        # elevation
        self.grav.saveElev(b"test/calc.elev")
        ownplt.plotElev("test/calc.elev", "test/rapids4_100.bath", ax)
            
    def plotHistZ(self):
        R0 = [self.R0.y, self.R0.z]
        R1 = [self.R1.y, self.R1.z]
        x = [p.y for p in self.H_free[0]]   
        axes = ownplt.plotHistZ('TMP', 3, R0, R1, 0.5, x)
        for ax in axes:
            ownplt.plotH(self.H, self.Hfd, ax, style=':', width=1.,color='white')

    def saveModel(self, filename):
        f = open(filename, 'a')
        for H in self.H_free:
            np.savetxt(f, [[p.z for p in H]], fmt='%.2f')
        for l in self.layers_free:
            np.savetxt(f, [l.phacomtxt[:,1]], fmt='%.2f')
        np.savetxt(f, [self.serpw_free], fmt='%.2f')
        f.write('\n')
        f.close()

    def saveRMS(self, filename):
        f = open(filename, 'a')
        np.savetxt(f, [[self.RMSs, self.RMSg, self.RMSe, 
                        self.phi_back, self.chainT[self.chainID]]], fmt='%g')           
        f.close()

    def loadModel(self, filename, offset):
        for H in self.H_free:
            a = np.genfromtxt(filename, skip_header=offset, max_rows=1)
            for i in range(len(H)):
                H[i].z = a[i]
            offset += 1                
        for l in self.layers_free:
            a = np.genfromtxt(filename, skip_header=offset, max_rows=1, unpack=True)
            l.phacomtxt[:,1] = a
            offset += 1
        self.serpw_free[:] = np.genfromtxt(filename, skip_header=offset, max_rows=1)

    def randModel(self):
        for H in self.H_free:
            for p in H:
                p.z = self.maxH / len(self.H_free) * np.random.rand() + self.minH
        for l in self.layers_free:
            c = np.random.rand(l.Nphases)
            l.phacomtxt[:,1] = c/sum(c)*100.                
        self.serpw_free[:] = self.maxC * np.random.rand(len(self.serpw_free)) + self.minC
        self.map2fd()

    def stepH(self):
        i = np.random.choice(len(self.H_free))
        j = np.random.choice(len(self.H_free[i]))
        z = self.H_free[i][j].z + self.sigmaH * np.random.randn()
        # NOTE: assuming all lists within H starting from sediments have equal length!
        s = sum(h[j].z for h in self.H[1:])        
        if z < self.minH or z + s - self.H_free[i][j].z > self.maxH:
#        if z < self.minH or z > self.maxH:
            return True
        self.H_free_back[:] = deepcopy(self.H_free)
        self.H_free[i][j].z = z
        self.map2fd()
        return False

    def stepC(self):
        i = np.random.choice(len(self.layers_free) + 1)
        if i == len(self.layers_free):
            j = np.random.choice(len(self.serpw_free))
            c = self.serpw_free[j] + self.sigmaC * np.random.randn()
            if c < self.minC or c > self.maxC:
                return True
            self.serpw_free_back[:] = deepcopy(self.serpw_free)
            self.serpw_free[j] = c
            self.map2fd(False)
        else:
            j1,j2 = np.random.choice(self.layers_free[i].Nphases, size=2, replace=False)
            dc = self.sigmaC * np.random.randn()
            c1 = self.layers_free[i].phacomtxt[j1,1] + dc
            c2 = self.layers_free[i].phacomtxt[j2,1] - dc
            if c1 < self.minC or c1 > self.maxC or c2 < self.minC or c2 > self.maxC:
                return True
            self.layers_free_back[:] = deepcopy(self.layers_free)
            self.layers_free[i].phacomtxt[j1,1] = c1
            self.layers_free[i].phacomtxt[j2,1] = c2
        return False
    
    def swap(self, rank, pair):
        comm = MPI.COMM_WORLD
        if rank == pair[0]:
            buf = np.array([self.chainID, self.phi_back])
            comm.Send(buf, dest=pair[1], tag=0)
        elif rank == pair[1]:
            buf = np.empty(2)
            comm.Recv(buf, source=pair[0], tag=0)               
            if np.random.rand() < np.exp((self.phi_back - buf[1]) * 
                              (1./self.chainT[self.chainID] - 1./self.chainT[int(buf[0])])):
                a = buf[0]
                buf[0] = self.chainID
                self.chainID = int(a)
#DEBUG: count the swap
                f=open('samp/swap','a')
                print(self.i, self.chainT[self.chainID], self.chainT[int(buf[0])], file=f)
                f.close()
##            
            comm.Send(buf[0], dest=pair[0], tag=1)
        if rank == pair[0]:
            buf = np.empty(1)
            comm.Recv(buf, source=pair[1], tag=1)
            self.chainID = int(buf)
    
    def inv(self, Niter, fname_ar='acceptance_rate', fname_samp='sampling', fname_rms='RMS'):
        Nchains = MPI.COMM_WORLD.Get_size()
        rank = MPI.COMM_WORLD.Get_rank()

        self.sigmaH = 5.
        self.minH = 0.
        self.maxH = 40.
        self.sigmaC = 3.75
        self.minC = 0.
        self.maxC = 100.
        
        chainTmax = 1e7
        chainTmin = 1e-7
#        Nchains_T1 = int(0.25*Nchains)
#        self.chainT = np.concatenate([np.ones(Nchains_T1-1), 
#                                      np.geomspace(chainTmin, chainTmax, Nchains-Nchains_T1+1)])
        self.chainT = np.geomspace(chainTmin, chainTmax, Nchains)
        self.chainID = rank
        
        self.Naccepted_H = 0
        self.Naccepted_C = 0
        self.randModel()
        self.phi_back = self.fwd(True) 
        for i in range(1,Niter+1):
            self.i = i
            f=open(fname_ar+str(self.chainID),'w')
            print(2.*self.Naccepted_H/i, 2.*self.Naccepted_C/i, i, file=f)
            f.close()
            
            if rank == 0:
                pair = np.random.choice(Nchains, size=2, replace=False)
            else:
                pair = np.empty(2, dtype=int)
            MPI.COMM_WORLD.Bcast(pair, root=0)
            self.swap(rank, pair)
            
            if i%2:
                if self.stepH():
                    continue
            else:
                if self.stepC():
                    continue
            self.phi = self.fwd(i%2)
            if np.random.rand() > np.exp((self.phi_back - self.phi)/self.chainT[self.chainID]):
                if i%2:
                    self.H_free[:] = deepcopy(self.H_free_back)
                else:
                    self.layers_free[:] = deepcopy(self.layers_free_back)
                    self.serpw_free[:] = deepcopy(self.serpw_free_back)
                self.phi = self.phi_back
                continue
            self.phi_back = self.phi  
            if i%2:
                self.Naccepted_H += 1
            else:
                self.Naccepted_C += 1

            self.saveModel(fname_samp+str(self.chainID))
            self.saveRMS(fname_rms+str(self.chainID))
                                
                
        
########################################################################
# inv
########################################################################
np.random.seed(MPI.COMM_WORLD.Get_rank())

z = compot(0.5)

z.map2fd()
z.setDenTmask()

z.inv(int(1e6), 'samp/acceptance_rate', 'samp/sampling', 'samp/RMS')

