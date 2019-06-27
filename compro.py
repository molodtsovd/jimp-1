from ctypes import *
import numpy as np

NUMPHA_TD =21
MAXENDMEM_TD =7
MAX_NAMLEN =10
MAX_NOX =15
MAX_NPHASES =10

class material():
    def __init__(self, name="vacuum", Q=0., k=0., Nphases=0, Noxides=0):
        self.name = name
        self.Q = Q
        self.k = k
        self.Nphases = Nphases
        self.Noxides = Noxides
        self.namphatxt = np.zeros((MAX_NPHASES, MAX_NAMLEN), dtype=c_char)
        self.phacomtxt = np.zeros((MAX_NPHASES, MAX_NOX+4))
        self.endmemrat = np.zeros((MAX_NPHASES, MAXENDMEM_TD))
        self.molref = np.zeros((MAX_NPHASES,))
        self.volref = np.zeros((MAX_NPHASES,))
        
    def load(self, fname, Nphases, Noxides):
        self.name = fname
        self.Nphases = Nphases
        self.Noxides = Noxides
        self.namphatxt = np.genfromtxt(fname, skip_header=9, 
                                       usecols = 0, dtype = c_char_p, 
                                       max_rows = Nphases)
        self.phacomtxt = np.genfromtxt(fname, skip_header=9, 
                                       usecols = range(1,1+Noxides+4), 
                                       max_rows = Nphases)        
        # parse endmem ratios
        f = open(fname, 'r')       
        for i in range(12 + Nphases):
            f.readline()
        for i in range(Nphases):
            byte = f.read(1)
            if byte == 'M':
                break
            j = 0
            while byte != '\n':
                byte = f.read(1)
                if(byte == ':'):
                    buf = ''
                    while byte != '\n':
                        byte = f.read(1)
                        if byte == ',':
                            break
                        buf += byte
                    self.endmemrat[i, j] = float(buf)
                    j += 1
        f.close()

        self.molref = self.phacomtxt[:,2].copy()
        self.volref = self.phacomtxt[:,1].copy()
        for j in range(self.Nphases):
            if np.all(self.endmemrat[j,:] == 0.):
                self.endmemrat[j,0] = 1.
    
    def print(self, filename=None):
        f = None
        if filename != None:
            f = open(filename, 'w')
        
        print('\n'+self.name, file=f)
        print('Q = %g W/m^3' % self.Q, file=f)
        print('k = %g W/(m*K)' % self.k, file=f)
        print('Nphases = %d' % self.Nphases, file=f)
        print('Noxides = %d' % self.Noxides, file=f)
        print('\nPhase\tvol %', file=f)
        for i in range(self.Nphases):
            print('%s\t%g' % (bytes(self.namphatxt[i]).decode('ascii'), 
                              self.phacomtxt[i, 1]), file=f)
#        print(self.namphatxt, file=f)
#        print(self.phacomtxt, file=f)
#        print(self.endmemrat, file=f)
#        print(self.molref, file=f)
#        print(self.volref, file=f)
        ox = self.getOxMolPerc()
        if self.Noxides == 8:
            labels = ['H2O','Al2O3','FeO','MgO','CaO','Na2O','K2O','SiO2']
        if self.Noxides == 7:      
            labels = ['SiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O']
        print('\nOxide\tmol %', file=f)
        for i in range(self.Noxides):
            print('%s\t%.2f' % (labels[i],ox[i]), file=f)
        if filename != None:
            f.close()
        
    def getOxMolPerc(self):
        phamol = (self.phacomtxt[:,1]-self.volref)*self.molref/self.volref+self.molref
        phamol /= np.sum(phamol) * np.sum(self.phacomtxt[:,4:4+self.Noxides], axis=1)
        return np.sum(self.phacomtxt[:,4:4+self.Noxides]*phamol[:,np.newaxis], axis=0)*100.
        
    def c_(self):
        mc = c_material()
#        for i in range(len(self.name)):
#            mc.name[i] = c_byte(self.name[i].encode())
        mc.name = 0
        if self.name == 'mantle':
            mc.name = 1
        mc.Q = self.Q
        mc.k = self.k
        mc.Nphases = self.Nphases
        mc.Noxides = self.Noxides
        
        # NOTE: pad namphatxt with whitespaces - 
        # on fortran side there is no null-termination and trim() is used
        for i in range(MAX_NPHASES * MAX_NAMLEN):
            mc.namphatxt[i] = 0x20
                
        for i in range(self.Nphases):
            for j in range(len(self.namphatxt[i])):
                mc.namphatxt[i*MAX_NAMLEN + j] = self.namphatxt[i][j]
            for j in range(self.Noxides+4):
                mc.phacomtxt[i*(self.Noxides+4) + j] = self.phacomtxt[i, j]
            for j in range(MAXENDMEM_TD):
                mc.endmemrat[i*MAXENDMEM_TD + j] = self.endmemrat[i, j]

        mc.molref = (c_double * MAX_NPHASES)(*self.molref)
        mc.volref = (c_double * MAX_NPHASES)(*self.volref)
        
        return mc

class c_material(Structure):
  _fields_ = [
#    ("name", c_byte * 128),
    ("name", c_int), # 1=mantle, 0=water/crust/...
    ("Q", c_double),
    ("k", c_double),
    ("Nphases", c_int),
    ("Noxides", c_int),
    ("namphatxt", c_byte * (MAX_NAMLEN * MAX_NPHASES)),
    ("phacomtxt", c_double * ((MAX_NOX+4) * MAX_NPHASES)),
    ("endmemrat", c_double * (MAXENDMEM_TD * MAX_NPHASES)),  
    ("molref", c_double * MAX_NPHASES),
    ("volref", c_double * MAX_NPHASES)
  ]

def cmprint(mc):
#    print(mc.name)
    print(mc.Q)
    print(mc.k)
    print(mc.Nphases)
    print(mc.Noxides)
    N = mc.Nphases
    Nox = mc.Noxides
    for i in range(N):
        for j in range(MAX_NAMLEN):
            print(chr(mc.namphatxt[i*MAX_NAMLEN + j]), end='')
        print('\n')
    for i in range(N):
        for j in range(Nox+4):
            print(mc.phacomtxt[i*(Nox+4) + j], end=' ')
        print('\n')
    for i in range(N):
        for j in range(MAXENDMEM_TD):
            print(mc.endmemrat[i*MAXENDMEM_TD + j], end=' ')
        print('\n')
    for i in range(N):
        print(mc.molref[i], end=' ')
    print('\n')
    for i in range(N):
        print(mc.volref[i], end=' ')
    print('\n')
   

lib = cdll.LoadLibrary('./compro/compro.so')

lib.getvel.argtypes = [POINTER(c_material), POINTER(c_int), POINTER(c_double),
                       POINTER(c_double), POINTER(c_double), 
                       c_double, c_int, c_int, c_int, POINTER(c_double),
                        # for built-in heat solver:
                	   c_int, c_double, c_double, POINTER(c_byte), c_double, c_double]
lib.getvel.restype = c_int

def getvel(layers, layerID, T, h, vp, den, serp,
           updateT, Ttop, Tbot, Tmask, relTol=1.e-7): 
    nex = vp.shape[2]
    ney = vp.shape[1]
    nez = vp.shape[0]
    Nele = nex*ney*nez


    l = (c_material * len(layers))()
    for i in range(len(layers)):
        l[i] = layers[i].c_()
        #cmprint(l[i])

    l = cast(l, POINTER(c_material))
    lid = (c_int * Nele)(*layerID.ravel('F'))
    t = (c_double * Nele)(*T.ravel('F'))
    v = (c_double * Nele)(*vp.ravel('F'))
    d = (c_double * Nele)(*den.ravel('F'))
    s = (c_double * Nele)(*serp.ravel('F'))

    omega = 2. / (1. + np.sin(np.pi/nez)) # 1.97

#    import time
#    start = time.time()

    # h: km->m
    lib.getvel(l, lid, t, v, d, h*1.e3, nex, ney, nez, s,
               updateT, Ttop, Tbot, Tmask, relTol, omega)

#    print("compro_1 took %f ms" % ((time.time() - start) * 1000.))
    
    return [np.reshape(v, vp.shape, 'F'), np.reshape(d, den.shape, 'F'),
            np.reshape(t, T.shape, 'F')]
