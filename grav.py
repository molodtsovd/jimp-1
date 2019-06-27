import ctypes as ct
#from ctypes import *
import numpy as np

class point3d(ct.Structure):
     _fields_ = [("x", ct.c_double), ("y", ct.c_double), ("z", ct.c_double)]


lib = ct.cdll.LoadLibrary('./grav/grav.so')


class Grav(object):

    def __init__(self, R0, R1, h, nez_air, 
                 d_filename, bath_filename, iso_filename, elev_filename):
        lib.grav_.argtypes = [ct.Structure, ct.Structure, ct.c_double, ct.c_int, 
                              ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.c_char_p]
        lib.grav_.restype = ct.c_void_p

#        lib.FP_mask_.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_byte)]
#        lib.FP_mask_.restype = ct.c_double

        lib.FP_2D_.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), 
                               ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
        lib.FP_2D_.restype = ct.c_double

        lib.FP_SOR_.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), 
                                ct.c_double, ct.c_double]
        lib.FP_SOR_.restype = ct.c_double

        lib.saveGD_.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.saveGD_.restype = ct.c_void_p

        lib.saveElev_.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.saveElev_.restype = ct.c_void_p

        self.obj = lib.grav_(R0, R1 ,h, nez_air, 
                             d_filename, bath_filename, iso_filename, elev_filename)

#    def FP_mask(self, den, mask):
#      v = (ct.c_double * den.size)(*den.ravel('F'))
#      m = (ct.c_byte * mask.size)(*mask.ravel('F'))
#      return lib.FP_mask_(self.obj, v, m)

    def FP_2D(self, den):
      v = (ct.c_double * den.size)(*den.ravel('F'))
      rmsg = ct.c_double()
      rmse = ct.c_double()
      phi = lib.FP_2D_(self.obj, v, ct.byref(rmsg), ct.byref(rmse))
      return phi, rmsg.value, rmse.value
    
    
    def FP_SOR(self, den, h, relTol=1.e-7):
        
      nez = den.shape[0]
      omega = 1.97 #2. / (1. + np.sin(np.pi/nez)) # 1.97
      # gravitational constant: g[mGal] = G * density[g/cm3] * distance[km]
      c = -h**2 * 4. * np.pi * 6.674;
      v = (ct.c_double * den.size)(*(den.ravel('F')*c))
      return lib.FP_SOR_(self.obj, v, relTol, omega)


    def saveGD(self, fname):
      return lib.saveGD_(self.obj, fname)

    def saveElev(self, fname):
      return lib.saveElev_(self.obj, fname)

