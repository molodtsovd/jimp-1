import ctypes as ct
#from ctypes import *

class point3d(ct.Structure):
     _fields_ = [("x", ct.c_double), ("y", ct.c_double), ("z", ct.c_double)]


lib = ct.cdll.LoadLibrary('./tt3d/tt3d.so')


class TT3D(object):

    def __init__(self, R0, R1, h, ft_filename, bath_filename):
        lib.TT3D_.argtypes = [ct.Structure, ct.Structure,
                              ct.c_double, ct.c_char_p, ct.c_char_p]
        lib.TT3D_.restype = ct.c_void_p

        lib.FP_.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
        lib.FP_.restype = ct.c_double

        lib.saveSD_.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.saveSD_.restype = ct.c_void_p

        lib.saveTcal_.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double]
        lib.saveTcal_.restype = ct.c_void_p

        lib.setVref_.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double)]
        lib.setVref_.restype = ct.c_void_p

        self.obj = lib.TT3D_(R0, R1 ,h, ft_filename, bath_filename)

    def FP(self, vel):
      v = (ct.c_double * vel.size)(*vel.ravel('F'))
      rms = ct.c_double()
      phi = lib.FP_(self.obj, v, ct.byref(rms))
      return phi, rms.value

    def saveSD(self, fname):
      return lib.saveSD_(self.obj, fname)

    def saveTcal(self, fname, err):
      return lib.saveTcal_(self.obj, fname, err)

    def setVref(self, vel):
      
      v = (ct.c_double * vel.size)(*vel.ravel('F'))
      return lib.setVref_(self.obj, v)

