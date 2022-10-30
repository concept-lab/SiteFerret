import os
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

libname = os.path.abspath("libCfunc.so")

Cfunc = ctypes.CDLL(libname)
pdist_C = Cfunc.pdist
pdist_C.argtypes=[ctypes.POINTER(ctypes.c_void_p),ctypes.POINTER(ctypes.c_void_p),ctypes.POINTER(ctypes.c_void_p),ctypes.c_int, ctypes.c_int,ctypes.c_int]
pdist_C.restype = None

volumeArvo_C=Cfunc.ARVOLUME_C
volumeArvo_C.argtypes=[ctypes.POINTER(ctypes.c_void_p),ctypes.POINTER(ctypes.c_void_p),ctypes.c_int]
volumeArvo_C.restype = None #error message: 0 success , -1 fail

def Pdist_C(X1,X2):
    """
    Pdist wrapper using shared C library for improved performance
    """
    m1,n1 = X1.shape
    m2,n2 = X2.shape
    if ((n1!=n2) and (n1!=3)):
        raise ValueError("Both imput arrays must have dimension equal to 3")
    # Convert to pointers for C
    X1 = np.ascontiguousarray(X1,dtype = ctypes.c_double)
    X2 = np.ascontiguousarray(X2,dtype = ctypes.c_double)
    pX1 = X1.ctypes.data_as(ctypes.POINTER(ctypes.c_void_p))
    pX2 = X2.ctypes.data_as(ctypes.POINTER(ctypes.c_void_p))
    ######

    if(np.array_equal(X1,X2)):
        flag = 1
        #scenario for comparison among same radius probes
        m = m1
        #fix correct dimension for output
        d = np.empty(m*(m-1)//2)
        #Again, broadcast to correct format for ctypes
        d = np.ascontiguousarray(d,dtype = ctypes.c_double)
        #pointer to d, the return variable
        pd = d.ctypes.data_as(ctypes.POINTER(ctypes.c_void_p))
        m = ctypes.c_int(m1)
        ####
        # Call to C func 
        pdist_C(pd,pX1,pX1,ctypes.c_int(1),m,m)
    else:
        flag = 0
        d = np.empty(m1*m2)
        d = np.ascontiguousarray(d,dtype = ctypes.c_double)
        pd = d.ctypes.data_as(ctypes.POINTER(ctypes.c_void_p))
        m1=ctypes.c_int(m1)
        m2=ctypes.c_int(m2)
        pdist_C(pd,pX1,pX2,ctypes.c_int(0),m1,m2)
    
    return d,flag;


def getIndex(flag,index,m):
    """
    Returns indexes of original observations corresponding to queried distance index
    Careful, if flag 0 m = size of second set of observations (m2 above)
    """
    if(flag==1):
        #pairwise distances among same set
        if (index>=m*(m-1)//2):
            exit("got too large index")
        b = 1 -2*m 
        i = int(np.floor((-b - np.sqrt(b**2 - 8*index))/2))
        j = int(index + i*(b + i + 2)/2 + 1)
        # _static_previous_m = m #not very nice since it depends on ordering in Pdist..
    elif(flag==0):
        #pairwise distances among different sets
        # if (index>=m*_static_previous_m):
        #     exit("got too large index")
        i = index//m 
        j = index%m 
    else:
        raise ValueError("flag must be provided with values 1 or 0")

    return (i,j);