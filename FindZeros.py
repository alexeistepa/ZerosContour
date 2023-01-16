import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from sympy import *

def ver_cont(z0,h,g):
    """
    Contour integral of complex function g from the complex float z0 to the point z0 + 1j*h where h is a real float
    """
    def realfixreal(y): return g(z0 + 1j*y).real
    def imagfixreal(y): return g(z0 + 1j*y).imag
    return 1j*(quad(realfixreal,0,h,epsabs = 0.05, limit=200)[0] + 1j*quad(imagfixreal,0,h,epsabs = 0.1, limit=200 )[0]) # SPEED: Can modify tolerence limit here

def hor_cont(z0,h,g):
    """
    Contour integral of complex function g from the complex float z0 to the point z0 + h where h is a real float
    """
    def realfiximag(x): return g(z0 + x).real
    def imagfiximag(x): return g(z0 + x).imag
    return quad(realfiximag,0,h,epsabs = 0.05, limit=200)[0] + 1j*quad(imagfiximag,0,h,epsabs = 0.1, limit=200)[0] # SPEED: Can modify tolerence limit here

def Nzeros(b,w,h,g):
    '''
    Given complex float b, w,h > 0 and g = f'/f calculates N - P
    where N is the float of zeros of f in (Re(b),Re(b)+w)xi(Im(b),Im(b)+h) and P is the float of poles in that region.
    '''
    cont = hor_cont(b,w,g) + ver_cont(b + w,h,g) + hor_cont(b + w + 1j*h,-w,g) + ver_cont(b + 1j*h,-h,g)
    return round((cont/(2*np.pi*1j)).real)

def FindZeros(b,w,h,g,f2,minit = 10,maxit = 100):
    '''
    Consider a holomorphic function f.
    Given complex float b, positive floats w and h, the function g = f'/f and the function f2 = |f|^2,
    FindZeros returns approximations for the zeros of f in the region (Re(b),Re(b)+w)xi(Im(b),Im(b)+h).
    '''
    corners = [[b]]
    n = 1
    stop = False
    while n <= maxit and stop == False:
        tmp = []
        nums = []
        for btmp in corners[n - 1]:
            N1 = Nzeros(btmp,w/(2**n),h/(2**n),g)
            if N1 >= 1:
                tmp.append(btmp)
                nums.append(N1)
            N2 = Nzeros(btmp + w/(2**n),w/(2**n),h/(2**n),g)
            if N2 >= 1:
                tmp.append(btmp + w/(2**n))
                nums.append(N2)
            N3 = Nzeros(btmp + 1j*h/(2**n),w/(2**n),h/(2**n),g)
            if N3 >= 1:
                tmp.append(btmp + 1j*h/(2**n))
                nums.append(N3)
            N4 = Nzeros(btmp + w/(2**n) + 1j*h/(2**n),w/(2**n),h/(2**n),g)
            if N4 >= 1:
                tmp.append(btmp + w/(2**n) + 1j*h/(2**n))
                nums.append(N4)
        if sum(nums) == len(nums) and n >= minit:
            stop = True
        corners.append(tmp)
        n = n + 1
    guesses = np.array(corners[-1]) + w/(2**n) + 1j*h/(2**n)
    return guesses

