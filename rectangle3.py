import sys
import getopt
import numpy as np
from igakit.cad import *
from igakit.io import *
from math import sqrt
 
try:
	fileName=sys.argv[1]
except:
	fileName="geometry3"

fileName=fileName+".dat"

try:
	LargoX=float(sys.argv[2])
except:
	LargoX=1.0

try:
	LargoY=float(sys.argv[3])
except:
	LargoY=1.0

try:
	nx=int(sys.argv[4])
except:
	nx=101

try:
	ny=int(sys.argv[5])
except:
	ny=101

Lx = LargoX     #  Horizontal length of the rectangle
Ly = LargoY     #  Vertical length of the rectangle   
Nx = nx    		#  number of elements in the horizontal direction (x direction)
Ny = ny    		#  number of elements in the vertical direction (y direction)
p = 3         	#  Order of the basis functions
C = 0         	#  Inter-element continuity

#geom = bilinear([[[ 0.,  0.],[ 0., Ly]],[[ Lx,  0.],[ Lx, Ly]]])
geom = bilinear([[[ -Lx/2.0,  -Ly/2.0],[ -Lx/2.0, Ly/2.0]],[[ Lx/2.0, -Ly/2.0],[ Lx/2.0, Ly/2.0]]])

geom = refine( geom, factor=(1, 1), degree=(p, p))

geom = refine( geom, factor=(Nx, Ny) )

#geom.unclamp(0)
#geom.unclamp(1)

#nds seems to be number of spatial dimensions
PetIGA().write(fileName, geom, nsd=2)

if 0:
    from igakit.plot import plt
    plt.figure()
    plt.cplot(geom)
    plt.kplot(geom)
    plt.show()
