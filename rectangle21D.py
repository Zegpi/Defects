import sys
import getopt
import numpy as np
from igakit.cad import *
from igakit.io import *
from math import sqrt
 
try:
	fileName=sys.argv[1]
except:
	fileName="geometry21D"

fileName=fileName+".dat"

try:
	LargoX=float(sys.argv[2])
except:
	LargoX=1.0

try:
	nx=int(sys.argv[4])
except:
	nx=101

Lx = LargoX      	#  Horizontal length of the rectangle
Nx = nx    			#  number of elements in the horizontal direction (x direction)
p = 2         		#  Order of the basis functions
C = 0         		#  Inter-element continuity


points = np.zeros((2,2), dtype='d')
points[0,0] = -Lx/2.0
points[1,0] = +Lx/2.0

geom = linear(points)

geom = refine( geom, factor=1, degree=p)

geom = refine( geom, factor=Nx )

#geom.unclamp(0)
#geom.unclamp(1)


PetIGA().write(fileName, geom, nsd=1)

if 0:
    from igakit.plot import plt
    plt.figure()
    plt.cplot(geom)
    plt.kplot(geom)
    plt.show()
