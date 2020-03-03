"""
This python script shows the usage of igakit
(https://bitbucket.org/dalcinl/igakit) to post-process results
obtained using the demo code:

When running the C-code in PetIGA, two data files are
generated. One contains the geometry and discretization information,
the other is the solution vector.
"""

from igakit.io import PetIGA,VTK
import numpy
from numpy import linspace
import glob
import sys

# read in discretization info and potentially geometry
try:
	fileNameInput=sys.argv[1]
except:
	fileNameInput="AL-2d*"

fileName1="./results/"+fileNameInput+".dat"
fileName2="../Results/"+fileNameInput+".dat"

Lista1=glob.glob(fileName1)
Lista2=glob.glob(fileName2)

if len(Lista2)>0:
	fileName=fileName2
else:
	fileName=fileName1

nombre=fileNameInput.split("*")[0]

try:
	meshName=sys.argv[2]
except:
	meshName="geometry"

meshName="./"+meshName+".dat"
# read in discretization info and potentially geometry
nrb = PetIGA().read(meshName)

#for infile in glob.glob("results/poisson2d*.dat"):
for infile in glob.glob(fileName):
# read in solution vector as a numpy array
	try:
		sol = PetIGA().read_vec(infile,nrb)
	except:
		print('Mismatch in vector lengths, wrong mesh?, file='+sys.argv[1])
		sys.exit()

	outfile = infile[:-4] + ".vtk"				#string[:-4] cuts the last 4 characters from the string, in this case ".dat"

	num = numpy.round(numpy.sqrt(sol.size/2.0),decimals=0)
	#print(num)

	# write a function to sample the nrbs object
	uniform = lambda U: linspace(U[0], U[-1], 1*num)


# write a binary VTK file
	VTK().write(outfile,       # output filename
            nrb,                    # igakit NURBS object
            fields=sol,             # sol is the numpy array to plot 
            sampler=uniform,        # specify the function to sample points
            scalars={nombre+'(1)':0,nombre+'(2)':1}) # adds a scalar plot to the VTK file