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
	fileNameInput="sigma-2d*"

fName1="./results/"+fileNameInput+".dat"
fName2="../Results/"+fileNameInput+".dat"

Lista1=glob.glob(fName1)
Lista2=glob.glob(fName2)

if len(Lista2)>0:
	fileName1=fName2
else:
	fileName1=fName1

nombre=fileNameInput.split("*")[0]

try:
	fileNameInput2=sys.argv[2]
except:
	fileNameInput2="stress-2d*"

fName3="./results/"+fileNameInput2+".dat"
fName4="../Results/"+fileNameInput2+".dat"

Lista1=glob.glob(fName3)
Lista2=glob.glob(fName4)

if len(Lista2)>0:
	fileName2=fName4
else:
	fileName2=fName3

try:
	meshName=sys.argv[3]
except:
	meshName="geometry"

meshName="./"+meshName+".dat"
# read in discretization info and potentially geometry
nrb = PetIGA().read(meshName)

Lista1=glob.glob(fileName1)
Lista2=glob.glob(fileName2)

num1=len(Lista1)
num2=len(Lista2)

matA=numpy.ndarray([2,2])
matB=numpy.ndarray([2,2])

if num1==num2:

	for x in range(0, int(num1)):
		try:
			sol1=PetIGA().read_vec(Lista1[x],nrb)		#sol is of type numpy.ndarray, la forma del vector es sol.shape
			sol2=PetIGA().read_vec(Lista2[x],nrb)
		except:
			print('Mismatch in vector lengths, wrong mesh?, file1='+sys.argv[1]+' file2='+sys.argv[2])
			sys.exit()

		diff=numpy.ndarray([sol1.shape[0],sol1.shape[1]])

		for i in range(sol1.shape[0]):
			for j in range(sol1.shape[1]):

				matA[0][0]=sol2[i][j][0]-sol1[i][j][0]							#Sigma-Sigma_exact
				matA[0][1]=sol2[i][j][1]-sol1[i][j][1]
				matA[1][0]=sol2[i][j][2]-sol1[i][j][2]
				matA[1][1]=sol2[i][j][3]-sol1[i][j][3]

				#matB[0][0]=matA[0][0]*matA[0][0]+matA[1][0]*matA[1][0]
				#matB[0][1]=matA[0][0]*matA[0][1]+matA[1][0]*matA[1][1]
				#matB[1][0]=matA[0][1]*matA[0][0]+matA[1][1]*matA[1][0]
				#matB[1][1]=matA[0][1]*matA[0][1]+matA[1][1]*matA[1][1]

				val1=numpy.sqrt(matA[0][0]*matA[0][0]+matA[0][1]*matA[0][1]+matA[1][0]*matA[1][0]+matA[1][1]*matA[1][1])

				#lam,v=numpy.linalg.eig(matB);
				#eig1=numpy.max(lam);

				matA[0][0]=sol2[i][j][0]										#Sigma_exact
				matA[0][1]=sol2[i][j][1]
				matA[1][0]=sol2[i][j][2]
				matA[1][1]=sol2[i][j][3]

				#matB[0][0]=matA[0][0]*matA[0][0]+matA[1][0]*matA[1][0]
				#matB[0][1]=matA[0][0]*matA[0][1]+matA[1][0]*matA[1][1]
				#matB[1][0]=matA[0][1]*matA[0][0]+matA[1][1]*matA[1][0]
				#matB[1][1]=matA[0][1]*matA[0][1]+matA[1][1]*matA[1][1]

				#lam,v=numpy.linalg.eig(matB);
				#eig2=numpy.max(lam);

				val2=numpy.sqrt(matA[0][0]*matA[0][0]+matA[0][1]*matA[0][1]+matA[1][0]*matA[1][0]+matA[1][1]*matA[1][1])

				#diff[i][j]=eig1/eig2
				diff[i][j]=val1/val2

		outfile = Lista1[x][:-4] + " - Diff.vtk"

		num = numpy.round(numpy.sqrt(sol1.size/4.0),decimals=0)

		uniform = lambda U: linspace(U[0], U[-1], 1*num)

		# write a binary VTK file
		VTK().write(outfile,       			# output filename
	            nrb,                    	# igakit NURBS object
	            fields=diff,             	# field is the numpy array to plot 
	            sampler=uniform,        	# specify the function to sample points
	            scalars={nombre+'Diff(1)':0}) # adds a scalar plot to the VTK file

else:
	print('Number of files to compare is different')