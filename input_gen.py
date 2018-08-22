import os
import sys
import numpy
from numpy import *
import math
import string
#import pymatgen as mg
n = input("enter no of atoms:")
ifile = open("positions","r")
ifile1 = open("cell_params","r")
mass_dict = {'Ba':'137.327', 'Mo':'95.94', 'O':'15.9994', 'Ag':'107.8682', 'Nb':'92.90638', 'Ti':'47.867', 'Sr':'87.62', 'Hf':'178.49', 'V':'50.9415', 'Ta':'180.94788', 'Ca':'40.078', 'Zr':'91.224', 'K':'39.0983', 'Pb':'207.2', 'Na':'22.98976928', 'S':'32.065', 'Se':'78.96', 'Te':'127.6'}
pp_dict = {'Ba':'Ba.pbe-mt_fhi.UPF', 'Mo':'Mo.pbe-mt_fhi.UPF', 'O':'O.pbe-mt.UPF', 'Ag':'Ag.pbe-mt_fhi.UPF', 'Nb':'Nb.pbe-mt_fhi.UPF', 'Ti':'Ti.pbe-mt_fhi.UPF', 'Sr':'Sr.pbe-mt_fhi.UPF', 'Hf':'Hf.pbe-mt_fhi.UPF', 'V':'V.pbe-mt_fhi.UPF', 'Ta':'Ta.pbe-mt_fhi.UPF', 'Ca':'Ca.pbe-mt_fhi.UPF', 'Zr':'Zr.pbe-mt_fhi.UPF', 'K':'K.pbe-mt_fhi.UPF', 'Pb':'Pb.pbe-mt_fhi.UPF', 'Na':'Na.pbe-mt_fhi.UPF', 'S':'S.pbe-mt_fhi.UPF', 'Se':'Se.pbe-mt_fhi.UPF', 'Te':'Te.pbe-mt_fhi.UPF'}
A_ele_list = ["Ca", "Sr", "Ba", "Pb", "Ag", "Na", "K"]
B_ele_list = ["Ti", "Zr", "Hf", "Nb", "Ta", "Mo", "V"]
X_ele_list = ["S", "Se", "Te"]
x = zeros(n,float)
y = zeros(n,float)
z = zeros(n,float)
cellparamx = zeros(3,float)
cellparamy = zeros(3,float)
cellparamz = zeros(3,float)
charge = zeros(n,float)
mass = zeros(n,float)
k = 0
print A_ele_list[0], pp_dict[A_ele_list[0]], mass_dict[X_ele_list[0]]
for line in ifile:
    items = str.split(line)
    x[k] = float(items[0])
    y[k] = float(items[1])
    z[k] = float(items[2])
    k = k + 1
file.close(ifile)
k = 0
for line in ifile1:
    items = str.split(line)
    cellparamx[k] = float(items[0])
    cellparamy[k] = float(items[1])
    cellparamz[k] = float(items[2])
    k = k + 1
file.close(ifile1)
k = 0
for i in range(7):
    for j in range(7):
        for k in range(3):
            ofile = open(A_ele_list[i]+B_ele_list[j]+X_ele_list[k]+str(3)+".pwscf.in","w")
            fname = A_ele_list[i]+B_ele_list[j]+X_ele_list[k]+str(3)
            file.write(ofile, "&CONTROL"+"\n")
            file.write(ofile, "  calculation = 'vc-relax',"+"\n")
            file.write(ofile, "  forc_conv_thr = 0.001,"+"\n")
            file.write(ofile, "  tstress = TRUE,"+"\n")
            file.write(ofile, "  outdir = '/scratch/conte/k/kgudavis/from_python_code/quantum_espresso/ABX3_chalcogens/cubic/out/"+A_ele_list[i]+B_ele_list[j]+X_ele_list[k]+str(3)+"_cubic/',"+"\n")
            file.write(ofile, "  pseudo_dir = '/scratch/conte/a/azadoks/sq2qe/pseudo/pbe/',"+"\n")
            file.write(ofile, "  verbosity = 'high',"+"\n")
            file.write(ofile, "/"+"\n")
            file.write(ofile, "&SYSTEM"+"\n")
            file.write(ofile, "  degauss = 0.003,"+"\n")
            file.write(ofile, "  ecutwfc = 70.0,"+"\n")
            file.write(ofile, "  ibrav = 0,"+"\n")
            file.write(ofile, "  nat = 5,"+"\n")
            file.write(ofile, "  ntyp = 3,"+"\n")
            file.write(ofile, "  occupations = 'smearing',"+"\n")
            file.write(ofile, "  smearing = 'gauss',"+"\n")
            file.write(ofile, "/"+"\n")
            file.write(ofile, "&ELECTRONS"+"\n")
            file.write(ofile, "  conv_thr = 1e-08,"+"\n")
            file.write(ofile, "  electron_maxstep = 500,"+"\n")
            file.write(ofile, "/"+"\n")
            file.write(ofile, "&IONS"+"\n")
            file.write(ofile, "  bfgs_ndim = 1,"+"\n")
            file.write(ofile, "/"+"\n")
            file.write(ofile, "&CELL"+"\n")
            file.write(ofile, "  cell_dofree = 'all',"+"\n")
            file.write(ofile, "  press_conv_thr = 0.1D0,"+"\n")
            file.write(ofile, "/"+"\n")
            file.write(ofile, "ATOMIC_SPECIES"+"\n")
            file.write(ofile, A_ele_list[i] + '  ' + mass_dict[A_ele_list[i]] + '  ' + pp_dict[A_ele_list[i]] + '\n')
            file.write(ofile, B_ele_list[j] + '  ' + mass_dict[B_ele_list[j]] + '  ' + pp_dict[B_ele_list[j]] + '\n')
            file.write(ofile, X_ele_list[k] + '  ' + mass_dict[X_ele_list[k]] + '  ' + pp_dict[X_ele_list[k]] + '\n')
            file.write(ofile, "ATOMIC_POSITIONS crystal"+"\n")
            file.write(ofile, A_ele_list[i] + '  ' + str(x[0]) + '  ' + str(y[0]) + '  ' + str(z[0]) + '  ' + '\n')
            file.write(ofile, B_ele_list[j] + '  ' + str(x[1]) + '  ' + str(y[1]) + '  ' + str(z[1]) + '  ' + '\n')
            file.write(ofile, X_ele_list[k] + '  ' + str(x[2]) + '  ' + str(y[2]) + '  ' + str(z[2]) + '  ' + '\n')
            file.write(ofile, X_ele_list[k] + '  ' + str(x[3]) + '  ' + str(y[3]) + '  ' + str(z[3]) + '  ' + '\n')
            file.write(ofile, X_ele_list[k] + '  ' + str(x[4]) + '  ' + str(y[4]) + '  ' + str(z[4]) + '  ' + '\n')
            file.write(ofile, "K_POINTS automatic"+"\n")
            file.write(ofile, "  12 12 12 0 0 0"+"\n")
            file.write(ofile, "CELL_PARAMETERS bohr"+"\n")
            file.write(ofile, str(cellparamx[0]) + '  ' + str(cellparamy[0]) + '  ' + str(cellparamz[0]) + '\n')
            file.write(ofile, str(cellparamx[1]) + '  ' + str(cellparamy[1]) + '  ' + str(cellparamz[1]) + '\n')
            file.write(ofile, str(cellparamx[2]) + '  ' + str(cellparamy[2]) + '  ' + str(cellparamz[2]) + '\n')
            file.close(ofile)
"""for line in ifile:
	items = str.split(line)
	atom = items[1]
	mass[k] = float(mass_dict[atom])
#	print float(val_dict[atom])
	x[k] = float(items[2])
	y[k] = float(items[3])
	z[k] = float(items[4])
	charge[k] = (float(val_dict[atom]) - float(items[5]))
#	print charge[k], mass[k]
	k = k + 1
file.close(ifile)
print numpy.sum(charge)
#print len(mass_dict)
#print len(val_dict)
xc = 0.0
yc = 0.0
zc = 0.0
for i in range(n):
	xc = xc + (mass[i]*x[i])
	yc = yc + (mass[i]*y[i])
	zc = zc + (mass[i]*z[i])
xc = xc/numpy.sum(mass) 
yc = yc/numpy.sum(mass)
zc = zc/numpy.sum(mass)
#print xc,yc,zc
dipolex = 0.0
dipoley = 0.0
dipolez = 0.0
for i in range(n):
	dipolex += (charge[i]*x[i])
	dipoley += (charge[i]*y[i])
	dipolez += (charge[i]*z[i])
dipolex = dipolex - (xc*numpy.sum(charge))
dipoley = dipoley - (yc*numpy.sum(charge))
dipolez = dipolez - (zc*numpy.sum(charge))
print ("Dipole moment in X is %.6f" %dipolex) 
print ("Dipole moment in Y is %.6f" %dipoley)
print ("Dipole moment in Z is %.6f" %dipolez)
dipole = numpy.array([dipolex,dipoley,dipolez])
print ("Dipole moment magnitude is %.6f" %(numpy.linalg.norm(dipole)))
print dipole"""
