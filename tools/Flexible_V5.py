'''
Translation and Rotation of EFP parameters using 3X3 Rotation Matrix.
(Zero Order Approximation)
=======================================================================
Calculate rotation matrix between conformation A and B, 
in Bohr (.efp) and Angs unit (.inp) using three atoms. 
The order of the atoms "must" be the same for both conformations.
Usage: python Flexible.py known.efp unknown.inp --> unknown.efp
'''
import argparse,sys,os,shutil,re,math,fnmatch
import numpy as np
from numpy import *
def add_hydrogen(back_xyz,side_xyz):
	convfactor    = 1/0.529177249 
	x1,y1,z1      = back_xyz[0],back_xyz[1],back_xyz[2]
	x2,y2,z2      = side_xyz[0],side_xyz[1],side_xyz[2]
	actual_length = math.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))
	desire_length = 1.07886*convfactor

	x3 = ((x2-x1)*desire_length/actual_length)+x1
	y3 = ((y2-y1)*desire_length/actual_length)+y1
	z3 = ((z2-z1)*desire_length/actual_length)+z1

	return np.array([x3,y3,z3])

def read_efp_param(reflines1):
	parm_file,find_file = [],[]
	for obj in reflines1:
		parm_file.append(obj)
		find_file.append(obj.strip())
#	print parm_file.index(" COORDINATES (BOHR)\n")
	return parm_file,find_file

def read_new_struc(reflines2):
	convfactor  = 1/0.529177249 
#	convfactor  = 1
	temp1,temp2 = [],[]
	lineN       = 0
	for obj in reflines2:
		line  = obj.strip()
		lineN += 1
		if '$data' in line:
			coord1 = lineN+2
		elif "H000" in line:
			coord2 = lineN+1
#		elif "$end" in line:
#			coord2 = lineN
		else:
			continue
	lineN = 0
	for obj in reflines2:
		line  = obj.strip()
		lineN += 1
		if lineN > coord1 and lineN < coord2:
			temp2.append(float(line.split()[2])*convfactor)
			temp2.append(float(line.split()[3])*convfactor)
			temp2.append(float(line.split()[4])*convfactor)
			temp1.append(temp2)
			temp2 = []
	return np.array(temp1)

def save_efp_coord_and_ele(parm_file,find_file,save_pol,save_disp,save_xr):
#       =======================================================================
#       Save reference coordinates and bond-mid points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	for i in range(find_file.index("COORDINATES (BOHR)")+1,find_file.index("MONOPOLES")-1,1):
		if find_file[i].split()[0][0:2] != "BO":
			temp2.append(find_file[i].split()[0])
			temp2.append(float(find_file[i].split()[1]))
			temp2.append(float(find_file[i].split()[2]))
			temp2.append(float(find_file[i].split()[3]))
			temp1.append(temp2)
			temp2 = []
		else:
			temp4.append(find_file[i].split()[0])
			temp4.append(float(find_file[i].split()[1]))
			temp4.append(float(find_file[i].split()[2]))
			temp4.append(float(find_file[i].split()[3]))
			temp3.append(temp4)
			temp4 = []
	coor_xyz1 = np.array(temp1)
	bmid_xyz1 = np.array(temp3)
#       =======================================================================
#       Save reference monopole points
#       =======================================================================
	ref_monop = []
	for i in range(find_file.index("MONOPOLES"),find_file.index("DIPOLES"),1):
		ref_monop.append(parm_file[i])
#       =======================================================================
#       Save reference dipole points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	for i in range(find_file.index("DIPOLES")+1,find_file.index("QUADRUPOLES")-1,1):
		if find_file[i].split()[0][0:2] != "BO":
			temp2.append(find_file[i].split()[0])
			temp2.append(float(find_file[i].split()[1]))
			temp2.append(float(find_file[i].split()[2]))
			temp2.append(float(find_file[i].split()[3]))
			temp1.append(temp2)
			temp2 = []
		else:
			temp4.append(find_file[i].split()[0])
			temp4.append(float(find_file[i].split()[1]))
			temp4.append(float(find_file[i].split()[2]))
			temp4.append(float(find_file[i].split()[3]))
			temp3.append(temp4)
			temp4 = []
	ref_coor_dip = np.array(temp1)
	ref_bmid_dip = np.array(temp3)
#       =======================================================================
#       Save reference quadrupole points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	for i in range(find_file.index("QUADRUPOLES")+1,find_file.index("OCTUPOLES")-1,1):
		if len(find_file[i].split()) == 6:
			if find_file[i].split()[0][0:2] != "BO":
				temp2.append(find_file[i].split()[0])
				temp2.append(float(find_file[i].split()[1])) # XX
				temp2.append(float(find_file[i].split()[2])) # YY
				temp2.append(float(find_file[i].split()[3])) # ZZ
				temp2.append(float(find_file[i].split()[4])) # XY
			else:
				temp4.append(find_file[i].split()[0])
				temp4.append(float(find_file[i].split()[1])) # XX
				temp4.append(float(find_file[i].split()[2])) # YY 
				temp4.append(float(find_file[i].split()[3])) # ZZ
				temp4.append(float(find_file[i].split()[4])) # XY
		else:
			if len(temp2) > 0:
				temp2.append(float(find_file[i].split()[0])) # XZ
				temp2.append(float(find_file[i].split()[1])) # YZ
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				temp4.append(float(find_file[i].split()[0])) # XZ
				temp4.append(float(find_file[i].split()[1])) # YZ
				temp3.append(temp4)
				temp4 = []
			else:
				print("WHAT? check quadrupole")
				exit()
	ref_coor_qup = np.array(temp1)
	ref_bmid_qup = np.array(temp3)
#       =======================================================================
#       Save reference octupole points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]


	if not save_pol:
		end_of_ocpol = find_file.index("POLARIZABLE POINTS")-1
	else:
		if not save_disp:
			end_of_ocpol = find_file.index("DYNAMIC POLARIZABLE POINTS")-1
		else:
			if not save_xr:
				end_of_ocpol = find_file.index("PROJECTION BASIS SET")-1
			else:
				end_of_ocpol = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1

#	if not save_pol:
#		end_of_ocpol = find_file.index("POLARIZABLE POINTS")-1
#	else:
		end_of_ocpol = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1
	for i in range(find_file.index("OCTUPOLES")+1,end_of_ocpol,1):
		if len(find_file[i].split()) == 6:
			if find_file[i].split()[0][0:2] != "BO":
				temp2.append(find_file[i].split()[0])
				temp2.append(float(find_file[i].split()[1])) # XXX
				temp2.append(float(find_file[i].split()[2])) # YYY
				temp2.append(float(find_file[i].split()[3])) # ZZZ
				temp2.append(float(find_file[i].split()[4])) # XXY
			else:
				temp4.append(find_file[i].split()[0])
				temp4.append(float(find_file[i].split()[1])) # XXX
				temp4.append(float(find_file[i].split()[2])) # YYY
				temp4.append(float(find_file[i].split()[3])) # ZZZ
				temp4.append(float(find_file[i].split()[4])) # XXY
		elif len(find_file[i].split()) == 5:
			if len(temp2) > 0:
				temp2.append(float(find_file[i].split()[0])) # XXZ
				temp2.append(float(find_file[i].split()[1])) # XYY
				temp2.append(float(find_file[i].split()[2])) # YYZ
				temp2.append(float(find_file[i].split()[3])) # XZZ
			elif len(temp4) > 0:
				temp4.append(float(find_file[i].split()[0])) # XXZ
				temp4.append(float(find_file[i].split()[1])) # XYY
				temp4.append(float(find_file[i].split()[2])) # YYZ
				temp4.append(float(find_file[i].split()[3])) # XZZ
			else:
				print("WHAT? check octupole")
				exit()
		elif len(find_file[i].split()) == 2:
			if len(temp2) > 0:
				temp2.append(float(find_file[i].split()[0])) # YZZ
				temp2.append(float(find_file[i].split()[1])) # XYZ
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				temp4.append(float(find_file[i].split()[0])) # YZZ
				temp4.append(float(find_file[i].split()[1])) # XYZ
				temp3.append(temp4)
				temp4 = []
			else:
				print("WHAT? check octupole")
				exit()
		else:
			print("WHAT? check octupole")
			exit()
	ref_coor_ocp = np.array(temp1)
	ref_bmid_ocp = np.array(temp3)

	screen = []
	for i in range(find_file.index("SCREEN2      (FROM VDWSCL=   0.700)"),find_file.index("$END")+1,1):
		screen.append(parm_file[i])

	return [coor_xyz1,bmid_xyz1,ref_monop,ref_coor_dip,ref_bmid_dip,ref_coor_qup,ref_bmid_qup,ref_coor_ocp,ref_bmid_ocp,screen]

def save_efp_pol(parm_file,find_file,save_pol,save_disp,save_xr):
#       =======================================================================
#       Save reference static polarizable points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	if not save_disp:
		end_of_pol = find_file.index("DYNAMIC POLARIZABLE POINTS")-1
	else:
		if not save_xr:
			end_of_pol = find_file.index("PROJECTION BASIS SET")-1
		else:
			end_of_pol = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1

	for i in range(find_file.index("POLARIZABLE POINTS")+1,end_of_pol,1):
		if find_file[i].split()[0][0:2] == "CT":
			temp2.append(find_file[i].split()[0])
			temp2.append(float(find_file[i].split()[1])) # X
			temp2.append(float(find_file[i].split()[2])) # Y
			temp2.append(float(find_file[i].split()[3])) # Z
			temp1.append(temp2)
			temp2 = []
		else:
			if len(find_file[i].split()) == 5:
				temp4.append(float(find_file[i].split()[0])) # XX XZ
				temp4.append(float(find_file[i].split()[1])) # YY YZ
				temp4.append(float(find_file[i].split()[2])) # ZZ YX
				temp4.append(float(find_file[i].split()[3])) # XY ZX
			elif len(find_file[i].split()) == 1:
				temp4.append(float(find_file[i].split()[0])) # ZY
				temp3.append(temp4)
				temp4 = []
			else:
				print("WHAT? check static polarizability")
				exit()
	ref_lmo = np.array(temp1)
	ref_pol = np.array(temp3)

	return ref_lmo,ref_pol

def save_efp_dyn(parm_file,find_file,save_pol,save_disp,save_xr):
#       =======================================================================
#       Save reference dynamic polarizable points
#       =======================================================================
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	if not save_xr:
		end_of_dyn = find_file.index("PROJECTION BASIS SET")-1
	else:
		end_of_dyn = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1
	for i in range(find_file.index("DYNAMIC POLARIZABLE POINTS")+1,end_of_dyn,1):
		if find_file[i].split()[0] == "CT":
			temp2.append(find_file[i].split()[0]+find_file[i].split()[1])
			temp2.append(float(find_file[i].split()[2]))
			temp2.append(float(find_file[i].split()[3]))
			temp2.append(float(find_file[i].split()[4]))
			temp1.append(temp2)
			temp2 = []
		else:
			if len(find_file[i].split()) == 5:
				temp4.append(float(find_file[i].split()[0])) # XX XZ
				temp4.append(float(find_file[i].split()[1])) # YY YZ
				temp4.append(float(find_file[i].split()[2])) # ZZ YX
				temp4.append(float(find_file[i].split()[3])) # XY ZX
			elif len(find_file[i].split()) == 1:
				temp4.append(float(find_file[i].split()[0])) # ZY
				temp3.append(temp4)
				temp4 = []
			else:
				print("WHAT? check dispersion")
				exit()
	ref_dyn = np.array(temp3)
	return temp1,ref_dyn

def save_efp_exr(parm_file,find_file,save_pol,save_disp,save_xr):
#       =======================================================================
#       Save reference exchange-repulsion parameters 
#       =======================================================================
	prj_basis       = [] # need to change multiplicity if it is more than two
	orbital,n_basis = [],[]
	for i in range(find_file.index("PROJECTION BASIS SET"),find_file.index("MULTIPLICITY    1")+3,1):
		prj_basis.append(parm_file[i])
		line = find_file[i].split()
		if len(line) > 0:
			if line[0] == "S":
				orbital.append(line[0])
				n_basis.append(1)
			elif line[0] == "P":
				orbital.append(line[0])
				n_basis.append(3)
			elif line[0] == "L":
				orbital.append(line[0])
				n_basis.append(4)
			elif line[0] == "D":
				orbital.append(line[0])
				n_basis.append(6)
			elif line[0] == "F":
				orbital.append(line[0])
				n_basis.append(10)
			elif line[0] == "G":
				print("G orbital rotation is not implemented!!")
				print("Quit")
				quit()
			else:
				continue
		else:
			continue
	temp1 = []
	for i in range(find_file.index("MULTIPLICITY    1")+3,find_file.index("FOCK MATRIX ELEMENTS"),1):
		temp1.append(parm_file[i][5:20])
		temp1.append(parm_file[i][20:35])
		temp1.append(parm_file[i][35:50])
		temp1.append(parm_file[i][50:65])
		temp1.append(parm_file[i][65:80])
	fock_mat = []
	for i in range(find_file.index("FOCK MATRIX ELEMENTS"),find_file.index("LMO CENTROIDS"),1):
		fock_mat.append(parm_file[i])
	temp2,temp3 = [],[]
	for i in range(find_file.index("LMO CENTROIDS")+1,find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1,1):
		temp3.append(find_file[i].split()[0])
		temp3.append(float(find_file[i].split()[1]))
		temp3.append(float(find_file[i].split()[2]))
		temp3.append(float(find_file[i].split()[3]))
		temp2.append(temp3)
		temp3 = []
	xr_lmo = np.array(temp2)

	temp4 = []
	for i in range(len(temp1)):
		if len(temp1[i]) < 10:
			continue
		else:
			temp4.append(float(temp1[i]))

	n_lmo,n_bas = int(find_file[find_file.index("MULTIPLICITY    1")+2].split()[2]),int(find_file[find_file.index("MULTIPLICITY    1")+2].split()[3])
	ref_wvf     = np.zeros((n_lmo,n_bas))
	lineN       = 0
	for i in range(len(ref_wvf)):
		for j in range(len(ref_wvf[i])):
			ref_wvf[i][j] = temp4[lineN]
			lineN        += 1

	return orbital,n_basis,ref_wvf,fock_mat,xr_lmo

def organize_connectivity(C):
	temp = []
	for i in range(0,len(C)-1,1):
		v = C[i][0]
		for j in range(i+1,len(C),1):
			for k in C[i]:
				s = C[i]+C[j]
				s = remove_duplicates(s)
				if s[0] == 0 or s[1] == 0 or s[2] == 0:
					continue
				else:
					temp.append(s)
	return temp

def remove_duplicates(li):
	j,temp = 0,[0,0,0]

	if len(list(set(li))) == 3:
		for i in li:
			if li.count(i) == 2:
				temp[1] = i
			else:
				temp[j] = i
				j += 2
		return temp
	else:
		return temp

def reverse_lists(li):
	temp = list()
	for i in range(len(li)):
		a = li[i]
		a.reverse()
		temp.append(a)
	return temp

def finalize_connectivity(li):
	temp = list()
	for i in li:
		if i not in temp:
			temp.append(i)
	return temp

def define_topology(midorder,nearest): # [atom1,atom2,atom3,midorder1-2,midorder2-3] atom2 is the center
	for i in range(len(nearest)):
		k,temp = 0,[]
		if nearest[i][0]-nearest[i][1] > 0:
			m1,m2 = nearest[i][0],nearest[i][1]
		else:
			m1,m2 = nearest[i][1],nearest[i][0]
		if nearest[i][2]-nearest[i][1] > 0:
			m3,m4 = nearest[i][2],nearest[i][1]
		else:
			m3,m4 = nearest[i][1],nearest[i][2]
		temp.append(str(m1)+str(m2))
		temp.append(str(m3)+str(m4))

		for j in range(len(temp)):
			for k in range(len(midorder)):
				if str(temp[j]) == str(midorder[k][0])+str(midorder[k][1]):
					nearest[i] += [k+1]
				else:
					continue
	return nearest

def save_bond_connection(bmid,coord_xyz2):
	convfactor        = 1/0.529177249 
	temp1,temp2,temp3 = [],[],[]
	for i in range(len(bmid)):
		if len(bmid[i][0].split("BO")[1]) == 2:
			atom1,atom2 = bmid[i][0].split("BO")[1][0:1],bmid[i][0].split("BO")[1][1:2]
		elif len(bmid[i][0].split("BO")[1]) == 3:
			atom1,atom2 = bmid[i][0].split("BO")[1][0:2],bmid[i][0].split("BO")[1][2:3]
		else:
			if int(bmid[i][0].split("BO")[1][0:2]) > int(bmid[i][0].split("BO")[1][2:4]):
				atom1,atom2 = bmid[i][0].split("BO")[1][0:2],bmid[i][0].split("BO")[1][2:4]
			else:
				atom1,atom2 = bmid[i][0].split("BO")[1][0:3],bmid[i][0].split("BO")[1][3:4]

		temp2.append(bmid[i][0])
		temp2.append((coord_xyz2[int(atom1)-1][0]+coord_xyz2[int(atom2)-1][0])/2)
		temp2.append((coord_xyz2[int(atom1)-1][1]+coord_xyz2[int(atom2)-1][1])/2)
		temp2.append((coord_xyz2[int(atom1)-1][2]+coord_xyz2[int(atom2)-1][2])/2)
		temp1.append(temp2)
		temp2 = []
		temp3.append([int(atom1),int(atom2)])

	bmid_xyz2   = np.array(temp1)
	connection =  finalize_connectivity(organize_connectivity(temp3)+organize_connectivity(reverse_lists(temp3)))

	return bmid_xyz2,define_topology(reverse_lists(temp3),connection)

def figure_atoms(atom,r,j):
	if "O" in atom and r < 0.8:
		return j
	elif "N" in atom and r < 0.8:
		return j
	elif "S" in atom and r < 1.15:
		return j
	else:
		return "Nolone"

def lmo_position(coor,bmid,ref_lmo,anum,connection):
	outputs = []
	for i in range(len(connection)):
		temp1,temp2,temp3,temp4,temp5 = [],[],[],[],[]
		for j in range(len(ref_lmo)):
			x1,y1,z1    = ref_lmo[j][1],ref_lmo[j][2],ref_lmo[j][3]
			x2,y2,z2    = bmid[connection[i][3]-1][1],bmid[connection[i][3]-1][2],bmid[connection[i][3]-1][3]
			x3,y3,z3    = bmid[connection[i][4]-1][1],bmid[connection[i][4]-1][2],bmid[connection[i][4]-1][3]
			dx1,dy1,dz1 = float(x1)-float(x2),float(y1)-float(y2),float(z1)-float(z2)
			dx2,dy2,dz2 = float(x1)-float(x3),float(y1)-float(y3),float(z1)-float(z3)
			r1          = dx1**2+dy1**2+dz1**2
			r2          = dx2**2+dy2**2+dz2**2
			temp2.append(r1)
			temp4.append(r2)
		L1        = temp2.index(min(temp2))
		R1        = temp4.index(min(temp4))
		temp1.append(L1) # left
		temp3.append(R1) # right
		r1_min    = temp2[L1]
		r2_min    = temp4[R1]
		temp2[L1] = 1000
		temp4[R1] = 1000
		if math.sqrt(r1_min)/math.sqrt(min(temp2)) >= 0.6:
			L2 = temp2.index(min(temp2))
			temp1.append(L2)
		if math.sqrt(r2_min)/math.sqrt(min(temp4)) >= 0.6:
			R2 = temp4.index(min(temp4))
			temp3.append(R2)
		
		for j in range(len(ref_lmo)):
			x1,y1,z1    = ref_lmo[j][1],ref_lmo[j][2],ref_lmo[j][3]
			x2,y2,z2    = coor[connection[i][0]-1][1],coor[connection[i][0]-1][2],coor[connection[i][0]-1][3]
			x3,y3,z3    = coor[connection[i][1]-1][1],coor[connection[i][1]-1][2],coor[connection[i][1]-1][3]
			x4,y4,z4    = coor[connection[i][2]-1][1],coor[connection[i][2]-1][2],coor[connection[i][2]-1][3]
			dx1,dy1,dz1 = float(x1)-float(x2),float(y1)-float(y2),float(z1)-float(z2)
			dx2,dy2,dz2 = float(x1)-float(x3),float(y1)-float(y3),float(z1)-float(z3)
			dx3,dy3,dz3 = float(x1)-float(x4),float(y1)-float(y4),float(z1)-float(z4)
			r1          = math.sqrt(dx1**2+dy1**2+dz1**2)
			r2          = math.sqrt(dx2**2+dy2**2+dz2**2)
			r3          = math.sqrt(dx3**2+dy3**2+dz3**2)
			lone1       = figure_atoms(coor[connection[i][0]-1][0],r1,j)
			lone2       = figure_atoms(coor[connection[i][1]-1][0],r2,j)
			lone3       = figure_atoms(coor[connection[i][2]-1][0],r3,j)
			if lone1 != "Nolone":
				temp1.append(lone1)
			if lone2 != "Nolone": 
				temp5.append(lone2)
			if lone3 != "Nolone":
				temp3.append(lone3)
		left  = list(set(temp1))
		cent  = list(set(temp5))
		right = list(set(temp3))

		outputs.append([left,cent,right])

	return outputs

def define_three_points(info,xyz,a,b,c,d,e,f):
	x1,y1,z1 = float(xyz[info[a]-1][d]),float(xyz[info[a]-1][e]),float(xyz[info[a]-1][f])
	x2,y2,z2 = float(xyz[info[b]-1][d]),float(xyz[info[b]-1][e]),float(xyz[info[b]-1][f])
	x3,y3,z3 = float(xyz[info[c]-1][d]),float(xyz[info[c]-1][e]),float(xyz[info[c]-1][f])
	return np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])

def centroid(X):
	C1 = (X[0][0]+X[1][0]+X[2][0])/3
	C2 = (X[0][1]+X[1][1]+X[2][1])/3
	C3 = (X[0][2]+X[1][2]+X[2][2])/3
	return [C1,C2,C3]
#	C = sum(X)/len(X)
#	return C

def kabsch(P,Q):
	C = np.dot(np.transpose(P), Q)
	V, S, W = np.linalg.svd(C)
	d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
	
	if d:
		S[-1] = -S[-1]
		V[:, -1] = -V[:, -1]
	
	# Create Rotation matrix U
	U = np.dot(V, W)
	
	return U

def special_three_points(any_crd):
#	Just for center atom!!!
#	[atom1+atom2/2,center atom,normal_vector]
	a,b,c = any_crd[0],any_crd[1],any_crd[2]
	t = np.empty((3,3))
	for i in range(3):
		t[0][i] = 0.0
		t[1][i] = (a[i]+b[i])/2
	t[2][0] = any_crd[0][1]*any_crd[1][2] - any_crd[0][2]*any_crd[1][1]
	t[2][1] = any_crd[0][2]*any_crd[1][0] - any_crd[0][0]*any_crd[1][2]
	t[2][2] = any_crd[0][0]*any_crd[1][1] - any_crd[0][1]*any_crd[1][0]
	return t

def compute_rotation_matrix(known,unknown,string):
	if string:
#		rot = kabsch(known,unknown)
		t1 = matrix_normalization(special_three_points(matrix_normalization(unknown)))
		t2 = matrix_normalization(special_three_points(matrix_normalization(known)))
		r1  = matrix_finalization(t1)
		r2  = matrix_finalization(t2)
		rot = np.matmul(r1.T,r2)
	else:
		t1  = matrix_normalization(unknown)
		t2  = matrix_normalization(known)
		r1  = matrix_finalization(t1)
		r2  = matrix_finalization(t2)
		rot = np.matmul(r1.T,r2)

	return rot

def matrix_normalization(any_crd):
	a,b,c         = any_crd[0],any_crd[1],any_crd[2]
	t             = np.empty((3,3))
	t1norm,t2norm = 0.0,0.0
	for i in range(3):
		t[0][i] = a[i]-b[i]
		t[1][i] = c[i]-b[i]
		t1norm += pow(t[0][i],2)
		t2norm += pow(t[1][i],2)
	t1norm,t2norm = 1/math.sqrt(t1norm),1/math.sqrt(t2norm)
	for i in range(3):
		t[0][i] = t[0][i]*t1norm
		t[1][i] = t[1][i]*t2norm
	return t

def matrix_finalization(t):
	dotprd  = t[0][0]*t[1][0]+t[0][1]*t[1][1]+t[0][2]*t[1][2]
	t[1][0] = t[1][0]-dotprd*t[0][0]
	t[1][1] = t[1][1]-dotprd*t[0][1]
	t[1][2] = t[1][2]-dotprd*t[0][2]
	t[2][0] = t[0][1]*t[1][2]-t[0][2]*t[1][1]
	t[2][1] = t[0][2]*t[1][0]-t[0][0]*t[1][2]
	t[2][2] = t[0][0]*t[1][1]-t[0][1]*t[1][0]
	for i in range(3):
		tnorm = 0.0
		for j in range(3):
			tnorm += pow(t[i][j],2)
		if tnorm == float(0.0):
			tnorm = 1.0
		else:
			tnorm = tnorm
		for j in range(3):
			t[i][j] = t[i][j]/math.sqrt(tnorm)
	return t

def compute_multipole(tensor,rot):
	temp = np.empty((3,3))
	temp[0][0] = tensor[0]*rot[0][0] + tensor[1]*rot[0][1] + tensor[3]*rot[0][2]    # XX*Rot[0][0] + XY*Rot[1][0] + YZ*Rot[2][0]
	temp[0][1] = tensor[1]*rot[0][0] + tensor[2]*rot[0][1] + tensor[4]*rot[0][2]    # XY*Rot[0][0] + YY*Rot[1][0] + XZ*Rot[2][0]
	temp[0][2] = tensor[3]*rot[0][0] + tensor[4]*rot[0][1] + tensor[5]*rot[0][2]    # YZ*Rot[0][0] + XZ*Rot[1][0] + ZZ*Rot[2][0]
	temp[1][0] = tensor[0]*rot[1][0] + tensor[1]*rot[1][1] + tensor[3]*rot[1][2]    # For OCTUPOLE, LOOK AT THE GAMESS CODE
	temp[1][1] = tensor[1]*rot[1][0] + tensor[2]*rot[1][1] + tensor[4]*rot[1][2]
	temp[1][2] = tensor[3]*rot[1][0] + tensor[4]*rot[1][1] + tensor[5]*rot[1][2]
	temp[2][0] = tensor[0]*rot[2][0] + tensor[1]*rot[2][1] + tensor[3]*rot[2][2]
	temp[2][1] = tensor[1]*rot[2][0] + tensor[2]*rot[2][1] + tensor[4]*rot[2][2]
	temp[2][2] = tensor[3]*rot[2][0] + tensor[4]*rot[2][1] + tensor[5]*rot[2][2]
	return temp

def rotate_dipole(nefdip,rot,points,m):
	efdip         = np.empty((1,3))
	efdip[0]      = np.array([float(points[1]),float(points[2]),float(points[3])])
	efdip[0]      = np.matmul(rot,efdip.T).T
	nefdip[m][0] += efdip[0][0] #x
	nefdip[m][1] += efdip[0][1] #y
	nefdip[m][2] += efdip[0][2] #z
	nefdip[m][3] += 1.0

	return nefdip

def rotate_qupole(nefqua,rot,points,m):
	temp,efqua = np.empty((6,1)),np.empty((3,3))
	temp       = np.array([[float(points[1])],[float(points[4])],[float(points[2])],
			[float(points[5])],[float(points[6])],[float(points[3])]])
	wrk        = compute_multipole(temp,rot)
	for n in range(3):
		efqua[n] = np.matmul(rot,wrk[n])
	nefqua[m][0] += efqua[0][0] #xx
	nefqua[m][1] += efqua[1][1] #yy
	nefqua[m][2] += efqua[2][2] #zz
	nefqua[m][3] += efqua[1][0] #xy
	nefqua[m][4] += efqua[2][0] #xz
	nefqua[m][5] += efqua[2][1] #yz
	nefqua[m][6] += 1.0

	return nefqua

def rotate_ocpole(nefoct,rot,points,m):
	temp01,temp02,temp03 = np.empty((6,1)),np.empty((6,1)),np.empty((6,1))
	temp04,temp05,temp06 = np.empty((3,3)),np.empty((3,3)),np.empty((3,3))
	temp07,temp08,temp09 = np.empty((2,3)),np.empty((2,3)),np.empty((2,3))
	temp10,temp11,efoct  = np.empty((3,3)),np.empty((3,3)),np.empty((6,3))

	temp01 = np.array([[float(points[1])],[float(points[4])],[float(points[6])],
			[float(points[5])],[float(points[10])],[float(points[8])]])
	temp02 = np.array([[float(points[4])],[float(points[6])],[float(points[2])],
			[float(points[10])],[float(points[7])],[float(points[9])]])
	temp03 = np.array([[float(points[5])],[float(points[10])],[float(points[7])],
			[float(points[8])],[float(points[9])],[float(points[3])]])

	wrk1   = compute_multipole(temp01,rot)
	wrk2   = compute_multipole(temp02,rot)
	wrk3   = compute_multipole(temp03,rot)

	for n in range(3):
		temp04[n] = np.matmul(rot,wrk1[n])
		temp05[n] = np.matmul(rot,wrk2[n])
		temp06[n] = np.matmul(rot,wrk3[n])
	temp07[0][0],temp08[0][0],temp09[0][0] = temp04[0][0],temp05[0][0],temp06[0][0]
	temp07[0][1],temp08[0][1],temp09[0][1] = temp04[1][0],temp05[1][0],temp06[1][0]
	temp07[0][2],temp08[0][2],temp09[0][2] = temp04[1][1],temp05[1][1],temp06[1][1]
	temp07[1][0],temp08[1][0],temp09[1][0] = temp04[2][0],temp05[2][0],temp06[2][0]
	temp07[1][1],temp08[1][1],temp09[1][1] = temp04[2][1],temp05[2][1],temp06[2][1]
	temp07[1][2],temp08[1][2],temp09[1][2] = temp04[2][2],temp05[2][2],temp06[2][2]

	for p in range(1):
		for q in range(3):
			temp10[p][q]   = temp07[p][q]
			temp10[p+1][q] = temp08[p][q]
			temp10[p+2][q] = temp09[p][q]
			temp11[p][q]   = temp07[p+1][q]
			temp11[p+1][q] = temp08[p+1][q]
			temp11[p+2][q] = temp09[p+1][q]
	temp07 = np.matmul(temp10.T,rot.T).T
	temp08 = np.matmul(temp11.T,rot.T).T

	v = 0
	for p in range(3):
		efoct[v]   = temp07[p]
		efoct[v+1] = temp08[p]
		v += 2

	nefoct[m][0]  += efoct[0][0] #xxx
	nefoct[m][1]  += efoct[2][2] #yyy
	nefoct[m][2]  += efoct[5][2] #zzz
	nefoct[m][3]  += efoct[0][1] #xxy
	nefoct[m][4]  += efoct[1][0] #xxz
	nefoct[m][5]  += efoct[0][2] #xyy
	nefoct[m][6]  += efoct[3][1] #yyz
	nefoct[m][7]  += efoct[1][2] #xzz
	nefoct[m][8]  += efoct[3][2] #yzz
	nefoct[m][9]  += efoct[1][1] #xyz
	nefoct[m][10] += 1.0

	return nefoct

def prepare_lmo(xyz,vector):
	position = np.array([[float(xyz[1]),float(xyz[2]),float(xyz[3])]])
	tensor   = np.array([[float(vector[0]),float(vector[3]),float(vector[4])],
				[float(vector[6]),float(vector[1]),float(vector[5])],
				[float(vector[7]),float(vector[8]),float(vector[2])]])
	return position,tensor

def prepare_cmo(xyz):
	position = np.array([[float(xyz[1]),float(xyz[2]),float(xyz[3])]])
	return position

def rotate_lmo_tensor(neflmo,nefpol,rot,xyz1,xyz2,points,ref_lmo,ref_pol):
	for pt in points:
		lmo,tensor = prepare_lmo(ref_lmo[pt],ref_pol[pt])
		temp       = [np.empty((1,3))]

		for n in range(3):
			temp[0][0][n] = lmo[0][n]-xyz1[n]
		efp   = np.matmul(rot,temp[0].T).T
		efpol = np.matmul(rot,tensor)
		efpol = np.matmul(efpol,rot.T)

		neflmo[pt][0] += xyz2[0]+efp[0][0] #x
		neflmo[pt][1] += xyz2[1]+efp[0][1] #y
		neflmo[pt][2] += xyz2[2]+efp[0][2] #z
		neflmo[pt][3] += 1.0
		nefpol[pt][0] += efpol[0][0] #xx
		nefpol[pt][1] += efpol[1][1] #yy
		nefpol[pt][2] += efpol[2][2] #zz
		nefpol[pt][3] += efpol[0][1] #xy
		nefpol[pt][4] += efpol[0][2] #xz
		nefpol[pt][5] += efpol[1][2] #yz
		nefpol[pt][6] += efpol[1][0] #yx
		nefpol[pt][7] += efpol[2][0] #zx
		nefpol[pt][8] += efpol[2][1] #zy
		nefpol[pt][9] += 1.0

	return neflmo,nefpol

def rotate_dyn_tensor(nefdmo,nefdol,rot,xyz1,xyz2,points,dyn_lmo,ref_dyn):
	pNum = len(dyn_lmo)/12
	for k in range(12):
		rep = pNum*k
		for pt in points:
			lmo,tensor = prepare_lmo(dyn_lmo[pt+rep],ref_dyn[pt+rep])
			temp       = [np.empty((1,3))]

			for n in range(3):
				temp[0][0][n] = lmo[0][n]-xyz1[n]
			efp   = np.matmul(rot,temp[0].T).T
			efpol = np.matmul(rot,tensor)
			efpol = np.matmul(efpol,rot.T)
			nefdmo[pt+rep][0] += xyz2[0]+efp[0][0] #x
			nefdmo[pt+rep][1] += xyz2[1]+efp[0][1] #y
			nefdmo[pt+rep][2] += xyz2[2]+efp[0][2] #z
			nefdmo[pt+rep][3] += 1.0
			nefdol[pt+rep][0] += efpol[0][0] #XX
			nefdol[pt+rep][1] += efpol[1][1] #YY
			nefdol[pt+rep][2] += efpol[2][2] #ZZ
			nefdol[pt+rep][3] += efpol[0][1] #XY
			nefdol[pt+rep][4] += efpol[0][2] #XZ
			nefdol[pt+rep][5] += efpol[1][2] #YZ
			nefdol[pt+rep][6] += efpol[1][0] #YX
			nefdol[pt+rep][7] += efpol[2][0] #ZX
			nefdol[pt+rep][8] += efpol[2][1] #ZY
			nefdol[pt+rep][9] += 1.0

	return nefdmo,nefdol

def rotate_projection_wf(nefcmo,nefcof,rot,xyz1,xyz2,points,xr_lmo,orbital,n_basis,ref_wvf):
	for pt in points:
		tot = 0
		for n in range(len(orbital)):
			if orbital[n] == "S":
				tot += n_basis[n]
				nefcof = s_orbital(nefcof,ref_wvf[pt],pt,tot-1)
			elif orbital[n] == "P":
				tot  += n_basis[n]
				nefcof = pl_orbital(nefcof,rot,ref_wvf[pt],pt,tot-1,True)
			elif orbital[n] == "L":
				tot += n_basis[n]
				nefcof = pl_orbital(nefcof,rot,ref_wvf[pt],pt,tot-1,False)
			elif orbital[n] == "D":
				tot += n_basis[n]
				nefcof = d_orbital(nefcof,rot,ref_wvf[pt],pt,tot-1)
			elif orbital[n] == "F":
				tot += n_basis[n]
				nefcof = f_orbital(nefcof,rot,ref_wvf[pt],pt,tot-1)
			else:
				print("I said G orbital is not implemented!!")
				quit()
		nefcof[pt][tot] += 1.0

		lmo  = prepare_cmo(xr_lmo[pt])
		temp = [np.empty((1,3))]

		for n in range(3):
			temp[0][0][n] = lmo[0][n]-xyz1[n]
		efp   = np.matmul(rot,temp[0].T).T
		nefcmo[pt][0] += xyz2[0]+efp[0][0] #x
		nefcmo[pt][1] += xyz2[1]+efp[0][1] #y
		nefcmo[pt][2] += xyz2[2]+efp[0][2] #z
		nefcmo[pt][3] += 1.0

	return nefcmo,nefcof

def s_orbital(nefcof,rotate,pt,end):
	nefcof[pt][end] += rotate[end]

	return nefcof

def pl_orbital(nefcof,rot,rotate,pt,end,string):
	wf    = np.empty((1,3))
	wf[0] = np.array([rotate[end-2],rotate[end-1],rotate[end]])
	wf[0] = np.matmul(rot,wf[0].T).T

	if string:
		nefcof[pt][end-2] += wf[0][0]
		nefcof[pt][end-1] += wf[0][1]
		nefcof[pt][end]   += wf[0][2]
	else: # L orbital: does not touch the first coefficient
		nefcof[pt][end-3] += rotate[end-3]
		nefcof[pt][end-2] += wf[0][0]
		nefcof[pt][end-1] += wf[0][1]
		nefcof[pt][end]   += wf[0][2]

	return nefcof

def d_orbital(nefcof,rot,rotate,pt,end):
	const   = 1.154700538
	temp,wf = np.empty((6,1)),np.empty((3,3))
	temp    = np.array([[rotate[end-5]],[rotate[end-2]/const],[rotate[end-4]],
			[rotate[end-1]/const],[rotate[end]/const],[rotate[end-3]]])

	wrk = compute_multipole(temp,rot)
	for q in range(3):
		wf[q] = np.matmul(rot,wrk[q])
	nefcof[pt][end-5] += wf[0][0]        #xx 
	nefcof[pt][end-4] += wf[1][1]        #yy
	nefcof[pt][end-3] += wf[2][2]        #zz
	nefcof[pt][end-2] += wf[1][0]*const  #xy
	nefcof[pt][end-1] += wf[2][0]*const  #xz
	nefcof[pt][end-0] += wf[2][1]*const  #yz

	return nefcof

def f_orbital(nefcof,rot,rotate,pt,end):
	two,three,sqrt3,sqrt5 = 2.0,3.0,1.732050808,2.236067978
	temp01,temp02,temp03  = np.empty((6,1)),np.empty((6,1)),np.empty((6,1))
	temp04,temp05,temp06  = np.empty((3,3)),np.empty((3,3)),np.empty((3,3))
	temp07,temp08,temp09  = np.empty((2,3)),np.empty((2,3)),np.empty((2,3))
	temp10,temp11,wf      = np.empty((3,3)),np.empty((3,3)),np.empty((6,3))

	temp01 = np.array([[rotate[end-9]],[rotate[end-6]*sqrt5/three],[rotate[end-4]*sqrt5/three],
			[rotate[end-5]*sqrt5/three],[rotate[end]*(sqrt5/three)*(sqrt3/two)],[rotate[end-2]*sqrt5/three]])
	temp02 = np.array([[rotate[end-6]*sqrt5/three],[rotate[end-4]*sqrt5/three],[rotate[end-8]],
			[rotate[end]*(sqrt5/three)*(sqrt3/two)],[rotate[end-3]*sqrt5/three],[rotate[end-1]*sqrt5/three]])
	temp03 = np.array([[rotate[end-5]*sqrt5/three],[rotate[end]*(sqrt5/three)*(sqrt3/two)],[rotate[end-3]*sqrt5/three],
			[rotate[end-2]*sqrt5/three],[rotate[end-1]*sqrt5/three],[rotate[end-7]]])
	
	wrk1 = compute_multipole(temp01,rot)
	wrk2 = compute_multipole(temp02,rot)
	wrk3 = compute_multipole(temp03,rot)

	for q in range(3):
		temp04[q] = np.matmul(rot,wrk1[q])
		temp05[q] = np.matmul(rot,wrk2[q])
		temp06[q] = np.matmul(rot,wrk3[q])
	temp07[0][0],temp08[0][0],temp09[0][0] = temp04[0][0],temp05[0][0],temp06[0][0]
	temp07[0][1],temp08[0][1],temp09[0][1] = temp04[1][0],temp05[1][0],temp06[1][0]
	temp07[0][2],temp08[0][2],temp09[0][2] = temp04[1][1],temp05[1][1],temp06[1][1]
	temp07[1][0],temp08[1][0],temp09[1][0] = temp04[2][0],temp05[2][0],temp06[2][0]
	temp07[1][1],temp08[1][1],temp09[1][1] = temp04[2][1],temp05[2][1],temp06[2][1]
	temp07[1][2],temp08[1][2],temp09[1][2] = temp04[2][2],temp05[2][2],temp06[2][2]

	for r in range(1):
		for s in range(3):
			temp10[r][s]   = temp07[r][s]
			temp10[r+1][s] = temp08[r][s]
			temp10[r+2][s] = temp09[r][s]
			temp11[r][s]   = temp07[r+1][s]
			temp11[r+1][s] = temp08[r+1][s]
			temp11[r+2][s] = temp09[r+1][s]
	temp07 = np.matmul(temp10.T,rot.T).T
	temp08 = np.matmul(temp11.T,rot.T).T

	v = 0
	for r in range(3):
		wf[v]   = temp07[r]
		wf[v+1] = temp08[r]
		v += 2

	nefcof[pt][end-9] += wf[0][0]                       #xxx
	nefcof[pt][end-8] += wf[2][2]                       #yyy
	nefcof[pt][end-7] += wf[5][2]                       #zzz
	nefcof[pt][end-6] += wf[0][1]/sqrt5*three           #xxy
	nefcof[pt][end-5] += wf[1][0]/sqrt5*three           #xxz
	nefcof[pt][end-4] += wf[0][2]/sqrt5*three           #xyy
	nefcof[pt][end-3] += wf[3][1]/sqrt5*three           #yyz
	nefcof[pt][end-2] += wf[1][2]/sqrt5*three           #xzz
	nefcof[pt][end-1] += wf[3][2]/sqrt5*three           #yzz
	nefcof[pt][end-0] += wf[1][1]/sqrt5*three/sqrt3*two #xyz

	return nefcof

def wrt_ele(find_file,parm_file,xyz,nefdip,nefqua,nefoct,no_pol,no_disp,no_xr,fragname):
	string = ''

	for i in range(find_file.index("COORDINATES (BOHR)")+1):
		if '$' in parm_file[i]:
			string += ' $'+fragname+'\n'
		else:
			string += parm_file[i]

	lineN = 0
	for i in range(find_file.index("COORDINATES (BOHR)")+1,find_file.index("MONOPOLES")-1,1):
		x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
		x1 = '{:{align}{width}}'.format('%.10f'%float(xyz[lineN][0]),align='>',width=23-len(find_file[i].split()[0]))
		x2 = '{:{align}{width}}'.format('%.10f'%float(xyz[lineN][1]),align='>',width=15)
		x3 = '{:{align}{width}}'.format('%.10f'%float(xyz[lineN][2]),align='>',width=15)
		x4 = '{:{align}{width}}'.format('%s'%find_file[i].split()[4],align='>',width=12)
		x5 = '{:{align}{width}}'.format('%s'%find_file[i].split()[5],align='>',width=5)
		lineN += 1
		string += (x0+x1+x2+x3+x4+x5+'\n')
	string += ' STOP\n'

	for i in range(find_file.index("MONOPOLES"),find_file.index("DIPOLES")+1,1):
		string += parm_file[i]

	lineN = 0
	for i in range(find_file.index("DIPOLES")+1,find_file.index("QUADRUPOLES")-1,1):
		v1,v2,v3 = nefdip[lineN][0]/nefdip[lineN][3],nefdip[lineN][1]/nefdip[lineN][3],nefdip[lineN][2]/nefdip[lineN][3]
		x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
		x1 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=24-len(find_file[i].split()[0]))
		x2 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=16)
		x3 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=16)
		lineN += 1
		string += (x0+x1+x2+x3+'\n')
	string += ' STOP\n QUADRUPOLES\n'

	lineN = 0
	y1    = 1
	for i in range(find_file.index("QUADRUPOLES")+1,find_file.index("OCTUPOLES")-1,1):
		if y1%2 == 1:
			v1,v2 = nefqua[lineN][0]/nefqua[lineN][6],nefqua[lineN][1]/nefqua[lineN][6]
			v3,v4 = nefqua[lineN][2]/nefqua[lineN][6],nefqua[lineN][3]/nefqua[lineN][6]
			x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
			x1 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=24-len(find_file[i].split()[0]))
			x2 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=16)
			x4 = '{:{align}{width}}'.format('%.10f'%v4,align='>',width=16)
			string += (x0+x1+x2+x3+x4+' >'+'\n')
		else:
			v5,v6 = nefqua[lineN][4]/nefqua[lineN][6],nefqua[lineN][5]/nefqua[lineN][6]
			x1 = '{:{align}{width}}'.format('%.10f'%v5,align='>',width=24)
			x2 = '{:{align}{width}}'.format('%.10f'%v6,align='>',width=16)
			lineN += 1
			string += (x1+x2+'\n')
		y1 += 1
	string += ' STOP\n OCTUPOLES\n'

	lineN = 0
	y2    = 1
	if not no_pol:
		end_of_ocpol = find_file.index("POLARIZABLE POINTS")-1
	else:
		if not no_disp:
			end_of_ocpol = find_file.index("DYNAMIC POLARIZABLE POINTS")-1
		else:
			if not no_xr:
				end_of_ocpol = find_file.index("PROJECTION BASIS SET")-1
			else:
				end_of_ocpol = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1
	for i in range(find_file.index("OCTUPOLES")+1,end_of_ocpol,1):
		if y2%3 == 1:
			v1,v2 = nefoct[lineN][0]/nefoct[lineN][10],nefoct[lineN][1]/nefoct[lineN][10]
			v3,v4 = nefoct[lineN][2]/nefoct[lineN][10],nefoct[lineN][3]/nefoct[lineN][10]
			x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
			x1 = '{:{align}{width}}'.format('%.9f'%v1,align='>',width=24-len(find_file[i].split()[0]))
			x2 = '{:{align}{width}}'.format('%.9f'%v2,align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.9f'%v3,align='>',width=16)
			x4 = '{:{align}{width}}'.format('%.9f'%v4,align='>',width=16)
			string += (x0+x1+x2+x3+x4+' >'+'\n')
		elif y2%3 == 2:
			v5,v6 = nefoct[lineN][4]/nefoct[lineN][10],nefoct[lineN][5]/nefoct[lineN][10]
			v7,v8 = nefoct[lineN][6]/nefoct[lineN][10],nefoct[lineN][7]/nefoct[lineN][10]
			x1 = '{:{align}{width}}'.format('%.9f'%v5,align='>',width=24)
			x2 = '{:{align}{width}}'.format('%.9f'%v6,align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.9f'%v7,align='>',width=16)
			x4 = '{:{align}{width}}'.format('%.9f'%v8,align='>',width=16)
			string += (x1+x2+x3+x4+' >'+'\n')
		else:
			v9,v10 = nefoct[lineN][8]/nefoct[lineN][10],nefoct[lineN][9]/nefoct[lineN][10]
			x1 = '{:{align}{width}}'.format('%.9f'%v9,align='>',width=24)
			x2 = '{:{align}{width}}'.format('%.9f'%v10,align='>',width=16)
			lineN += 1
			string += (x1+x2+'\n')
		y2 += 1
	string += ' STOP\n'

	return string

def wrt_pol(find_file,neflmo,nefpol,no_disp,no_xr):
	string = ' POLARIZABLE POINTS\n'
	if not no_disp:
		end_of_pol = find_file.index("DYNAMIC POLARIZABLE POINTS")-1
	else:
		if not no_xr:
			end_of_pol = find_file.index("PROJECTION BASIS SET")-1
		else:
			end_of_pol = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1

	lineN = 0
	t     = 1
	for i in range(find_file.index("POLARIZABLE POINTS")+1,end_of_pol,1):
		if t%4 == 1:
			v1,v2,v3 = neflmo[lineN][0]/neflmo[lineN][3],neflmo[lineN][1]/neflmo[lineN][3],neflmo[lineN][2]/neflmo[lineN][3]
			x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
			x1 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=21-len(find_file[i].split()[0]))
			x2 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=15)
			string += (x0+x1+x2+x3+'\n')
		elif t%4 == 2:
			v1,v2 = nefpol[lineN][0]/nefpol[lineN][9],nefpol[lineN][1]/nefpol[lineN][9]
			v3,v4 = nefpol[lineN][2]/nefpol[lineN][9],nefpol[lineN][3]/nefpol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=21)
			x2 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=15)
			x4 = '{:{align}{width}}'.format('%.10f'%v4,align='>',width=15)
			string += (x1+x2+x3+x4+' >'+'\n')
		elif t%4 == 3:
			v5,v6 = nefpol[lineN][4]/nefpol[lineN][9],nefpol[lineN][5]/nefpol[lineN][9]
			v7,v8 = nefpol[lineN][6]/nefpol[lineN][9],nefpol[lineN][7]/nefpol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v5,align='>',width=21)
			x2 = '{:{align}{width}}'.format('%.10f'%v6,align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%v7,align='>',width=15)
			x4 = '{:{align}{width}}'.format('%.10f'%v8,align='>',width=15)
			string += (x1+x2+x3+x4+' >'+'\n')
		else:
			v9 = nefpol[lineN][8]/nefpol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v9,align='>',width=21)
			lineN += 1
			string += (x1+'\n')
		t+= 1
	string += ' STOP\n'

	return string

def wrt_disp(find_file,nefdmo,nefdol,no_xr):
	string = ' DYNAMIC POLARIZABLE POINTS\n'
	if not no_xr:
		end_of_dyn = find_file.index("PROJECTION BASIS SET")-1
	else:
		end_of_dyn = find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1

	lineN = 0
	t     = 1
	for i in range(find_file.index("DYNAMIC POLARIZABLE POINTS")+1,end_of_dyn,1):
		if t%4 == 1:
			v1,v2,v3 = nefdmo[lineN][0]/nefdmo[lineN][3],nefdmo[lineN][1]/nefdmo[lineN][3],nefdmo[lineN][2]/nefdmo[lineN][3]
			x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
			x1 = '{:{align}{width}}'.format('%s'%find_file[i].split()[1],align='>',width=5-len(find_file[i].split()[0]))
			x2 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=23-len(x0)-len(x1))
			x3 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=15)
			x4 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=15)
			if int(find_file[i].split()[1]) == 1:
				line2 = find_file[i].split('--')
				x5 = '{:{align}{width}}'.format('%s'%line2[1],align='>',width=22)
				string += (x0+x1+x2+x3+x4+' --'+x5+'\n')
			else:
				string += (x0+x1+x2+x3+x4+'\n')
		elif t%4 == 2:
			v1,v2 = nefdol[lineN][0]/nefdol[lineN][9],nefdol[lineN][1]/nefdol[lineN][9]
			v3,v4 = nefdol[lineN][2]/nefdol[lineN][9],nefdol[lineN][3]/nefdol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=16)
			x2 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=16)
			x4 = '{:{align}{width}}'.format('%.10f'%v4,align='>',width=16)
			string += (x1+x2+x3+x4+' >'+'\n')
		elif t%4 == 3:
			v5,v6 = nefdol[lineN][4]/nefdol[lineN][9],nefdol[lineN][5]/nefdol[lineN][9]
			v7,v8 = nefdol[lineN][6]/nefdol[lineN][9],nefdol[lineN][7]/nefdol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v5,align='>',width=16)
			x2 = '{:{align}{width}}'.format('%.10f'%v6,align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%v7,align='>',width=16)
			x4 = '{:{align}{width}}'.format('%.10f'%v8,align='>',width=16)
			string += (x1+x2+x3+x4+' >'+'\n')
		else:
			v9 = nefdol[lineN][8]/nefdol[lineN][9]
			x1 = '{:{align}{width}}'.format('%.10f'%v9,align='>',width=16)
			lineN += 1
			string += (x1+'\n')
		t+= 1
	string += ' STOP\n'

	return string

def wrt_xr(find_file,parm_file,nefcmo,nefcof,xyz,fock_mat):
	string = ' PROJECTION BASIS SET\n'
	lineN  = 0
	c,w    = 0,0

	for i in range(find_file.index("PROJECTION BASIS SET")+1,find_file.index("MULTIPLICITY    1")-1,1):
		if "A" in find_file[i]:
			x0 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
			x1 = '{:{align}{width}}'.format('%.10f'%xyz[c][0],align='>',width=25-len(find_file[i].split()[0]))
			x2 = '{:{align}{width}}'.format('%.10f'%xyz[c][1],align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%xyz[c][2],align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%find_file[i].split()[4],align='>',width=7)
			c += 1
			string += (x0+x1+x2+x3+x4+'\n')
		else:
			string += parm_file[i]
	string += ' STOP\n MULTIPLICITY    1\n STOP\n'
	string += parm_file[find_file.index("MULTIPLICITY    1")+2]

	for i in range(find_file.index("MULTIPLICITY    1")+3,find_file.index("FOCK MATRIX ELEMENTS"),1):
		x0 = '{:{align}{width}}'.format('%s'%parm_file[i][0:2],align='>',width=2)
		x1 = '{:{align}{width}}'.format('%s'%parm_file[i][2:5],align='>',width=3)
		seek1,seek2,division = int(parm_file[i][0:2])-1,int(parm_file[i][2:5])*5,len(nefcof[0])-1

		if seek2-5 < len(nefcof[0])-1:
			v  = nefcof[seek1][seek2-5]/nefcof[seek1][division]
			x2 = '{:{align}{width}}'.format('%.8e'%v,align='>',width=15)
			string += x0+x1+x2
		if seek2-4 < len(nefcof[0])-1:
			v  = nefcof[seek1][seek2-4]/nefcof[seek1][division]
			x3 = '{:{align}{width}}'.format('%.8e'%v,align='>',width=15)
			string += x3
		if seek2-3 < len(nefcof[0])-1:
			v  = nefcof[seek1][seek2-3]/nefcof[seek1][division]
			x4 = '{:{align}{width}}'.format('%.8e'%v,align='>',width=15)
			string += x4
		if seek2-2 < len(nefcof[0])-1:
			v  = nefcof[seek1][seek2-2]/nefcof[seek1][division]
			x5 = '{:{align}{width}}'.format('%.8e'%v,align='>',width=15)
			string += x5
		if seek2-1 < len(nefcof[0])-1:
			v  = nefcof[seek1][seek2-1]/nefcof[seek1][division]
			x6 = '{:{align}{width}}'.format('%.8e'%v,align='>',width=15)
			string += x6
		string += '\n'

	for i in range(len(fock_mat)):
		string += fock_mat[i]
	string += ' LMO CENTROIDS\n'

	for i in range(find_file.index("LMO CENTROIDS")+1,find_file.index("SCREEN2      (FROM VDWSCL=   0.700)")-1,1):
		v1,v2,v3 = nefcmo[w][0]/nefcmo[w][3],nefcmo[w][1]/nefcmo[w][3],nefcmo[w][2]/nefcmo[w][3]
		x1 = '{:{align}{width}}'.format('%s'%find_file[i].split()[0],align='<',width=len(find_file[i].split()[0]))
		x2 = '{:{align}{width}}'.format('%.10f'%v1,align='>',width=23-len(find_file[i].split()[0]))
		x3 = '{:{align}{width}}'.format('%.10f'%v2,align='>',width=15)
		x4 = '{:{align}{width}}'.format('%.10f'%v3,align='>',width=15)
		w += 1
		string += (x1+x2+x3+x4+'\n')
	string += ' STOP\n'

	return string

def rotate_parameter(reflines,coord_xyz2,no_pol,no_disp,no_xr,fragname):
	parm_file,find_file                             = read_efp_param(reflines)
	# coor_xyz1,bmid_xyz1,ref_monop,ref_coor_dip,ref_bmid_dip,ref_coor_qup,ref_bmid_qup,ref_coor_ocp,ref_bmid_ocp,screen
	# Default, electrostatics
	default                                         = save_efp_coord_and_ele(parm_file,find_file,no_pol,no_disp,no_xr)
	atom_bmid                                       = len(default[0])+len(default[1])
	nefdip,nefqua,nefoct                            = np.zeros((atom_bmid,4)),np.zeros((atom_bmid,7)),np.zeros((atom_bmid,11))
	bmid_xyz2,connection                            = save_bond_connection(default[1],coord_xyz2)

	if not no_pol:
		ref_lmo,ref_pol                         = save_efp_pol(parm_file,find_file,no_pol,no_disp,no_xr)
		lmo_pts                                 = lmo_position(default[0],default[1],ref_lmo,len(default[0]),connection)
		neflmo,nefpol                           = np.zeros((len(ref_lmo),4)),np.zeros((len(ref_lmo),10))
	if not no_disp:
		dyn_lmo,ref_dyn                         = save_efp_dyn(parm_file,find_file,no_pol,no_disp,no_xr)
		temp = []
		for i in range(len(dyn_lmo)/12):
			temp.append(dyn_lmo[i])
		dyn_pts                                 = lmo_position(default[0],default[1],np.array(temp),len(default[0]),connection)
		nefdmo,nefdol                           = np.zeros((len(dyn_lmo),4)),np.zeros((len(dyn_lmo),10))
	if not no_xr:
		orbital,n_basis,ref_wvf,fock_mat,xr_lmo = save_efp_exr(parm_file,find_file,no_pol,no_disp,no_xr)
		xr_pts                                  = lmo_position(default[0],default[1],xr_lmo,len(default[0]),connection)
		nefcmo,nefcof                           = np.zeros((len(xr_lmo),4)),np.zeros((len(xr_lmo),sum(n_basis)+1))

	#=============================================================================#
	# Rotation Starts!                                                            #
	# Obtain three rotation matrices to transfer parameters during each iteration #
	#=============================================================================#

	for i in range(len(connection)):
		# Define Three Connected Atoms
		P1  = define_three_points(connection[i],default[0],0,1,2,1,2,3)
		P2  = define_three_points(connection[i],default[0],2,1,0,1,2,3)
		Q1  = define_three_points(connection[i],coord_xyz2,0,1,2,0,1,2)
		Q2  = define_three_points(connection[i],coord_xyz2,2,1,0,0,1,2)

		rot_left  = compute_rotation_matrix(P1,Q1,False)
		rot_right = compute_rotation_matrix(P2,Q2,False)
#		rot_cent  = compute_rotation_matrix(P1,Q1,True)

		nefdip    = rotate_dipole(nefdip,rot_left,default[3][connection[i][0]-1],connection[i][0]-1)
		nefdip    = rotate_dipole(nefdip,rot_left,default[4][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
		nefdip    = rotate_dipole(nefdip,rot_left,default[3][connection[i][1]-1],connection[i][1]-1)
		nefdip    = rotate_dipole(nefdip,rot_right,default[3][connection[i][1]-1],connection[i][1]-1)
		nefdip    = rotate_dipole(nefdip,rot_right,default[3][connection[i][2]-1],connection[i][2]-1)
		nefdip    = rotate_dipole(nefdip,rot_right,default[4][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))
		nefqua    = rotate_qupole(nefqua,rot_left,default[5][connection[i][0]-1],connection[i][0]-1)
		nefqua    = rotate_qupole(nefqua,rot_left,default[6][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
		nefqua    = rotate_qupole(nefqua,rot_left,default[5][connection[i][1]-1],connection[i][1]-1)
		nefqua    = rotate_qupole(nefqua,rot_right,default[5][connection[i][1]-1],connection[i][1]-1)
		nefqua    = rotate_qupole(nefqua,rot_right,default[5][connection[i][2]-1],connection[i][2]-1)
		nefqua    = rotate_qupole(nefqua,rot_right,default[6][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))
		nefoct    = rotate_ocpole(nefoct,rot_left,default[7][connection[i][0]-1],connection[i][0]-1)
		nefoct    = rotate_ocpole(nefoct,rot_left,default[8][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
		nefoct    = rotate_ocpole(nefoct,rot_left,default[7][connection[i][1]-1],connection[i][1]-1)
		nefoct    = rotate_ocpole(nefoct,rot_right,default[7][connection[i][1]-1],connection[i][1]-1)
		nefoct    = rotate_ocpole(nefoct,rot_right,default[7][connection[i][2]-1],connection[i][2]-1)
		nefoct    = rotate_ocpole(nefoct,rot_right,default[8][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))

#		nefdip    = rotate_dipole(nefdip,rot_left,default[3][connection[i][0]-1],connection[i][0]-1)
#		nefdip    = rotate_dipole(nefdip,rot_left,default[4][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
#		nefdip    = rotate_dipole(nefdip,rot_cent,default[3][connection[i][1]-1],connection[i][1]-1)
#		nefdip    = rotate_dipole(nefdip,rot_right,default[3][connection[i][2]-1],connection[i][2]-1)
#		nefdip    = rotate_dipole(nefdip,rot_right,default[4][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))
#
#		nefqua    = rotate_qupole(nefqua,rot_left,default[5][connection[i][0]-1],connection[i][0]-1)
#		nefqua    = rotate_qupole(nefqua,rot_left,default[6][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
#		nefqua    = rotate_qupole(nefqua,rot_cent,default[5][connection[i][1]-1],connection[i][1]-1)
#		nefqua    = rotate_qupole(nefqua,rot_right,default[5][connection[i][2]-1],connection[i][2]-1)
#		nefqua    = rotate_qupole(nefqua,rot_right,default[6][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))
#
#		nefoct    = rotate_ocpole(nefoct,rot_left,default[7][connection[i][0]-1],connection[i][0]-1)
#		nefoct    = rotate_ocpole(nefoct,rot_left,default[8][connection[i][3]-1],connection[i][3]-1+len(coord_xyz2))
#		nefoct    = rotate_ocpole(nefoct,rot_cent,default[7][connection[i][1]-1],connection[i][1]-1)
#		nefoct    = rotate_ocpole(nefoct,rot_right,default[7][connection[i][2]-1],connection[i][2]-1)
#		nefoct    = rotate_ocpole(nefoct,rot_right,default[8][connection[i][4]-1],connection[i][4]-1+len(coord_xyz2))

		temp1,temp2 = [],[]
		for j in range(3):
			x1,y1,z1 = default[0][connection[i][j]-1][1],default[0][connection[i][j]-1][2],default[0][connection[i][j]-1][3]
			x2,y2,z2 = coord_xyz2[connection[i][j]-1][0],coord_xyz2[connection[i][j]-1][1],coord_xyz2[connection[i][j]-1][2]
			temp1.append([float(x1),float(y1),float(z1)])
			temp2.append([float(x2),float(y2),float(z2)])
		
		if not no_pol:
			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_left,temp1[0],temp2[0],lmo_pts[i][0],ref_lmo,ref_pol)
			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_left,temp1[1],temp2[1],lmo_pts[i][1],ref_lmo,ref_pol)
			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_right,temp1[1],temp2[1],lmo_pts[i][1],ref_lmo,ref_pol)
			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_right,temp1[2],temp2[2],lmo_pts[i][2],ref_lmo,ref_pol)
		
		if not no_disp:
			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_left,temp1[0],temp2[0],dyn_pts[i][0],dyn_lmo,ref_dyn)
			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_left,temp1[1],temp2[1],dyn_pts[i][1],dyn_lmo,ref_dyn)
			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_right,temp1[1],temp2[1],dyn_pts[i][1],dyn_lmo,ref_dyn)
			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_right,temp1[2],temp2[2],dyn_pts[i][2],dyn_lmo,ref_dyn)
		
		if not no_xr:
			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_left,temp1[0],temp2[0],xr_pts[i][0],xr_lmo,orbital,n_basis,ref_wvf)
			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_left,temp1[1],temp2[1],xr_pts[i][1],xr_lmo,orbital,n_basis,ref_wvf)
			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_right,temp1[1],temp2[1],xr_pts[i][1],xr_lmo,orbital,n_basis,ref_wvf)
			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_right,temp1[2],temp2[2],xr_pts[i][2],xr_lmo,orbital,n_basis,ref_wvf)

#		if not no_pol:
#			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_left,temp1[0],temp2[0],lmo_pts[i][0],ref_lmo,ref_pol)
#			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_cent,temp1[1],temp2[1],lmo_pts[i][1],ref_lmo,ref_pol)
#			neflmo,nefpol = rotate_lmo_tensor(neflmo,nefpol,rot_right,temp1[2],temp2[2],lmo_pts[i][2],ref_lmo,ref_pol)
#		
#		if not no_disp:
#			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_left,temp1[0],temp2[0],dyn_pts[i][0],dyn_lmo,ref_dyn)
#			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_cent,temp1[1],temp2[1],dyn_pts[i][1],dyn_lmo,ref_dyn)
#			nefdmo,nefdol = rotate_dyn_tensor(nefdmo,nefdol,rot_right,temp1[2],temp2[2],dyn_pts[i][2],dyn_lmo,ref_dyn)
#		
#		if not no_xr:
#			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_left,temp1[0],temp2[0],xr_pts[i][0],xr_lmo,orbital,n_basis,ref_wvf)
#			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_cent,temp1[1],temp2[1],xr_pts[i][1],xr_lmo,orbital,n_basis,ref_wvf)
#			nefcmo,nefcof = rotate_projection_wf(nefcmo,nefcof,rot_right,temp1[2],temp2[2],xr_pts[i][2],xr_lmo,orbital,n_basis,ref_wvf)

#	outputfile = 'test1.efp'
#	ofile      = open(outputfile,'w')

	xyz = []
	for i in range(len(coord_xyz2)):
		xyz.append(coord_xyz2[i])
	for i in range(len(bmid_xyz2)):
		xyz.append(list(bmid_xyz2[i][1:]))
	
	new_parm = ''
	new_parm += wrt_ele(find_file,parm_file,xyz,nefdip,nefqua,nefoct,no_pol,no_disp,no_xr,fragname)

	if not no_pol and not no_disp and not no_xr:
		print("Rotate Electrostatic+Polarization+Dispersion+Exchange-Repulsion")
		new_parm += wrt_pol(find_file,neflmo,nefpol,no_disp,no_xr)
		new_parm += wrt_disp(find_file,nefdmo,nefdol,no_xr)
		new_parm += wrt_xr(find_file,parm_file,nefcmo,nefcof,xyz,fock_mat)
	elif not no_pol and not no_disp and no_xr:
		print("Rotate Electrostatic+Polarization+Dispersion")
		new_parm += wrt_pol(find_file,neflmo,nefpol,no_disp,no_xr)
		new_parm += wrt_disp(find_file,nefdmo,nefdol,no_xr)
	elif not no_pol and no_disp and not no_xr:
		print("Rotate Electrostatic+Polarization+Exchange-Repulsion")
		new_parm += wrt_pol(find_file,neflmo,nefpol,no_disp,no_xr)
		new_parm += wrt_xr(find_file,parm_file,nefcmo,nefcof,xyz,fock_mat)
	elif not no_pol and no_disp and no_xr:
		print("Rotate Electrostatic+Polarization")
		new_parm += wrt_pol(find_file,neflmo,nefpol,no_disp,no_xr)
	elif no_pol and not no_disp and not no_xr:
		print("Rotate Electrostatic+Dispersion+Exchange-Repulsion")
		new_parm += wrt_disp(find_file,nefdmo,nefdol,no_xr)
		new_parm += wrt_xr(find_file,parm_file,nefcmo,nefcof,xyz,fock_mat)
	elif no_pol and not no_disp and no_xr:
		print("Rotate Electrostatic+Dispersion")
		new_parm += wrt_disp(find_file,nefdmo,nefdol,no_xr)
	elif no_pol and no_disp and not no_xr:
		print("Rotate Electrostatic+Exchange-Repulsion")
		new_parm += wrt_xr(find_file,parm_file,nefcmo,nefcof,xyz,fock_mat)
	else:
		print("Rotate Electrostatic")

	for i in range(len(default[9])):
		new_parm += default[9][i]
	
	return new_parm
#=============================================================#
# combine two efp files ( should be written efficiently !! )  #
#=============================================================#
def find_line_number(parm):
	lineN = 0
	for i in range(len(parm.split('\n'))):
		line  = parm.split('\n')[i].strip()
		lineN += 1
		if 'COORDINATES (BOHR)' in line:
			coord1 = lineN
		elif 'MONOPOLES' in line:
			coord2 = lineN-1
			mon1   = lineN
		elif 'DIPOLES' in line:
			mon2   = lineN-1
			dip1   = lineN
		elif 'QUADRUPOLES' in line:
			dip2   = lineN-1
			qup1   = lineN
		elif 'OCTUPOLES' in line:
			qup2   = lineN-1
			ocp1   = lineN
		elif 'POLARIZABLE POINTS' in line:
			if 'DYNAMIC' not in line:
				ocp2 = lineN-1
				pol1 = lineN
			else:
				pol2 = lineN-1
				dyn1 = lineN
		elif 'PROJECTION BASIS SET' in line:
			dyn2 = lineN-1
			bas1 = lineN
		elif 'MULTIPLICITY' in line:
			bas2 = lineN-1
		elif 'PROJECTION WAVEFUNCTION' in line:
			wav1 = lineN
			lmon,basn = int(line.split()[2]),int(line.split()[3])
		elif 'FOCK MATRIX ELEMENTS' in line:
			wav2 = lineN
			foc1 = lineN
		elif 'LMO CENTROIDS' in line:
			foc2 = lineN
			lmo1 = lineN
		elif 'SCREEN' in line:
			if line.split()[0] == 'SCREEN2':
				lmo2  = lineN-1
				scr21 = lineN
			else:
				scr22 = lineN-1
				scr11 = lineN
		elif '$END' in line:
			scr12 = lineN-1
		else:
			continue

	return coord1,coord2,mon1,mon2,dip1,dip2,qup1,qup2,ocp1,ocp2,pol1,pol2,dyn1,dyn2,bas1,bas2,wav1,wav2,foc1,foc2,lmo1,lmo2,scr21,scr22,scr11,scr12,lmon,basn

def combine_parm(back_parm,side_parm,back_anum):

	b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b_lmo,b_bas = find_line_number(back_parm)
	s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s_lmo,s_bas = find_line_number(side_parm)

	parm_out = ''
	parm_out += combine_coordinate(back_parm.split('\n')[b1:b2-1],side_parm.split('\n')[s1:s2-1],back_anum)
	parm_out += combine_monopole(back_parm.split('\n')[b3:b4-1],side_parm.split('\n')[s3:s4-1],back_anum)
	parm_out += combine_dipole(back_parm.split('\n')[b5:b6-1],side_parm.split('\n')[s5:s6-1],back_anum)
	parm_out += combine_quadrupole(back_parm.split('\n')[b7:b8-1],side_parm.split('\n')[s7:s8-1],back_anum)
	parm_out += combine_octupole(back_parm.split('\n')[b9:b10-1],side_parm.split('\n')[s9:s10-1],back_anum)
	static_pol,lmo_cent = combine_polarizability(back_parm.split('\n')[b11:b12-1]+side_parm.split('\n')[s11:s12-1])
	parm_out += static_pol
	parm_out += combine_dynamic_polarizability(back_parm.split('\n')[b13:b14-1],side_parm.split('\n')[s13:s14-1])
	parm_out += combine_projection_basis_set(back_parm.split('\n')[b15:b16-1],side_parm.split('\n')[s15:s16-1],back_anum)
	parm_out += combine_wave_function(back_parm.split('\n')[b17:b18-1],side_parm.split('\n')[s17:s18-1],b_lmo,b_bas,s_lmo,s_bas)
	parm_out += combine_fock_matrix(back_parm.split('\n')[b19:b20-1],side_parm.split('\n')[s19:s20-1],b_lmo,s_lmo)
	parm_out += lmo_cent
	parm_out += combine_screen(back_parm.split('\n')[b23:b24-1],side_parm.split('\n')[s23:s24-1],back_anum,'SCREEN2')
	parm_out += combine_screen(back_parm.split('\n')[b25:b26-1],side_parm.split('\n')[s25:s26-1],back_anum,'SCREEN ')

	ofile = open('test4.efp','w')
	ofile.write(parm_out+' $END')

def rewrite_atom_number(atom,num):
	return 'A'+str(int(atom[1:3])+num).zfill(2)+atom[3:]

def rewrite_bmid_number(bmid,num):

	if len(bmid.split("BO")[1]) == 2:
		atom1,atom2 = bmid.split("BO")[1][0:1],bmid.split("BO")[1][1:2]
	elif len(bmid.split("BO")[1]) == 3:
		atom1,atom2 = bmid.split("BO")[1][0:2],bmid.split("BO")[1][2:3]
	else:
		if int(bmid.split("BO")[1][0:2]) > int(bmid.split("BO")[1][2:4]):
			atom1,atom2 = bmid.split("BO")[1][0:2],bmid.split("BO")[1][2:4]
		else:
			atom1,atom2 = bmid.split("BO")[1][0:3],bmid.split("BO")[1][3:4]
	return "BO"+str(int(atom1)+num)+str(int(atom2)+num)

def combine_coordinate(backbone,sidechain,backnum):
	string1,string2 = '',''
	string1 += '          RUNTYP=MAKEFP EFFECTIVE FRAGMENT POTENTIAL DATA FOLLOWS...\n'
	string1 += '          FRAGNAMEEFP GENERATED by ZEROTH ORDER APPROXIMATION\n'
	string1 += ' $FRAGNAME\n'
	string1 += 'EFP DATA FOR FRAGNAME SCFTYP=RHF     ... GENERATED WITH BASIS SET=XXX\n'
	string1 += ' COORDINATES (BOHR)\n'
	for i in range(len(backbone)):
		if 'BO' in backbone[i].split()[0]:
			bmid = rewrite_bmid_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%backbone[i].split()[4],align='>',width=12)
			x5 = '{:{align}{width}}'.format('%s'%backbone[i].split()[5],align='>',width=5)
			string2 += x0+x1+x2+x3+x4+x5+'\n'
		else:
			atom = rewrite_atom_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%backbone[i].split()[4],align='>',width=12)
			x5 = '{:{align}{width}}'.format('%s'%backbone[i].split()[5],align='>',width=5)
			string1 += x0+x1+x2+x3+x4+x5+'\n'
	for i in range(len(sidechain)):
		if 'BO' in sidechain[i].split()[0]:
			bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%sidechain[i].split()[4],align='>',width=12)
			x5 = '{:{align}{width}}'.format('%s'%sidechain[i].split()[5],align='>',width=5)
			string2 += x0+x1+x2+x3+x4+x5+'\n'
		else:
			atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%sidechain[i].split()[4],align='>',width=12)
			x5 = '{:{align}{width}}'.format('%s'%sidechain[i].split()[5],align='>',width=5)
			string1 += x0+x1+x2+x3+x4+x5+'\n'

	return string1+string2+' STOP\n'

def combine_monopole(backbone,sidechain,backnum):
	chg             = 0.0
	string1,string2 = '',''
	string1 += ' MONOPOLES\n'
	for i in range(len(backbone)):
		if 'BO' in backbone[i].split()[0]:
			bmid = rewrite_bmid_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.5f'%float(backbone[i].split()[2]),align='>',width=10)
			string2 += x0+x1+x2+'\n'
			chg += (float(backbone[i].split()[1])+float(backbone[i].split()[2]))
		else:
			atom = rewrite_atom_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.5f'%float(backbone[i].split()[2]),align='>',width=10)
			string1 += x0+x1+x2+'\n'
			chg += (float(backbone[i].split()[1])+float(backbone[i].split()[2]))
	for i in range(len(sidechain)):
		if 'BO' in sidechain[i].split()[0]:
			bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.5f'%float(sidechain[i].split()[2]),align='>',width=10)
			string2 += x0+x1+x2+'\n'
			chg += (float(sidechain[i].split()[1])+float(sidechain[i].split()[2]))
		else:
			atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.5f'%float(sidechain[i].split()[2]),align='>',width=10)
			string1 += x0+x1+x2+'\n'
			chg += (float(sidechain[i].split()[1])+float(sidechain[i].split()[2]))
	print('TOTAL FRAGMENT\'S CHARGE:'+'{:{align}{width}}'.format('%.5f'%float(chg),align='>',width=8))
	return string1+string2+' STOP\n'

def combine_dipole(backbone,sidechain,backnum):
	string1,string2 = '',''
	string1 += ' DIPOLES\n'
	for i in range(len(backbone)):
		if 'BO' in backbone[i].split()[0]:
			bmid = rewrite_bmid_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=24-len(bmid))
			x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=16)
			string2 += x0+x1+x2+x3+'\n'
		else:
			atom = rewrite_atom_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=24-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=16)
			string1 += x0+x1+x2+x3+'\n'
	for i in range(len(sidechain)):
		if 'BO' in sidechain[i].split()[0]:
			bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=24-len(bmid))
			x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=16)
			string2 += x0+x1+x2+x3+'\n'
		else:
			atom = rewrite_atom_number(sidechain[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=24-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=16)
			x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=16)
			string1 += x0+x1+x2+x3+'\n'
	return string1+string2+' STOP\n'

def combine_quadrupole(backbone,sidechain,backnum):
	string  = ''
	string += ' QUADRUPOLES\n'
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	for i in range(len(backbone)):
		if len(backbone[i].split()) == 6:
			if 'BO' in backbone[i].split()[0]:
				bmid = rewrite_bmid_number(backbone[i].split()[0],0)
				x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
				x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=24-len(bmid))
				x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[4]),align='>',width=16)
				temp2.append(x0+x1+x2+x3+x4+' >\n')
			else:
				atom = rewrite_atom_number(backbone[i].split()[0],0)
				x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
				x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=24-len(atom))
				x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[4]),align='>',width=16)
				temp4.append(x0+x1+x2+x3+x4+' >\n')
		else:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=16)
				temp2.append(x1+x2+'\n')
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=16)
				temp4.append(x1+x2+'\n')
				temp3.append(temp4)
				temp4 = []
			else:
				print('wrong quadrupole printed')
				exit()

	for i in range(len(sidechain)):
		if len(sidechain[i].split()) == 6:
			if 'BO' in sidechain[i].split()[0]:
				bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
				x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
				x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=24-len(bmid))
				x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[4]),align='>',width=16)
				temp2.append(x0+x1+x2+x3+x4+' >\n')
			else:
				atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
				x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
				x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=24-len(atom))
				x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[4]),align='>',width=16)
				temp4.append(x0+x1+x2+x3+x4+' >\n')
		else:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=16)
				temp2.append(x1+x2+'\n')
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=16)
				temp4.append(x1+x2+'\n')
				temp3.append(temp4)
				temp4 = []
			else:
				print('wrong quadrupole printed')

	qupol = temp3+temp1

	for i in range(len(qupol)):
		string += str(qupol[i][0])+str(qupol[i][1])

	return string + ' STOP\n'

def combine_octupole(backbone,sidechain,backnum):
	string = ' OCTUPOLES\n'
	temp1,temp2 = [],[]
	temp3,temp4 = [],[]
	for i in range(len(backbone)):
		if len(backbone[i].split()) == 6:
			if 'BO' in backbone[i].split()[0]:
				bmid = rewrite_bmid_number(backbone[i].split()[0],0)
				x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=24-len(bmid))
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[4]),align='>',width=16)
				temp2.append(x0+x1+x2+x3+x4+' >\n')
			else:
				atom = rewrite_atom_number(backbone[i].split()[0],0)
				x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=24-len(atom))
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[4]),align='>',width=16)
				temp4.append(x0+x1+x2+x3+x4+' >\n')
		elif len(backbone[i].split()) == 5:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[3]),align='>',width=16)
				temp2.append(x1+x2+x3+x4+' >'+'\n')
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[3]),align='>',width=16)
				temp4.append(x1+x2+x3+x4+' >'+'\n')
			else:
				print('wrong octupole printed')
				exit()
		else:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=16)
				temp2.append(x1+x2+'\n')
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=16)
				temp4.append(x1+x2+'\n')
				temp3.append(temp4)
				temp4 = []
			else:
				print("wrong octupole printed")
				exit()
	for i in range(len(sidechain)):
		if len(sidechain[i].split()) == 6:
			if 'BO' in sidechain[i].split()[0]:
				bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
				x0 = '{:{align}{width}}'.format('%s'%bmid,align='<',width=len(bmid))
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=24-len(bmid))
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[4]),align='>',width=16)
				temp2.append(x0+x1+x2+x3+x4+' >\n')
			else:
				atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
				x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=24-len(atom))
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[3]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[4]),align='>',width=16)
				temp4.append(x0+x1+x2+x3+x4+' >\n')
		elif len(sidechain[i].split()) == 5:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[3]),align='>',width=16)
				temp2.append(x1+x2+x3+x4+' >'+'\n')
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=16)
				x3 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=16)
				x4 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[3]),align='>',width=16)
				temp4.append(x1+x2+x3+x4+' >'+'\n')
			else:
				print('wrong octupole printed')
				exit()
		else:
			if len(temp2) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=16)
				temp2.append(x1+x2+'\n')
				temp1.append(temp2)
				temp2 = []
			elif len(temp4) > 0:
				x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[0]),align='>',width=24)
				x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=16)
				temp4.append(x1+x2+'\n')
				temp3.append(temp4)
				temp4 = []
			else:
				print("wrong octupole printed")
				exit()

	ocpol = temp3+temp1
	
	for i in range(len(ocpol)):
		string += str(ocpol[i][0])+str(ocpol[i][1])+str(ocpol[i][2])
	
	return string + ' STOP\n'

def combine_polarizability(pol_tensor):
	string1,string2 = ' POLARIZABLE POINTS\n',' LMO CENTROIDS\n'
	t,ct   = 1,1
	for i in range(len(pol_tensor)):
		if t % 4 == 1:
			lmop = 'CT'+str(ct)
			x0 = '{:{align}{width}}'.format('%s'%lmop,align='<',width=len(lmop))
			x1 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[1]),align='>',width=18-len(lmop))
			x2 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[3]),align='>',width=15)
			string1 += x0+x1+x2+x3+'\n'
			string2 += x0+x1+x2+x3+'\n'
		elif t % 4 == 2:
			x1 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[0]),align='>',width=18)
			x2 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[1]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[2]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[3]),align='>',width=15)
			string1 += x1+x2+x3+x4+' >\n'
		elif t%4 == 3:
			x1 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[0]),align='>',width=18)
			x2 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[1]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[2]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[3]),align='>',width=15)
			string1 += x1+x2+x3+x4+' >\n'
		else:
			x1 = '{:{align}{width}}'.format('%.10f'%float(pol_tensor[i].split()[0]),align='>',width=18)
			string1 += x1+'\n'
			ct      += 1
		t  += 1
	string1 += ' STOP\n'
	string2 += ' STOP\n'

	return string1,string2

def combine_dynamic_polarizability(backbone,sidechain):
	string = ' DYNAMIC POLARIZABLE POINTS\n'

	temp1,temp2 = [],[]
	for i in range(0,len(backbone),len(backbone)/12):
		for j in range(i,i+len(backbone)/12,1):
			temp2.append(backbone[j])
		temp1.append(temp2)
		temp2 = []
	temp3,temp4 = [],[]
	for i in range(0,len(sidechain),len(sidechain)/12):
		for j in range(i,i+len(sidechain)/12,1):
			temp4.append(sidechain[j].split(' --')[0])
		temp3.append(temp4)
		temp4 = []

	for i in range(12):
		ct = 1
		for j in range(len(temp1[i])):
			if temp1[i][j].split()[0] == 'CT':
				x0 = '{:{align}{width}}'.format('%s'%temp1[i][j].split()[0],align='<',width=len(temp1[i][j].split()[0]))
				x1 = '{:{align}{width}}'.format('%s'%str(ct),align='>',width=5-len(temp1[i][j].split()[0]))
				x2 = '{:{align}{width}}'.format('%.10f'%float(temp1[i][j].split()[2]),align='>',width=20-len(x0)-len(x1))
				x3 = '{:{align}{width}}'.format('%.10f'%float(temp1[i][j].split()[3]),align='>',width=15)
				x4 = '{:{align}{width}}'.format('%.10f'%float(temp1[i][j].split()[4]),align='>',width=15)
				ct += 1
				if int(temp1[i][j].split()[1]) == 1:
					line2 = temp1[i][j].split('--')
					x5 = '{:{align}{width}}'.format('%s'%line2[1],align='>',width=22)
					string += (x0+x1+x2+x3+x4+' --'+x5+'\n')
				else:
					string += (x0+x1+x2+x3+x4+'\n')
			else:
				string += str(temp1[i][j])+'\n'
		for j in range(len(temp3[i])):
			if temp3[i][j].split()[0] == 'CT':
				x0 = '{:{align}{width}}'.format('%s'%temp3[i][j].split()[0],align='<',width=len(temp3[i][j].split()[0]))
				x1 = '{:{align}{width}}'.format('%s'%str(ct),align='>',width=5-len(temp3[i][j].split()[0]))
				x2 = '{:{align}{width}}'.format('%.10f'%float(temp3[i][j].split()[2]),align='>',width=20-len(x0)-len(x1))
				x3 = '{:{align}{width}}'.format('%.10f'%float(temp3[i][j].split()[3]),align='>',width=15)
				x4 = '{:{align}{width}}'.format('%.10f'%float(temp3[i][j].split()[4]),align='>',width=15)
				ct += 1
				string += (x0+x1+x2+x3+x4+'\n')
			else:
				string += str(temp3[i][j])+'\n'

	return string+' STOP\n'

def combine_projection_basis_set(backbone,sidechain,backnum):
	string = ' PROJECTION BASIS SET\n'
	for i in range(len(backbone)):
		if len(backbone[i].split()) > 0 and backbone[i].split()[0][0:1] == 'A':
			atom = rewrite_atom_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[1]),align='>',width=25-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(backbone[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%backbone[i].split()[4],align='>',width=7)
			string += (x0+x1+x2+x3+x4+'\n')
		else:
			string += backbone[i]+'\n'

	for i in range(len(sidechain)):
		if len(sidechain[i].split()) > 0 and sidechain[i].split()[0][0:1] == 'A':
			atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[1]),align='>',width=25-len(atom))
			x2 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[2]),align='>',width=15)
			x3 = '{:{align}{width}}'.format('%.10f'%float(sidechain[i].split()[3]),align='>',width=15)
			x4 = '{:{align}{width}}'.format('%s'%sidechain[i].split()[4],align='>',width=7)
			string += (x0+x1+x2+x3+x4+'\n')
		else:
			string += sidechain[i]+'\n'

	return string+' STOP\n MULTIPLICITY    1\n STOP\n'

def combine_wave_function(backbone,sidechain,b_lmo,b_bas,s_lmo,s_bas):
	tot_lmo,tot_bas = b_lmo+s_lmo,b_bas+s_bas
	wave_function   = np.zeros((tot_lmo,tot_bas))
	string          = ' PROJECTION WAVEFUNCTION     '+str(tot_lmo)+'     '+str(tot_bas)+'\n'

	temp1,temp2 = [],[]
	for i in range(len(backbone)):
		if backbone[i][5:20] != '':
			temp1.append(backbone[i][5:20])
		if backbone[i][20:35] != '':
			temp1.append(backbone[i][20:35])
		if backbone[i][35:50] != '':
			temp1.append(backbone[i][35:50])
		if backbone[i][50:65] != '':
			temp1.append(backbone[i][50:65])
		if backbone[i][65:80] != '':
			temp1.append(backbone[i][65:80])
	for i in range(len(sidechain)):
		if sidechain[i][5:20] != '':
			temp2.append(sidechain[i][5:20])
		if sidechain[i][20:35] != '':
			temp2.append(sidechain[i][20:35])
		if sidechain[i][35:50] != '':
			temp2.append(sidechain[i][35:50])
		if sidechain[i][50:65] != '':
			temp2.append(sidechain[i][50:65])
		if sidechain[i][65:80] != '':
			temp2.append(sidechain[i][65:80])

	lineN = 0
	for i in range(b_lmo):
		for j in range(b_bas):
			wave_function[i][j] = temp1[lineN]
			lineN              += 1
	lineN = 0
	for i in range(b_lmo,tot_lmo,1):
		for j in range(b_bas,tot_bas,1):
			wave_function[i][j] = temp2[lineN]
			lineN              += 1

	for i in range(len(wave_function)):
		lineN = 1
		x0    = '{:{align}{width}}'.format('%s'%str(i+1),align='>',width=2)
		x1    = '{:{align}{width}}'.format('%s'%str(lineN),align='>',width=3)
		for j in range(len(wave_function[i])):
			if j != 0 and j % 5 == 0:
				string += x0+x1+'\n'
				lineN += 1
				x1 = '{:{align}{width}}'.format('%s'%str(lineN),align='>',width=3)
				x1 += '{:{align}{width}}'.format('%.8e'%wave_function[i][j],align='>',width=15)
			else:
				x1 += '{:{align}{width}}'.format('%.8e'%wave_function[i][j],align='>',width=15)
		string += x0+x1+'\n'

	return string

def combine_fock_matrix(backbone,sidechain,b_lmo,s_lmo):
	tot_lmo     = b_lmo+s_lmo
	temp1,temp2 = np.zeros((b_lmo,b_lmo)),[]
	for i in range(len(backbone)):
		for j in range(len(backbone[i].split())):
			if backbone[i].split()[j] == '>':
				continue
			else:
				temp2.append(float(backbone[i].split()[j]))
	for i in range(b_lmo+1):
		length = 0 + (i*(i-1)/2)
		for j in range(i):
			temp1[i-1][j] = temp2[j+length]
			temp1[j][i-1] = temp2[j+length]


	temp3,temp4 = np.zeros((s_lmo,s_lmo)),[]
	for i in range(len(sidechain)):
		for j in range(len(sidechain[i].split())):
			if sidechain[i].split()[j] == '>':
				continue
			else:
				temp4.append(float(sidechain[i].split()[j]))
	for i in range(s_lmo+1):
		length = 0 + (i*(i-1)/2)
		for j in range(i):
			temp3[i-1][j] = temp4[j+length]
			temp3[j][i-1] = temp4[j+length]

	fock_mat = np.zeros((tot_lmo,tot_lmo))
	for i in range(len(temp1)):
		for j in range(len(temp1[i])):
			fock_mat[i][j] = temp1[i][j]
	for i in range(len(temp3)):
		for j in range(len(temp3[i])):
			fock_mat[i+len(temp1)][j+len(temp1[0])] = temp3[i][j]


	lineN,string = 0,''
	string += ' FOCK MATRIX ELEMENTS\n'
	for i in range(len(fock_mat)):
		for j in range(i+1):
			x1      = '{:{align}{width}}'.format('%.10f'%float(fock_mat[i][j]),align='>',width=16)
			lineN  += 1
			string += x1
			if lineN == 4:
				string += ' >\n'
				lineN = 0

	return string+'\n'

def combine_screen(backbone,sidechain,backnum,scr):
	string1 = scr+'      (FROM VDWSCL=   0.700)\n'
	string2 = ''

	for i in range(len(backbone)):
		if 'BO' in backbone[i].split()[0]:
			bmid = rewrite_bmid_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%' '+bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=14)
			string2 += x0+x1+x2+'\n'
		else:
			atom = rewrite_atom_number(backbone[i].split()[0],0)
			x0 = '{:{align}{width}}'.format('%s'%' '+atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.9f'%float(backbone[i].split()[2]),align='>',width=14)
			string1 += x0+x1+x2+'\n'
	for i in range(len(sidechain)):
		if 'BO' in sidechain[i].split()[0]:
			bmid = rewrite_bmid_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%' '+bmid,align='<',width=len(bmid))
			x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=23-len(bmid))
			x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=14)
			string2 += x0+x1+x2+'\n'
		else:
			atom = rewrite_atom_number(sidechain[i].split()[0],backnum)
			x0 = '{:{align}{width}}'.format('%s'%' '+atom,align='<',width=len(atom))
			x1 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[1]),align='>',width=23-len(atom))
			x2 = '{:{align}{width}}'.format('%.9f'%float(sidechain[i].split()[2]),align='>',width=14)
			string1 += x0+x1+x2+'\n'

	return string1+string2+' STOP\n'

def main():
	parser = argparse.ArgumentParser(
	usage  ='%(prog)s [options] conf_a.efp conf_b.inp\nDefalut: Rotate Everything!!!!',
	formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-v', '--version', action='version', version='Flexible EFP ' + __version__)
	parser.add_argument('conf_a', metavar='conf_a', type=str, help='Conformation in .efp format (backbone)')
	parser.add_argument('conf_b', metavar='conf_b', type=str, help='Conformation in .inp format (Unknown)')
	parser.add_argument('-p', '--no-pol', action='store_true', help='Turn off the rotation of polarization')
	parser.add_argument('-d', '--no-disp', action='store_true', help='Turn off the rotation of dispersion')
	parser.add_argument('-xr', '--no-xr', action='store_true', help='Turn off the rotation of exchange repulsion')
	args = parser.parse_args()

	ifile1     = open(args.conf_a,'r')
	ifile2     = open(args.conf_b,'r')
	reflines1  = ifile1.readlines()
	reflines2  = ifile2.readlines()
	coord_xyz2 = read_new_struc(reflines2)

	frag  = args.conf_b.split('.inp')[0]
	frag_parm = rotate_parameter(reflines1,np.array(coord_xyz2),args.no_pol,args.no_disp,args.no_xr,frag)

	ofile1 = open(frag+'.efp','w')
	ofile1.write(frag_parm)

	exit()
	# additional
	combine_parm(back_parm,side_parm,len(back_xyz2))

if __name__ == "__main__":
	main()
