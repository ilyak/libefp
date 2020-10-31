import sys,os

def main(inputfile1,inputfile2):
	pwd       = os.getcwd()+'/'
	ifile1    = open(inputfile1,'r') # for elec
	ifile2    = open(inputfile2,'r') # for rest
	reflines1 = ifile1.readlines()
	reflines2 = ifile2.readlines()

	lineN = 0
	for obj in reflines1:
		line  = obj.strip()
		lineN += 1
		if 'POLARIZABLE' in line and 'DYNAMIC' not in line:
			elec1 = lineN
		elif 'SCREEN2' in line:
			elec2 = lineN
		else:
			continue
	print(inputfile1,inputfile2)
	lineN = 0
	for obj in reflines2:
		line  = obj.strip()
		lineN += 1
		if 'POLARIZABLE' in line and 'DYNAMIC' not in line:
			rest1 = lineN
		elif 'SCREEN2' in line:
			rest2 = lineN
		else:
			continue

	lineN = 0
	parm1,parm3 = [],[]
	for obj in reflines1:
		lineN += 1
		if lineN < elec1:
			if lineN == 4:
				parm1.append('6-31G* elec, 6-31+xG(3df,2p) - else\n')
			else:
				parm1.append(obj)
		elif lineN >= elec2:
			parm3.append(obj)
		else:
			continue

	lineN = 0
	parm2 = []
	for obj in reflines2:
		lineN += 1
		if lineN >= rest1 and lineN < rest2:
			parm2.append(obj)
		else:
			continue


	efparm = parm1+parm2+parm3
	output = inputfile1.split('/')[1]
	ofile = open(output,'w')
	for i in range(len(efparm)):
		ofile.write(efparm[i])

main(sys.argv[1],sys.argv[2])
