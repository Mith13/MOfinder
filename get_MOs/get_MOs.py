#!/usr/bin/python3
import math
import sys
from collections import deque

class MolInfo:
    def __init__(self,n_elec=0,mult=0,dim=0,found_seg=False):
        self.n_elec = (int)(n_elec)
        self.mult   = (int)(mult)
        self.dim    = (int)(dim)
        self.found  = found_seg
    def printInfo(self):
        print ("N of electrons  ... " + str(self.n_elec))
        print ("Multiplicity    ... " + str(self.mult))
        print ("N of basis sets ... " + str(self.dim))


#Fuck Python
def parseOutput(filename, atom, orb_type, thrs):
    seg_header  = deque(maxlen=3)
    found_seg   = False
    orb_num     = ""
    inputfile   = open(filename, 'r') 
    line        = inputfile.readline()
    orb_list    = []
    loewdin_seg = []
    n_elec      = 0
    mult        = 0
    dim         = 0
    while line: 
        if "Multiplicity           Mult            ...." in line:
            list1  = line.split()
            mult = list1[-1]
        if "Number of Electrons    NEL             ...." in line:
            list1  = line.split()
            n_elec = list1[-1]
        if "Basis Dimension        Dim             ...." in line:
            list1  = line.split()
            dim = list1[-1]
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
           found_seg = True
        if "MAYER POPULATION ANALYSIS" in line:
            break
        if found_seg is True:
           search_str = atom + " " + orb_type
           while "--------  --------  --------  --------  --------  --------" not in line and line:
               #print (line)
               if search_str in line:
                   #print(line)
                   data = line.split()
                   loewdin_seg.append(line) 
                   for i in range(3,len(data)):
                       if (float)(data[i]) > thrs:
                           orb_data = []
                           orb_data.append(orb_num[i-3])
                           orb_data.append(data[0])
                           orb_data.append(data[1])
                           orb_data.append(data[2])
                           orb_data.append(data[i])
                           orb_list.append(orb_data)
               seg_header.append(line)
               line = inputfile.readline()
           orb_num = seg_header[0].split()
        line = inputfile.readline()
    orb_list.sort(key=lambda x: int(x[0]))
    mol_info = MolInfo(n_elec,mult,dim,found_seg)
    if found_seg is True:
       for dat in orb_list:
#           print("Orbital "+ '\033[91m' + dat[0] +  '\033[0m' + " has contribution from " + '\033[91m'+ dat[1]+dat[2] + " " + dat[3] + '\033[0m' + " equal " + '\033[92m' + dat[4] + '\033[0m')
           print("Orbital "+ dat[0] +  " has contribution from " + dat[1]+dat[2] + " " + dat[3] + " equal "  + dat[4])
    else:
        print ("There are no data to analyse");
    return orb_list,mol_info

def prepareCharmolOrbs(orb_list,mol_info):
    orbs = []
    for entry in orb_list:
        orbs.append(entry[0])

    unique_orbs  = set(orbs)
    orbs         = list(unique_orbs)
    orbs.sort(key=int)

    charmol_occ  = "#Occupied\nisovalue 0.06\nalpha "
    charmol_virt = "#Virtual\nisovalue 0.03\nalpha "
    max_occ      = (int)(mol_info.n_elec/2+mol_info.mult%2)
    iorb         = 0
    for iorb in range(0,len(orbs)):
        if (int)(orbs[iorb]) < max_occ:
            charmol_occ+=str(int(orbs[iorb])+1)+" "
        else:
            break

    for iorb2 in range(iorb,len(orbs)):
        charmol_virt+=str(int(orbs[iorb2])+1)+" "
    print ("\nOrbital windows for Charmol to export\n")
    charmol_occ+=" griddensity low"
    print (charmol_occ)
    charmol_virt+=" griddensity low\n"
    print(charmol_virt)



#if len(sys.argv) < 4:
#    print("There is a missing parameter \n Example get_MOs.py file atom_name orbital_type threshold(optional)")
#    exit()
thrs = 10.0
#filename = sys.argv[1]
#atom     = sys.argv[2]
#orb_type = sys.argv[3]
#thrs     = (float)(sys.argv[4])

filename = "cu2o2.0.0.loc000_dft_small.out.bcp"
atom="Cu"
orb_type="d"

if thrs < 5.0:
    thrs = 5.0
if len(atom) < 2:
    atom+=" " * (2 - len(atom))

print("Will printout orbitals with the population of " + orb_type + " on atom/s " + atom + " higher than " + str(thrs) +"%")
orbitals,mol_info  = parseOutput(filename,atom,orb_type,thrs)
prepareCharmolOrbs(orbitals,mol_info)
mol_info.printInfo()

