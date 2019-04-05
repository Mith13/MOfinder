#!/usr/bin/python3
import math
import sys
from collections import deque

################ MOLECULE INFO CLASS ###########################
class MolInfo:
    def __init__ (self,n_elec=0,mult=0,dim=0,found_seg=False):
        self.n_elec = (int)(n_elec)
        self.mult   = (int)(mult)
        self.dim    = (int)(dim)
        self.found  = found_seg

    def printInfo (self):
        print ("N of electrons  ... " + str(self.n_elec))
        print ("Multiplicity    ... " + str(self.mult))
        print ("N of basis sets ... " + str(self.dim))
        print ("")

################ MO CLASS #########################
class MO:
    def __init__ (self,header=None,pop=None,orb_number=-1,en=0):
        self.header = header
        self.pop    = [float(i) for i in pop]
        self.orb_n  = (int)(orb_number)


        if all(len(i) != len(self.header[0]) for i in self.header):
            print("There is some problem with creating MO. Different sizes for header: " + str(len(self.header[0])) + " " + str(len(self.header[1])) + " " + str(len(self.header[2])))
            exit()
        if len(self.pop) != len(self.header[0]):
            print("There is some problem with creating MO. Different sizes for lists;  header:" +  str(len(self.header[0])) + " C: " + str(len(self.pop)))
            exit()

    def printMO (self):
        if self.orb_n == 1:
            print("Problem iwth
        print ("Printing orbital " + '\033[91m' + str(self.orb_n) +  '\033[0m')
        for i in range(0,len(self.pop)):
            if self.pop[i] > 5:
                print (self.header[0][i] + " " + self.header[1][i] + " " + self.header[2][i] + " " + '\033[92m' + str(self.pop[i]) +  '\033[0m')
        print("")



#Fuck Python
def parseOutput(filename, atom, orb_type, thrs):
    seg_header  = deque(maxlen=3)
    empty_line  = False
    found_seg   = False
    orb_num     = ""
    inputfile   = open(filename, 'r')
    line        = inputfile.readline()
    orb_list    = list()
    MOs         = list()
    n_elec      = 0
    mult        = 0
    dim         = 0
    while line:
        if not line.strip():
            if empty_line is True:
                found_seg = False
            empty_line = True:
        else :
            empty_line = False
        # Parse basic info
        if "Multiplicity           Mult            ...." in line:
            list1  = line.split()
            mult = list1[-1]
        if "Number of Electrons    NEL             ...." in line:
            list1  = line.split()
            n_elec = list1[-1]
        if "Basis Dimension        Dim             ...." in line:
            list1  = line.split()
            dim = list1[-1]
        # Parse MO segment
        if "MOLECULAR ORBITALS" in line:
            found_seg = True

        if found_seg is True:
           idx2save = set()
           loewdin_seg = list()
           save_MOs = False
           search_str = atom + " " + orb_type
           while "--------  --------  --------  --------  --------  --------" not in line and line:
               #print (line)
               data = line.split()
               loewdin_seg.append(data)
               if search_str in line:
                   #print(line)
                   for i in range(3,len(data)):
                       if (float)(data[i]) > thrs:
                           save_MOs = True
#                           orb_data = []
#                           orb_data.append(orb_num[i-3])
#                           orb_data.append(data[0])
#                           orb_data.append(data[1])
#                           orb_data.append(data[2])
#                           orb_data.append(data[i])
                           orb_list.append([orb_num[i-3],data[0],data[1],data[2],data[i]])
                           idx2save.add(i)
               seg_header.append(line)
               line = inputfile.readline()

           if save_MOs is True:
               loewdin_seg = loewdin_seg[:len(loewdin_seg)-5]
               loewdin_seg = list(map(list,zip(*loewdin_seg)))
               MO_header = [loewdin_seg[0],loewdin_seg[1],loewdin_seg[2]]

               for i in idx2save:
                  MOs.append(MO(MO_header,loewdin_seg[i],orb_num[i-3]))

           orb_num = seg_header[0].split()
           loewdin_seg.clear()
        line = inputfile.readline()


        # Parse Loewdin segment
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
           found_seg = True

        if "MAYER POPULATION ANALYSIS" in line:
            break
        if found_seg is True:
           idx2save = set()
           loewdin_seg = list()
           save_MOs = False
           search_str = atom + " " + orb_type
           while "--------  --------  --------  --------  --------  --------" not in line and line:
               #print (line)
               data = line.split()
               loewdin_seg.append(data)
               if search_str in line:
                   #print(line)
                   for i in range(3,len(data)):
                       if (float)(data[i]) > thrs:
                           save_MOs = True
#                           orb_data = []
#                           orb_data.append(orb_num[i-3])
#                           orb_data.append(data[0])
#                           orb_data.append(data[1])
#                           orb_data.append(data[2])
#                           orb_data.append(data[i])
                           orb_list.append([orb_num[i-3],data[0],data[1],data[2],data[i]])
                           idx2save.add(i)
               seg_header.append(line)
               line = inputfile.readline()

           if save_MOs is True:
               loewdin_seg = loewdin_seg[:len(loewdin_seg)-5]
               loewdin_seg = list(map(list,zip(*loewdin_seg)))
               MO_header = [loewdin_seg[0],loewdin_seg[1],loewdin_seg[2]]

               for i in idx2save:
                  MOs.append(MO(MO_header,loewdin_seg[i],orb_num[i-3]))

           orb_num = seg_header[0].split()
           loewdin_seg.clear()
        line = inputfile.readline()
    # sort according to MO number
    orb_list.sort(key=lambda x: int(x[0]))
    MOs.sort(key=lambda x: x.orb_n)
    # create MolInfo
    mol_info = MolInfo(n_elec,mult,dim,found_seg)
    if found_seg is True:
       for dat in orb_list:
           print("Orbital "+ '\033[91m' + dat[0] +  '\033[0m' + " has contribution from " + '\033[91m'+ dat[1]+dat[2] + " " + dat[3] + '\033[0m' + " equal " + '\033[92m' + dat[4] + '\033[0m')
    else:
        print ("There are no data to analyse");
    return orb_list,mol_info,MOs

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



if len(sys.argv) < 4:
    print("There is a missing parameter \n Example get_MOs.py file atom_name orbital_type threshold(optional)")
    exit()
thrs = 10.0
filename = sys.argv[1]
atom     = sys.argv[2]
orb_type = sys.argv[3]
thrs     = (float)(sys.argv[4])
if thrs < 5.0:
    thrs = 5.0
if len(atom) < 2:
    atom+=" " * (2 - len(atom))

print("Will printout orbitals with the population of " + orb_type + " on atom/s " + atom + " higher than " + str(thrs) +"%")
orbitals,mol_info,MOs  = parseOutput(filename,atom,orb_type,thrs)
prepareCharmolOrbs(orbitals,mol_info)
mol_info.printInfo()
for orb in MOs:
    orb.printMO()

