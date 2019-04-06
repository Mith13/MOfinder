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

################ base MO CLASS #########################
class MO(object):
    def __init__ (self,header=None,pop=None,orb_number=-1,en=0):
        self.header = header
        self.pop    = [float(i) for i in pop]
        self.orb_n  = (int)(orb_number)
        self.energy = (float)(en)

        if all(len(i) != len(self.header[0]) for i in self.header):
            header_str = ""
            for item in self.header:
                header_str+=str(len(item))+" "
            print("There is some problem with creating MO. Different sizes for header: " + header_str)
            exit()
        if len(self.pop) != len(self.header[0]):
            print("There is some problem with creating MO. Different sizes for lists;  header:" +  str(len(self.header[0])) + " C: " + str(len(self.pop))) 
            exit()

    def printMO (self):   
        thres = 5
        if self.orb_n == -1:
            print("Problem with MO number")
            exit()
        print ("Printing orbital E(\033[91mv{:3d}v\033[0m): {:9.6f}".format(self.orb_n,self.energy))
        for i in range(0,len(self.pop)):
            if abs(self.pop[i]) > thres:
                header_str = ""
                for item in self.header:
                    header_str+=item[i]+" "

                print ("{:1.10s} \033[92m {:4.1f} \033[0m %".format(header_str,self.pop[i]))
        print("")


################ AO_MO CLASS #########################
class AO_MO(MO):

    def printMO (self):   
        thres = 0.001
        if self.orb_n == -1:
            print("Problem with MO number")
            exit()
        print ("Printing orbital E(\033[91m {:3d} \033[0m): {:9.6f}".format(self.orb_n,self.energy))
        for i in range(0,len(self.pop)):
            if abs(self.pop[i]) > thres:
                header_str = ""
                for item in self.header:
                    header_str+=item[i]+" "

                #print (header_str + " " + '\033[92m' + str(self.pop[i]) +  '\033[0m')
                print ("{:<10} \033[92m {:<12.6} \033[0m".format(header_str,self.pop[i]))
        print("")
    def orb_content (self, orb_name):
        orb_sum = 0
        total   = 0
        for i in range(0,len(self.pop)):
            total += self.pop[i]*self.pop[i]
            if orb_name in self.header[1][i]:
                orb_sum += self.pop[i]*self.pop[i]
        ratio = orb_sum/total*100
        if ratio > 10:
            print ("Orbital {} {} content is {:4.2f}".format(self.orb_n,orb_name,ratio))
        return ratio
# Auxialary functions 
def empty_lines(empty,line):
    empty_line   = False
    empty_2lines = False
    if not line.strip():
        if empty is True:
            empty_2lines = True 
        empty_line = True
    else:
            empty_line = False
            empty_2lines = False
    return empty_line,empty_2lines
    
def compare_lists(list1,list2):
    return list(set(list2)-set(list1))

#Fuck Python
def parseOutput(filename, atom, orb_type, thrs,thrsAOMO):
    seg_header   = deque(maxlen=3)
    empty_line  = False
    empty_2lines = False
    found_pop    = False
    found_mos    = False
    orb_num      = ""
    orb_en       = ""
    inputfile    = open(filename, 'r') 
    line         = inputfile.readline()
    orb_list     = list() 
    MOs_pop      = list() 
    MOs_AO_MO    = list() 
    n_elec       = 0
    mult         = 0
    dim          = 0

    while line:
        if "FINAL SINGLE POINT ENERGY" in line:
            break
        
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
        
        ######################################
        # Parse MO segment
        ######################################
        if "MOLECULAR ORBITALS" in line:
           print ("Starting to analyze the MO coefficients")
           found_mos = True

        if found_mos is True:
           idx2save = set() 
           read_seg = list()
           save_MOs = False
           #while "--------  --------  --------  --------  --------  --------" not in line and line:
           while "                  --------  " not in line and line :
               empty_line,empty_2lines=empty_lines(empty_line,line)
               if empty_2lines is True:
                   break
               if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
                   break
               data = line.split()
               read_seg.append(data)
               if len(data) > 0 and atom in data[0] and orb_type in data[1]:
                   for i in range(2,len(data)):
                       if (float)(data[i]) > thrsAOMO:
                           save_MOs = True
                           orb_list.append([orb_num[i-3],data[0],data[1],data[i]])
                           idx2save.add(i)
               seg_header.append(line)
               line = inputfile.readline()

           if save_MOs is True:
               if empty_2lines is True:
                   read_seg = read_seg[:len(read_seg)-1]
               else:
                   read_seg = read_seg[:len(read_seg)-3]
               read_seg = list(map(list,zip(*read_seg)))
               MO_header = [read_seg[0],read_seg[1]]
               for i in idx2save:
                  MOs_AO_MO.append(AO_MO(MO_header,read_seg[i],orb_num[i-2],orb_en[i-2]))

           orb_num = seg_header[0].split()
           orb_en  = seg_header[1].split()

        ######################################
        # Parse Loewdin segment
        ######################################
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
           print ("Starting to analyze the Loewdin population per MO")
           found_pop = True
           #reset lists and bools
           found_mos = False
           orb_list.clear()
        
        if "MAYER POPULATION ANALYSIS" in line:
            break
        if found_pop is True:
           idx2save = set() 
           read_seg = list()
           save_MOs = False
           search_str = atom + " " + orb_type
           while "                  --------  " not in line and line :
               empty_line,empty_2lines=empty_lines(empty_line,line)
               if empty_2lines is True:
                   break
               data = line.split()
               read_seg.append(data) 
               if search_str in line:
                   for i in range(3,len(data)):
                       if (float)(data[i]) > thrs:
                           save_MOs = True
                           orb_list.append([orb_num[i-3],data[0],data[1],data[2],data[i]])
                           idx2save.add(i)
               seg_header.append(line)
               line = inputfile.readline()

           if save_MOs is True:
               if empty_2lines is True:
                    read_seg = read_seg[:len(read_seg)-1]
               else:
                    read_seg = read_seg[:len(read_seg)-5]
               read_seg = list(map(list,zip(*read_seg)))
               MO_header = [read_seg[0],read_seg[1],read_seg[2]]

               for i in idx2save:
                  MOs_pop.append(MO(MO_header,read_seg[i],orb_num[i-3],orb_en[i-3]))

           orb_num = seg_header[0].split()
           orb_en  = seg_header[1].split()
           read_seg.clear() 
        line = inputfile.readline()

    # sort according to MO number
    orb_list.sort(key=lambda x: int(x[0]))
    MOs_pop.sort(key=lambda x: x.orb_n)
    MOs_AO_MO.sort(key=lambda x: x.orb_n)
    
    # create MolInfo 
    mol_info = MolInfo(n_elec,mult,dim,found_pop or found_mos)
    
    return orb_list,mol_info,MOs_pop,MOs_AO_MO

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

orbitals,mol_info,MOs,AOMOs  = parseOutput(filename,atom,orb_type,thrs,1e-5)

prepareCharmolOrbs(orbitals,mol_info)
mol_info.printInfo()

for orb in MOs:
    orb.printMO()

print ("You should also considerate these orbitals")
extra_MOs=compare_lists([item.orb_n for item in MOs], [item.orb_n for item in AOMOs])
for n in extra_MOs:
    orb=next((item for item in MOs if n == item.orb_n),None)
    if orb is not None:
        orb.printMO()
    orb=next((item for item in AOMOs if n == item.orb_n),None)
    if orb is not None:
        orb.orb_content('d')
