#!/usr/bin/python3
import math
import sys
import re
from collections import deque

################ MOLECULE INFO CLASS ###########################
# Some basic info about system from ORCA's output. Mainly good to 
# determine Fermi vacuum
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
# Base class, used for loewdin reduced MO analysis.
# Stores only orbitals with orb content higher than use specified threshold
class MO(object):
    def __init__ (self,header=None,pop=None,orb_number=-1,en=0,isBeta=False):
        self.header = header
        if pop is not None:
            self.pop    = [float(i) for i in pop]
        self.orb_n  = (int)(orb_number)
        self.energy = (float)(en)
        self.isBeta = isBeta

        if all(len(i) != len(self.header[0]) for i in self.header):
            header_str = ""
            for item in self.header:
                header_str+=str(len(item))+" "
            print("There is some problem with creating MO. Different sizes for header: " + header_str)
            exit()
        if len(self.pop) != len(self.header[0]):
            print("There is some problem with creating MO. Different sizes for lists;  header:" +  str(len(self.header[0])) + " C: " + str(len(self.pop))) 
            exit()

    def printMO (self,orb_atom,orb_type):   
        thres = 5
        if self.orb_n == -1:
            print("Problem with MO number")
            exit()
        print ("Printing orbital E({:3d}): {:9.6f}".format(self.orb_n,self.energy))
        for i in range(0,len(self.pop)):
            if abs(self.pop[i]) > thres:
                if self.isBeta is True:
                    header_str = "Beta "
                else:
                    header_str = ""
                for item in self.header:
                    header_str+=item[i]+" "
                if orb_atom in header_str and orb_type in header_str:
                    print ("{:1.10s} \033[92m{:4.1f}\033[0m %".format(header_str,self.pop[i]))
                else:
                    print ("{:1.10s} {:4.1f} %".format(header_str,self.pop[i]))
        print("")


################ AO_MO CLASS #########################
# Derived class from MO, stores AO_MO cooeficients higher
# than hardcored value, default is 0.00001
# Basically same as base class
class AO_MO(MO):

    def printMO (self):   
        thres = 0.001
        if self.orb_n == -1:
            print("Problem with MO number")
            exit()
        print ("Printing orbital E({:3d}): \033[91m{:9.6f}\033[0m ".format(self.orb_n,self.energy))
        for i in range(0,len(self.pop)):
            if abs(self.pop[i]) > thres:
                if self.isBeta is True:
                    header_str = "Beta "
                else:
                    header_str = ""
                for item in self.header:
                    header_str+=item[i]+" "

                print (header_str + " " + '\033[92m' + str(self.pop[i]) +  '\033[0m')
                #print ("{:<10} {:<12.6}".format(header_str,self.pop[i]))
        print("")

    def orb_content (self, orb_name, orb_atom, thres):
        orb_sum  = []
        total    = 0
        print_it = False 
        for i in range(0,len(self.pop)):
            total += self.pop[i]*self.pop[i]
            orb_sum.append(self.pop[i]*self.pop[i])
        
        ratio = [x/total*100 for x in orb_sum]
        for i in range(0,len(ratio)):
            if orb_name in self.header[1][i] and orb_atom in self.header[0][i] and ratio[i] > thres:
                print_it = True
        if print_it is True:
            print ("Printing orbital E({:3d}): {:9.6f} ".format(self.orb_n,self.energy))
            for i in range(0,len(ratio)):
                if ratio[i] > 5:
                    tmp = re.split('(\D+)',self.header[0][i])
                    if self.isBeta is True:
                        header_str = "\033[91m "
                    else:
                        header_str = " "
                    header_str = header_str + tmp[0] + " " + tmp[1] + " " + self.header[1][i]
                    if orb_atom in header_str and orb_type in header_str:
                        print ("{:1.10s} \033[92m{:4.1f}\033[0m %".format(header_str,ratio[i]))
                    else:
                        print ("{:1.10s} {:4.1f}\033[0m %".format(header_str,ratio[i]))
                    #print ("Orbital {self.heade} {} content is {:4.2f}".format(self.orb_n,orb_name,ratio))
        if print_it is True:
            if self.orb_n is None:
                print ("THE FUCK IS WRONG WITH YOU!!!!")
            return self.orb_n

#Find two empty lines after matrix print, is same in every ORCA output,
# unless they will break it in future
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
    
# Nice pythonic way to compare two lists :)
def compare_lists(list1,list2):
    return list(set(list2)-set(list1))

#def find_permutations(list1,list2):
#    compare_lists(list1,list2)
#    return a


#Fuck Python
def parseOutput(filename, atom, orb_type, thrs,thrsAOMO):
# Bunch of local variables and switches
# It works in such way, that it always stores last three lines
# in the memory. This is because coeeficient and data starts with
# a pattern ----- ---- etc. and orbital numbers are written out
# in 3 lines earlier
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
    isBeta       = False
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
           # Find where orbital coeffs starts
           while "                  --------  " not in line and line :
               empty_line,empty_2lines=empty_lines(empty_line,line)
               if empty_2lines is True:
                   break
               if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
                   break
               #Obviously ORCA is stupid and sometimes there is not space between two coeffs
               data = line.split()
               coeffs = re.findall(r'[0-9]+.[0-9]{6}',line)
               data = data[0:2]
               data.extend(coeffs)
               read_seg.append(data)
               if len(data) > 0 and atom in data[0] and orb_type in data[1]:
                   for i in range(2,len(data)):
                       if (float)(data[i]) > thrsAOMO:
                           save_MOs = True
                           # Save orbital number, atom, orb_type and coefficients
                           orb_list.append([orb_num[i-2],data[0],data[1],data[i]])
                           # also save an absolute column index
                           idx2save.add(i)
               seg_header.append(line)
               line = inputfile.readline()
           # ORCA prints in six columns blocks, here we take care of distinguishing between
           # last or not last block
           if save_MOs is True:
               if empty_2lines is True:
                   read_seg = read_seg[:len(read_seg)-1]
                   isBeta = True
               else:
                   read_seg = read_seg[:len(read_seg)-3]
               read_seg = list(map(list,zip(*read_seg)))
               MO_header = [read_seg[0],read_seg[1]]
               for i in idx2save:
                  # create instance for every MO above threshold 
                  MOs_AO_MO.append(AO_MO(MO_header,read_seg[i],orb_num[i-2],orb_en[i-2],isBeta))

           orb_num = seg_header[0].split()
           orb_en  = seg_header[1].split()

        ######################################
        # Parse Loewdin segment
        ######################################
        # Same things as in previous section
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
           print ("Starting to analyze the Loewdin population per MO")
           found_pop = True
           #reset lists and bools
           found_mos = False
           isBeta    = False
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
                    isBeta = True
               else:
                    read_seg = read_seg[:len(read_seg)-5]
               read_seg = list(map(list,zip(*read_seg)))
               MO_header = [read_seg[0],read_seg[1],read_seg[2]]

               for i in idx2save:
                  MOs_pop.append(MO(MO_header,read_seg[i],orb_num[i-3],orb_en[i-3],isBeta))

           orb_num = seg_header[0].split()
           orb_en  = seg_header[1].split()
           read_seg.clear() 
        line = inputfile.readline()
    inputfile.close()
    # sort according to MO number
    orb_list.sort(key=lambda x: int(x[0]))
    MOs_pop.sort(key=lambda x: x.orb_n)
    MOs_AO_MO.sort(key=lambda x: x.orb_n)
    
    # create MolInfo 
    mol_info = MolInfo(n_elec,mult,dim,found_pop or found_mos)
    
    return orb_list,mol_info,MOs_pop,MOs_AO_MO

# creates orbitals for charmol, needs it in 1...99999 format
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
        if (int)(orbs[iorb]) < max_occ-1:
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
    orb.printMO(atom,orb_type)

mos_in_AOMOs = []
print ("Printing all orbitals with content of " + orb_type + " higher than " + str(thrs/2))
for orb in AOMOs:
    #orb.printMO()
    #this is stupid way to do it, there will be a lot of 'None's but I am too lazy to do it correctly 
    mos_in_AOMOs.append(orb.orb_content(orb_type,atom,thrs/2))
print ("Orbitals not common for both methods")
extra_MOs=compare_lists([item.orb_n for item in MOs], [item for item in mos_in_AOMOs if item is not None]) # here goes the bug fix for None
extra_MOs.extend(compare_lists([item for item in mos_in_AOMOs if item is not None], [item.orb_n for item in MOs]))
extra_MOs.sort(key=lambda x: int(x))
print (extra_MOs)
#for n in extra_MOs:
#    orb=next((item for item in MOs if n == item.orb_n),None)
#    if orb is not None:
#        orb.printMO()
#    orb=next((item for item in AOMOs if n == item.orb_n),None)
#    if orb is not None:
#        orb.orb_content('d')

