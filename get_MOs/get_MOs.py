import math
import sys
from collections import deque


#Global variables ????  Fuck Python
def get_MOpop(filename, atom, orb_type, thrs):
    seg_header = deque(maxlen=3)
    found_seg = False
    orb_num = ""
    inputfile = open(filename, 'r') 
    line = inputfile.readline()
    while line: 
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
           found_seg = True
        if found_seg is True:
                #print (line)
                search_str = atom + " " + orb_type
                while "--------  --------  --------  --------  --------  --------" not in line and line:
                    #print (line)
                    if search_str in line:
                        #print(line)
                        data = line.split()
                        for i in range(3,len(data)):
                            if (float)(data[i]) > thrs:
                                print("Orbital " + orb_num[i-3] + " has contribution from " + data[0]+data[1] + " " + data[2] + " equal " + data[i])
                    seg_header.append(line)
                    line = inputfile.readline()
                orb_num = seg_header[0].split()
        line = inputfile.readline()


print(len(sys.argv))
if len(sys.argv) < 4:
    print("There is a missing parameter \n Example get_MOs.py file atom_name orbital_type threshold(optional)")
    #exit()
thrs = 10.0
#filename = sys.argv[1]
#atom = sys.argv[2]
#orb_type = sys.argv[3]
#thrs = (float)(sys.argv[4])
filename = "cu2o2.0.0.loc000_dft_small.out.bcp"
atom = "Cu"
orb_type = "d"
if thrs < 5.0:
    thrs = 5.0
if len(atom) < 2:
    atom+=" " * (2 - len(atom))

print("Will printout orbitals with the population of " + orb_type + " on atom " + atom + " higher than " + str(thrs))
get_MOpop(filename,atom,orb_type,thrs)

