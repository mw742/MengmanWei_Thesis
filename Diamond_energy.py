#################################################################
## Acyclic alkanes on diamond framework
##
## Input an InChI of an acyclic alkane with at least three carbons
##
## Input of just an InChI will do a conformation search which
## lists all low-energy conformations and an estimate of their
## energies, approximately in kJ/mol.
##
## Input of an InChI and a conformation number will report the
## energy and generate an .xyz file of the structure of that
## conformation
##
#################################################################

import sys
import array
import time

start = time.process_time()


def printatoms():
  print("                I#,At,Di,Tw,Nd, x, y, z, Co ")
  for i in range(number_of_carbons):
    print("%3d " % (i), atom[i*9:i*9+9])

def printxyzatoms():
  # this is an XYZ file
  print(number_of_carbons)
  print(sys.argv[1])
  scalefactor = 0.85
  for i in range(number_of_carbons):
    print("C   %8.4f  %8.4f  %8.4f" % (float(atom[i*9+5])*0.885,float(atom[i*9+6])*0.885,float(change_enantiomer*atom[i*9+7])*0.885))

def printmolatoms():
  # NB must be in .mol order, not InChI order, for a molfile
  # This subroutine needs more work
  print(sys.argv[1])
  print()
  print()
  print("%3d%3d  0  0  1  0  0  0  0  0999 V2000" % (number_of_carbons,number_of_carbons-1))
  scalefactor = 0.85
  for i in range(number_of_carbons):
    print("%8.4f  %8.4f  %8.4f C   0  0  0  0  0  0  0  0  0  0  0  0" % (float(atom[i*9+5])*0.885,float(atom[i*9+6])*0.885,float(atom[i*9+7])*0.885))
  for i in range(number_of_carbons-1):
    print("%3d%3d  1  0" % (inchi2wn[atom[i*9+12]],inchi2wn[atom[i*9+9]]))
  print("M  END")
  print()
  print("$$$$")

def atominfo(atomnumber,information):
  # Return data for atom 'atomnumber'
  # 0: InchiNumber
  # 1: Previous Attached Atom
  # 2: Was Direction; now energy of atom
  # 3: Twist
  # 4: New Direction
  # 5,6,7: x,y,z coordinates
  # 8: conformation - number of the relevant torsion angle
  # print("atominfo working",atomnumber,information,atomnumber*8+information)
  return atom[atomnumber*9+information]

def setatominfo(atomnumber,information,datapoint):
  atom[atomnumber*9+information]=datapoint
   
#################################################################
## From a conformation number, recalculate newdirection and
## calculate the x, y, z coordinates of the conformation
## The current conformation is define by list conformation[]
## NB: direction (atominfo(#, 2) is currently unused
## Rest of system is unchanged
## old_conformation[] has the list for the energies in the starting structure
## except first time around
def find_all_conf_number():
  # energy:
  # count number of overlaps, adjacent atoms, 1,5 diaxial, gauche
  num_overlap=0
  num_adjacent=0
  num_15diaxial=0
  num_gauche=0
  num_axeq=0
  energy_total=0
  # Energy components
  energy_axeq=-2
  energy_gauche=2
  energy_15diaxial=20
  energy_adjacent=200
  energy_overlap=energy_adjacent
  current_atom_energy=0
  current_atom_marker=0
  num_repeat_torsions=-1
  nt=0
  while nt < number_rotatable_bonds:
    if conformation[nt]==old_conformation[nt]:
      num_repeat_torsions+=1
    else:
      nt=number_rotatable_bonds
    nt+=1
  #print("Number of repeated torsions",num_repeat_torsions+1,conformation,old_conformation)
  #print("Conformation",conformation)
  ia=0
  

  while ia < number_of_carbons:
    current_atom_marker=ia
    current_torsion=atominfo(ia,8)
    ##if current_torsion > -1:
    if current_torsion > -1 and current_torsion >= num_repeat_torsions:
      if atominfo(inchi2wn[atominfo(inchi2wn[atominfo(ia,1)],1)],1)==0:
        total_twist=0
      else:
        # This is a torsion which can change
        #print(current_torsion)
        #print(ia)
        total_twist=(atominfo(ia,3)+conformation[current_torsion])%3
      #print("### total_twist, line 119,total_twist, atominfo(ia,3), conformation[current_torsion], ia, current_torsion")
      #print("### total_twist, line 119",total_twist, atominfo(ia,3), conformation[current_torsion], ia, current_torsion)
      # find direction of penultimate atom
      direction=atominfo(inchi2wn[atominfo(inchi2wn[atominfo(ia,1)],1)],4)
      #print("direction:", direction)
      # find direction of attached atom
      old_direction=atominfo(inchi2wn[atominfo(ia,1)],4)
      #print("old_direction:", old_direction)
      if twisting[total_twist*64+old_direction*8+direction]>8:
        #print("conf_num, conformation",conf_num,conformation)
        #print("current_torsion",current_torsion,"num_repeat_torsions",num_repeat_torsions)
        #print("atom, attach_atom",ia,atominfo(ia,0),atominfo(ia,3),inchi2wn[atominfo(ia,1)])
        #print("t-1, t, t+1",twisting[total_twist*64+old_direction*8+direction-1],twisting[total_twist*64+old_direction*8+direction],twisting[total_twist*64+old_direction*8+direction+1])
        #print("twisting",twisting[total_twist*64+old_direction*8+direction],"tt,od,d",total_twist,old_direction,direction,"combine",total_twist*64+old_direction*8+direction)
        #printatoms()
        #printxyzatoms()
        print("****** twisting too big ****** ","current_atom_marker:",current_atom_marker,"old_direction:",old_direction,"direction:",direction,"total_twist:",total_twist)
        print(twisting[total_twist*64+old_direction*8+direction])
      newdirection=twisting[total_twist*64+old_direction*8+direction]
      #print("newdirection:", newdirection)
      iattach_atom=inchi2wn[atominfo(ia,1)]
      #print("starting on atom",ia)
      #print(atominfo(ia,0),atominfo(ia,1),atominfo(ia,2),atominfo(ia,3),atominfo(ia,4),atominfo(ia,5),atominfo(ia,6),atominfo(ia,7),atominfo(ia,8))
      #print("dir, old, new, iattach_atom",direction,old_direction,newdirection,iattach_atom,"tt",total_twist)
      xcoord = atominfo(iattach_atom,5)+directions[newdirection*3]
      ycoord = atominfo(iattach_atom,6)+directions[newdirection*3+1]
      zcoord = atominfo(iattach_atom,7)+directions[newdirection*3+2]
      #print("old x,y,z",atominfo(iattach_atom,5),atominfo(iattach_atom,6),atominfo(iattach_atom,7),"directions",newdirection,";",directions[newdirection*3],directions[newdirection*3+1],directions[newdirection*3+2])
      #print("new x,y,z",xcoord,ycoord,zcoord)
      setatominfo(ia,4,newdirection)
      setatominfo(ia,5,xcoord)
      setatominfo(ia,6,ycoord)
      setatominfo(ia,7,zcoord)
      if ia > 2:
        # No energy contribution from first three atoms, which are joined to each other
        # Now calculate the contributions to the energy
        # compare distance with all previous atoms
        # However, can we avoid Pythagoras, because distances so close?
        # overlap - idential coordinates
        # adjacent - one bond length (root 3) away, and not bonded: dx+dy+dz=3, but so does 0,0,3
        # 1,5 - two bond lengths: 2root2, and not bonded: dx+dy+dz=4
        # gauche - three bonds (does bonding matter? probably not) (0,0,0 via 2,2,0 to 1,3,-1) root11
        # make a list of atoms sorted by x, y, z?
        # Should butane and methylpropane have the same energy? Give an energy for methyl, methylene, methine, C ?
        # Compile a list of connections (1,3; 1,4) for each atom?
        # favourable non-bonded interactions?
        # print("worrying about energy")
        current_atom_energy=0
        ea=0
        while ea < ia:
        #for ea in range(0, ia):
          #if ea>=ia:
          #  print("ea>=ia",ea,ia)
          xea = atominfo(ea,5)-xcoord
          yea = atominfo(ea,6)-ycoord
          zea = atominfo(ea,7)-zcoord
          distance2 = xea*xea+yea*yea+zea*zea
          #print("distance2",ia,ea,"d2",distance2,"x,y,z",xea,yea,zea)
          if distance2==0:
            num_overlap=num_overlap+1
            # no point in going on - energy too high
            energy_total=-1-ia
            test_atom=current_atom_marker
            while test_atom < number_of_carbons:
              if atominfo(test_atom,8) < current_torsion:
                ia=current_atom_marker
                break
              else:
                ia=number_of_carbons
                ea=ia
              test_atom+=1            
            current_atom_energy+=energy_overlap
          if distance2==3 and not inchi2wn[atominfo(ia,1)]==ea:
            #print("adjacent",ia,ea,atominfo(ia,1),inchi2wn[atominfo(ia,1)])
            num_adjacent=num_adjacent+1
            # no point in going on - energy too high
            energy_total=-1-ia
            test_atom=current_atom_marker
            while test_atom < number_of_carbons:
              if atominfo(test_atom,8) < current_torsion:
                ia=current_atom_marker
                break
              else:
                ia=number_of_carbons
                ea=ia
              test_atom+=1            
            current_atom_energy+=energy_adjacent
          if distance2==8:
            # two bonds away, or 1,5 diaxial?
            # if two bonds away, atominfo(ia,1) must be attached to ea
            if not (atominfo(ia,1)==atominfo(ea,1) or inchi2wn[atominfo(inchi2wn[atominfo(ia,1)],1)]==ea):
              # if they were to be bound to each other
              # print("1,3 connection",ia,ea,";",atominfo(ia,1),atominfo(ea,1),";",inchi2wn[atominfo(inchi2wn[atominfo(ia,1)],1)],ea)
              # Good moment to look for gauche interactions???
              # What is attached to ea, other than atominfo(ia,1) ?
              # Might be helpful to make a list in advance
              num_15diaxial=num_15diaxial+1
              current_atom_energy+=energy_15diaxial
          if distance2==11:
            # possible gauche interaction; only if connected??? Rather unusual if not connected; Risk it!
            num_gauche=num_gauche+1
            current_atom_energy+=energy_gauche
          if distance2==16:
            num_axeq=num_axeq+1
            current_atom_energy+=energy_axeq
          ea+=1
      if energy_total > -1:
        energy_total+=current_atom_energy
        #print("current_atom_marker,ia,current_atom_energy",current_atom_marker,ia,current_atom_energy,energy_total)
      setatominfo(current_atom_marker,2,current_atom_energy)
    else:
      energy_total+=atominfo(current_atom_marker,2)
      #print("Adding energy from previous structure",ia,current_atom_marker,atominfo(current_atom_marker,2),energy_total)
      #printatoms()
      #print("===")
    ia+=1
  #if energy_total > -1:
    # gauche interaction is about 2 kJ/mol (two butane conformations)
    # 1,5 diaxial is at least 20 kJ/mol - this is a relaxed diaxial cyclohexane, so should be more
    # overlap and adjacent are very large - no more detail needed
    #print("Energy_total",energy_total)
    #energy_total=energy_gauche*num_gauche+energy_15diaxial*num_15diaxial+energy_adjacent*(num_overlap+num_adjacent)
    print("Energy of",conf_num,":",num_overlap,num_adjacent,num_15diaxial,num_gauche,num_axeq,"overall",energy_total,conformation)
    #printatoms()
    #else:
    #print("*** Energy overload ***",conf_num,energy_total,-1-energy_total,conformation)
  print()
  return energy_total




#################################################################
## Start and calculate preliminary quantities:
## Number of carbons, stereogenic centres, rotatable bonds, etc
#################################################################

print("##########################################")
print("## Diamond Lattice Energy, version 1.0  ##")
print("##  Jonathan Goodman, jmg11@cam.ac.uk   ##")
print("##       Mengman Wei, mw742@cam.ac.uk   ##")
print("##   University of Cambridge, 2018-21   ##")
print("##########################################")
print(sys.argv[1])

number_of_carbons=int(sys.argv[1].split("/")[1].split("H")[0].split("C")[1])
#print("Number of carbons: ",number_of_carbons)

change_enantiomer=1
if sys.argv[1].find("/t") > -1:
  stereogenic_centres=sys.argv[1].split("/")[4].split("t")[1]
  number_stereogenic_centres=1+stereogenic_centres.count(",")
  list_stereogenic_centres=stereogenic_centres.split(",")
  list_numbers_stereogenic_centres=array.array("i")
  for i in range(0, len(list_stereogenic_centres)):
    list_numbers_stereogenic_centres.append(int(list_stereogenic_centres[i].replace("+","").replace("-","")))
  #print(number_stereogenic_centres,"stereogenic_centres:",stereogenic_centres)
  #print(list_stereogenic_centres)
  if sys.argv[1].find("/m1") > -1:
    change_enantiomer=-1
else:
  stereogenic_centres=''
  number_stereogenic_centres=0
connectivity=sys.argv[1].split("/")[2].split("c")[1]
print("connectivity ",connectivity)
list_of_connectivity=connectivity.replace(")",")-").replace("(","-(").replace(",","-,").split("-")
print("list_of_connectivity ",list_of_connectivity)

# inchi2workingnumber translates the inchi atom number to the position in the array
# InChI numbers are 1 - n; working numbers are 0 to n-1
# The atoms are numbered as InChI numbers, and as the position in the atom list
# For example, if an InChI has c1-6(2)7(3,4)5
# Then: inchi2wn[6]=1
# and:  wn2inchi[1]=6
# Then: inchi2wn[7]=3
# and:  wn2inchi[3]=7
inchi2wn=array.array('i')
wn2inchi=array.array('i')
tmp_inchi2workingnumber=connectivity.replace(")","-").replace("(","-").replace(",","-").split("-")
for workingposition in range(0, number_of_carbons):
  wn2inchi.append(int(tmp_inchi2workingnumber[workingposition]))
  inchi2wn.append(0)
inchi2wn.append(0)
for workingposition in range(0, number_of_carbons):
  inchi2wn[wn2inchi[workingposition]]= workingposition
tmp_inchi2workingnumber.clear()
#print("wn2inchi", wn2inchi)
#print("inchi2wn", inchi2wn[1: number_of_carbons+1])

# Number of methyl groups; may not need this, but never know...
#number_methyl=int(sys.argv[1].split("/")[3].split("h")[1].split(",").pop().split("-").pop().split("H")[0])
#print(number_methyl,"methyls")

# Number of rotatable torsions is number of bonds after methyl groups have been removed
# which is one less than the number of atoms after all methyl groups have gone
# This does not allow for symmetry of tBu, and other symmetrical structures
number_rotatable_bonds=len(connectivity.replace(")","-").split("-"))-3
number_conformations=3**number_rotatable_bonds
#print(connectivity.replace(")","-").split("-"))
# if number_rotatable_bonds == 0:
#     print("Only one conformation: nothing to do")
#     exit()
# else:
print(number_rotatable_bonds,"rotatable bonds;",number_conformations,"conformations")
conformation=array.array("i")
old_conformation=array.array("i")
conf_number=0
#print("len(sys.argv)",len(sys.argv),sys.argv)
if len(sys.argv) > 2:
  conf_number=int(sys.argv[2])
tmp_conf_number=conf_number
if number_rotatable_bonds == 0:
    conformation.append(0)
else:
    for i in range(0,number_rotatable_bonds):
      conformation.append(0)
old_conformation=conformation[:]
old_conformation[0]=-1
for i in range(0,number_rotatable_bonds):
  #print("conf_sort_out_setup",number_rotatable_bonds,i,tmp_conf_number,conformation)
  conformation[number_rotatable_bonds-i-1]=tmp_conf_number%3
  tmp_conf_number=int(tmp_conf_number/3)
#print("Conformation",conformation)

# The directions array gives the vector for each step between atoms
directions=array.array('i',
  [ 1, 1, 1,
   -1,-1, 1,
   -1, 1,-1,
    1,-1,-1,
   -1,-1,-1,
    1, 1,-1,
    1,-1, 1,
   -1, 1, 1 ])
#print(directions)

# The twistings array gives the twisted vectors
# The 9 should never occur
# With penultimate direction  pendir
# and attached atom direction attdir
# the new direction will be pendir (untwisted)
# twisting[attdir*8+pendir] to twist up
# twisting[64+attdir*8+pendir] to twist down
# newdirection=twisting[twist*64+old_direction*8+direction]
twisting=array.array('i',
  [ 9, 9, 9, 9,  9, 5, 6, 7,
    9, 9, 9, 9,  4, 9, 6, 7,
    9, 9, 9, 9,  4, 5, 9, 7,
    9, 9, 9, 9,  4, 5, 6, 9,
    9, 1, 2, 3,  9, 9, 9, 9,
    0, 9, 2, 3,  9, 9, 9, 9,
    0, 1, 9, 3,  9, 9, 9, 9,
    0, 1, 2, 9,  9, 9, 9, 9,

    9, 9, 9, 9,  9, 7, 5, 6,
    9, 9, 9, 9,  6, 9, 7, 4,
    9, 9, 9, 9,  7, 4, 9, 5,
    9, 9, 9, 9,  5, 6, 4, 9,
    9, 2, 3, 1,  9, 9, 9, 9,
    3, 9, 0, 2,  9, 9, 9, 9,
    1, 3, 9, 0,  9, 9, 9, 9,
    2, 0, 1, 9,  9, 9, 9, 9,

    9, 9, 9, 9,  9, 6, 7, 5,
    9, 9, 9, 9,  7, 9, 4, 6,
    9, 9, 9, 9,  5, 7, 9, 4,
    9, 9, 9, 9,  6, 4, 5, 9,
    9, 3, 1, 2,  9, 9, 9, 9,
    2, 9, 3, 0,  9, 9, 9, 9,
    3, 0, 9, 1,  9, 9, 9, 9,
    1, 2, 0, 9,  9, 9, 9, 9 ])
#print(twisting)



#################################################################
# Assemble molecule
# First three atoms are all connected as propane
# Each atom has the InChI atom number,
# the InChI# of the atom to which it is attached,
# the bond direction (from list of eight possibilities),
# twist (0,1,2) for multiple atoms on same last bond
# the overall new direction
# x,y,z coordinates (signed integers)
# NB: first atom is one, for easy array indexing
# I#,At,Di,Tw,Nd, x, y, z
#
# Need to index flexible torsion angles
# Each uniquely defined by terminal atom, but need a list
#################################################################
atom1=int(list_of_connectivity[0])
atom2=int(list_of_connectivity[1])
attach=[atom1,atom2]
atom=array.array('i',
  [int(list_of_connectivity[0]),0,0,0,5,0,0,0,-1,
   int(list_of_connectivity[1]),wn2inchi[0],0,0,0,1,1,1,-1])
if not list_of_connectivity[2].isnumeric():
  if list_of_connectivity[3].find(",") > -1:
      atom.extend([int(list_of_connectivity[2].replace("(","").replace(")","")),wn2inchi[1],0,1,7,0,2,2,-1])
      atom.extend([int(list_of_connectivity[3].replace("(","").replace(")","").replace(",","")),wn2inchi[1],0,2,6,2,0,2,-1])
      if list_of_connectivity[3].find(")") == -1:
          attach.append(int(list_of_connectivity[3].replace("(","").replace(")","")))
      else:
          atom.extend([int(list_of_connectivity[4]),wn2inchi[1],0,0,5,2,2,0,-1])
          attach.append(int(list_of_connectivity[4]))  
  else:
      atom.extend([int(list_of_connectivity[2].replace("(","").replace(")","")),wn2inchi[1],0,1,7,0,2,2,-1])
      if list_of_connectivity[2].find(")") > -1:
          atom.extend([int(list_of_connectivity[3]),wn2inchi[1],0,0,5,2,2,0,-1])
          attach.append(int(list_of_connectivity[3]))
      else:
          attach.append("(")
          attach.append(int(list_of_connectivity[2].replace("(","").replace(")","")))      
else:
  atom.extend([int(list_of_connectivity[2]),wn2inchi[1],0,0,5,2,2,0,-1])
  attach.append(int(list_of_connectivity[2]))


## current_position is the end of the fixed atoms and the start of the flexible molecule
## NB - 9 is the number of integers per atom
current_position=int(len(atom)/9)
#print("current_position ",current_position)

#################################################################
# Now add rest of structure.
# Assume, for the moment, structure zero, so all torsion angles are 180
# Probably hold conformation as base three string:
#################################################################
#print("attach",attach)

# To sort out each atom, need to find:
# old_direction = actual direction of attached atom: atominfo(inchi2wn[attach_atom],4)
# direction = actual direction of penultimate atom: atominfo(inchi2wn[atominfo(inchi2wn[attach_atom],1)],4)
# new direction = direction + affect of twist, so actual direction = newdirection=twisting[twist*64+old_direction*8+direction]

highest_torsion_angle=-1
#list_of_torsion_angles=[]
#list_of_torsion_angles.append(current_torsion_angle)
#print("list_of_torsion_angles",list_of_torsion_angles)
# The array conformation defines whether the angle should be 180, 60 or -60
# conformation[list_of_torsion_angles[current_torsion_angle]] will be zero, one or two

for workingposition in range(current_position, number_of_carbons):
  #print("workingposition loop", workingposition,wn2inchi[workingposition],list_of_connectivity[workingposition],"  ",attach)
  current_atom=int(list_of_connectivity[workingposition].replace("(","").replace(")","").replace(",",""))
  # Choose arbitrary conformation for base structure; number zero? As extended as possible
  # print("direction",direction,"x,y,z",directions[direction*3],directions[direction*3+1],directions[direction*3+2])
  # What atom (InChI#) is this atom attached to? Put in attach_atom
  if list_of_connectivity[workingposition].isnumeric():
    # This atom is attached to the top atom in list attach
    attach_atom = int(attach[-1])
    attach.append(current_atom)
    # not a branch, so twist is zero 0 in this construction phase - conformation is all zeros
    old_direction=atominfo(inchi2wn[attach_atom],4)
    penultimate_atom=atominfo(inchi2wn[attach_atom],1)
    direction=atominfo(inchi2wn[penultimate_atom],4)
    twist=0
    newdirection = direction
    #print("workingposition,direction, new, twist",workingposition,direction, newdirection)
  else:
    #print("not numeric", workingposition, current_atom)
    # This could be open bracket - so a branch
    # Close bracket, so end of a branch
    # Comma, so alternative branch
    # Or a combination of these
    if list_of_connectivity[workingposition].find("(") > -1:
      #print("Found '('")
      attach_atom = int(attach[-1])
      attach.append("(")
      attach.append(current_atom)
      penultimate_atom=atominfo(inchi2wn[attach_atom],1)
      direction=atominfo(inchi2wn[penultimate_atom],4)
      old_direction=atominfo(inchi2wn[attach_atom],4)
      twist=1
      # NB this arbitrarily assigns stereochemistry; sort out after assigning whole structure
      #print("workingposition; direction, old; twist",workingposition,"; ",direction, old_direction,"; ",twist)
      if list_of_connectivity[workingposition].find(")") > -1:
        #print(workingposition," attach ",attach)
        attach.pop()
        attach.pop()
      #print("attach_atom:", attach_atom)
    else:
      #print("'(' not found; try ',' amd ')'")
      if list_of_connectivity[workingposition].find(",") > -1 or list_of_connectivity[workingposition].find(")") > -1:
        # Go back to last "("
        #
        # Three situations:
        #  ) - go back to last branch
        #  , - attach to last branch but maintain current
        # ,) - attach to last branch and go back to last branch
        #
        # For example, if attach is [ 1, 13, 17, (, 14, (, 4 ] and current atom is 5
        #
        # If the latest character is 5) then should be attached to 4 and close branch
        #   so attach becomes [ 1, 13, 17, (, 14 ] ; attach_atom is 4
        #
        # If the latest character is ,5 then should be attached to 14 and maintain branch
        #   so attach becomes [ 1, 13, 17, (, 14, (, 5 ] ; attach_atom is 14
        #
        # If the latest character is ,5) then should be attached to 14 and close branch
        #   so attach becomes [ 1, 13, 17, (, 14 ] ; attach_atom is 14
        #
        attach_atom = int(attach[-1])
        top_item=attach.pop()
        while not top_item=="(":
          top_item=attach.pop()
          #print("workingposition,attach",workingposition,attach)
        if list_of_connectivity[workingposition].find(",") > -1:
          attach_atom = int(attach[-1])
          twist=2
          if list_of_connectivity[workingposition].find(")") == -1:
            # comma and not )
            attach.append("(")
            attach.append(current_atom)
        else:
          twist=0
        # attach should now be sorted out
        penultimate_atom=atominfo(inchi2wn[attach_atom],1)
        direction=atominfo(inchi2wn[penultimate_atom],4)
        old_direction=atominfo(inchi2wn[attach_atom],4)
      else:
        print("There is no 'else:'; should never be here")

  # Find current torsion angle
  if highest_torsion_angle < 0:
    current_torsion_angle=0
    highest_torsion_angle=0
  else:
    # New torsion, or re-use old one? Has attached atom been used before?
    #print("torsions",current_torsion_angle,highest_torsion_angle,attach_atom,"current_position,workingposition",current_position,workingposition)
    test_current_torsion_angle = -1
    for check_at in range(current_position, workingposition):
      if atominfo(check_at,1) == attach_atom:
        # reuse torsion
        test_current_torsion_angle = atominfo(check_at,8)
        current_torsion_angle = atominfo(check_at,8)
    if test_current_torsion_angle < 0:
        if atominfo(inchi2wn[attach_atom],8) == -1:
            current_torsion_angle = 0
        else:
          # Need new torsion angle
          highest_torsion_angle += 1
          current_torsion_angle = highest_torsion_angle
          #print("highest_torsion_angle",highest_torsion_angle,"current_torsion_angle",current_torsion_angle)

  # coordinates are those of attached atom, plus direction
  # printatoms()
  # print("twist,old_direction,direction",twist,old_direction,direction)
  newdirection=twisting[twist*64+old_direction*8+direction]
  #print("xyz coords",attach_atom,inchi2wn[attach_atom],"old, direction, new, twist",old_direction,direction,newdirection,twist)
  xcoord = atominfo(inchi2wn[attach_atom],5)+directions[newdirection*3]
  ycoord = atominfo(inchi2wn[attach_atom],6)+directions[newdirection*3+1]
  zcoord = atominfo(inchi2wn[attach_atom],7)+directions[newdirection*3+2]
  # we have a sidechain; find post-sidechain atom
  #atom.extend([current_atom, direction, twist, attach_atom, xcoord, ycoord, zcoord], conformation)
  #if not attach_atom in list_of_torsion_angles:
  #  list_of_torsion_angles.append(attach_atom)
  #current_torsion_angle = list_of_torsion_angles.index(attach_atom)
  #print("current_torsion_angle",current_torsion_angle,list_of_torsion_angles)
  #print("workingposition",workingposition,"; atom extend",current_atom, direction, attach_atom, xcoord, ycoord, zcoord)
  #print(current_atom, attach_atom, direction, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle)
  #atom.extend([current_atom, attach_atom, direction, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle])
  atom.extend([current_atom, attach_atom, 0, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle])


#################################################################
# Check Chiral Centres
#print(list_stereogenic_centres)
#print("list_numbers_stereogenic_centres",list_numbers_stereogenic_centres)
if number_stereogenic_centres>0:
  #print("Stereocentres to check")
  #printatoms()
  sca=array.array("i")
  for s_centre in range(0, len(list_numbers_stereogenic_centres)):
    del sca
    sca=array.array("i")
    #print("stereocentre: ",list_numbers_stereogenic_centres[s_centre], s_centre)
    #print("list_numbers_stereogenic_centres",list_numbers_stereogenic_centres[s_centre], list_stereogenic_centres[s_centre])
    #printatoms()
    sca.append(atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],1))
    smallest = sca[-1]
    largest = sca[-1]
    #for workingposition in range(0, number_of_carbons):
    for workingposition in range(1, number_of_carbons+1):
      #print("growing sca",list_numbers_stereogenic_centres[s_centre],workingposition,atominfo(inchi2wn[workingposition],1),sca)
      if atominfo(inchi2wn[workingposition],1)==list_numbers_stereogenic_centres[s_centre]:
        sca.append(workingposition)
        if sca[-1] < smallest:
          smallest = sca[-1]
        if sca[-1] > largest:
          largest=sca[-1]
      #print("largest",workingposition,largest,sca,"lnsc",list_numbers_stereogenic_centres)
    second_smallest = largest
    second_largest = smallest
    for workingposition in range(0, len(sca)):
      if sca[workingposition] < second_smallest and not sca[workingposition] == smallest:
        second_smallest = sca[workingposition]
      if sca[workingposition] > second_largest and not sca[workingposition] == largest:
        second_largest = sca[workingposition]
    #print("sca, smallest, second_smallest, second_largest, largest",sca, smallest, second_smallest, second_largest, largest)
    if len(sca) > 3:
      smallest = second_smallest

    #print("stereocentre",atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],5),atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],6),atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],7))
    #print("smallest",atominfo(inchi2wn[smallest],5),atominfo(inchi2wn[smallest],6),atominfo(inchi2wn[smallest],7))
    #print("smallest",atominfo(inchi2wn[second_largest],5),atominfo(inchi2wn[second_largest],6),atominfo(inchi2wn[second_largest],7))
    #print("largest",atominfo(inchi2wn[largest],5),atominfo(inchi2wn[largest],6),atominfo(inchi2wn[largest],7))

    # Transfer three largest coordinates to matrix: smallest, second largest, largest
    m=array.array("i")
    m.append(atominfo(inchi2wn[largest],5)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],5))
    m.append(atominfo(inchi2wn[largest],6)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],6))
    m.append(atominfo(inchi2wn[largest],7)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],7))

    m.append(atominfo(inchi2wn[second_largest],5)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],5))
    m.append(atominfo(inchi2wn[second_largest],6)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],6))
    m.append(atominfo(inchi2wn[second_largest],7)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],7))

    m.append(atominfo(inchi2wn[smallest],5)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],5))
    m.append(atominfo(inchi2wn[smallest],6)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],6))
    m.append(atominfo(inchi2wn[smallest],7)-atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],7))

    #print("m",m)

    # We now have three points, translated so stereogenic centre is at (0,0,0)
    # These are the points with the larges InChI#
    # Construct a vector from the origin, to the average position
    # Construct a vector from the highest priority to the second highest priority
    # (priority increases with InChI#)
    # Cross these vectors together, than dot with vector to smallest point
    # Sign of the outcome should correlate with InChI +/- stereochemical indicator

    m2=array.array("i")
    m2.append( m[6] )
    m2.append( m[7] )
    m2.append( m[8] )

    m2.append( m[0]+m[3]+m[6] )
    m2.append( m[1]+m[4]+m[7] )
    m2.append( m[2]+m[5]+m[8] )

    m2.append( m[3]-m[0] )
    m2.append( m[4]-m[1] )
    m2.append( m[5]-m[2] )

    sign_stereo = m2[0]*(m2[4]*m2[8]-m2[5]*m2[7])-m2[1]*(m2[3]*m2[8]-m2[5]*m2[6])+m2[2]*(m2[3]*m2[7]-m2[4]*m2[6])

    #print("centre: ",list_numbers_stereogenic_centres[s_centre],sign_stereo,list_stereogenic_centres[s_centre])
    #printatoms()

    inchi_sign=1
    if list_stereogenic_centres[s_centre].find("+")==-1:
      inchi_sign=-1
    if inchi_sign*sign_stereo > 0:
      # need to invert centre
      #print("Inverting centre",list_stereogenic_centres[s_centre],sca)
      # If the stereocenter is a methine, may need to alter just one twist
      for workingposition in range(0, len(sca)):
        i = (3-atominfo(inchi2wn[sca[workingposition]],3))%3
        if i>0:
          setatominfo(inchi2wn[sca[workingposition]], 3, i)
          #print("check setatominfo",i,inchi2wn[sca[workingposition]],(3-atom[inchi2wn[sca[workingposition]]*8+3])%3,atominfo(inchi2wn[sca[workingposition]],3))
          #atom[inchi2wn[sca[workingposition]]*9+3]=(3-atom[inchi2wn[sca[workingposition]]*9+3])%3
          twist=atominfo(inchi2wn[sca[workingposition]],3)
          # why use direction, not newdirection?
          # could this be replaced by atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4) ??
          #direction=atominfo(inchi2wn[sca[workingposition]],2)
          #if direction == atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4):
          #  print("OK!", workingposition, direction)
          #else:
          #  print("Oh dear!", workingposition, direction, atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4))
          # direction (penultimate direction)
          # old_direction (direction of attached atom)
          direction=atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4)
          old_direction=atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],4)
          newdirection=twisting[twist*64+old_direction*8+direction]
          #print("about to change newdirection",sca[workingposition],inchi2wn[sca[workingposition]],"o,d,n",old_direction,direction,newdirection)
          setatominfo(inchi2wn[sca[workingposition]],4,newdirection)
    #printatoms()

  #Coordinates are now inconsistent with twists; and with conformation[]
  # This should set it right - but perhaps not needed at this stage
  #printatoms()
  for workingposition in range(current_position, number_of_carbons):
    attach_atom=atominfo(workingposition,1)
    newdirection=atominfo(workingposition,4)
    #print("prep ",current_position,workingposition,attach_atom,newdirection)
    xcoord = atominfo(inchi2wn[attach_atom],5)+directions[newdirection*3]
    ycoord = atominfo(inchi2wn[attach_atom],6)+directions[newdirection*3+1]
    zcoord = atominfo(inchi2wn[attach_atom],7)+directions[newdirection*3+2]
    setatominfo(workingposition,5,xcoord)
    setatominfo(workingposition,6,ycoord)
    setatominfo(workingposition,7,zcoord)
    #print("reset",workingposition,";",xcoord,ycoord,zcoord,";",attach_atom,newdirection)


#################################################################
# Everything should now be OK for connectivity, etc
# work out energy and geometry for different conformations
# conf_number=0
#################################################################

energy_list=[]
lowest_energy=999999
lowest_energy_conformation=0
if len(sys.argv) > 2:
  conf_num=conf_number
  energy=find_all_conf_number()
  if energy < 0:
    print("Energy of structure",sys.argv[2],"is > 200")
  else:
    print("Energy of structure",sys.argv[2],"is",energy)
else:
  # need to do conformation search
  conf_num=0
  while conf_num < number_conformations:
    tmp_conf_number=conf_num
    num_rot=0
    for i in range(0,number_rotatable_bonds):
      #print("conformation",number_rotatable_bonds,i,number_rotatable_bonds-i,conformation)
      conformation[number_rotatable_bonds-i-1]=tmp_conf_number%3
      tmp_conf_number=int(tmp_conf_number/3)
    #print("Conformation",conf_num,"Energy: ",find_all_conf_number(),conformation)
    energy=find_all_conf_number()
    if energy < 0:
      stop_atom=-1-energy
      # Should now be able to skip all structures with same motif
      # need to skip by torsion number not atom number...
      # What is relevant torsion angle from atom number?
      # atominfo(stop_atom,8)
      # comment out following lines to stop skipping:
      skip_number=int(3**(number_rotatable_bonds-atominfo(stop_atom,8)-1))-1
      conf_num+=skip_number
    else:
      if energy<lowest_energy:
        lowest_energy=energy
        lowest_energy_conformation=conf_num
      if energy < 200:
        energy_list.append([conf_num,energy])
        #print("Conformation",conf_num,"Energy: ",energy,conformation)
    old_conformation=conformation[:]
    conf_num+=1

  #printatoms()
  print()
  print("Dihedral angles, using InChI numbering")
  current_dihedral=0
  for i in range(number_of_carbons):
    #print("i",i,atom[i*9])
    a1=atom[i*9]
    a2=0
    a3=0
    a4=0
    if atom[inchi2wn[a1]*9+1]>0:
      a2=atom[inchi2wn[a1]*9+1]
      if atom[inchi2wn[a2]*9+1]>0:
        a3=atom[inchi2wn[a2]*9+1]
        if atom[inchi2wn[a3]*9+1]>0:
          a4=atom[inchi2wn[a3]*9+1]
          if atom[inchi2wn[a1]*9+8]==current_dihedral:
            current_dihedral +=1
            print("Dihedral",current_dihedral,"atoms:",a1,a2,a3,a4)
    #print(i,"Dihedral tmp",dihedral_number,"atoms:",a1,a2,a3,a4)
  print()
  print("Number of conformations with accessible energies (< 200):",len(energy_list))
  for i in range(0,len(energy_list)):
    tmp_conf_number=energy_list[i][0]
    #print("tmp_conf_number",energy_list[i],energy_list[i][0],tmp_conf_number)
    float_conformation=[]
    for j in range(0,number_rotatable_bonds):
      #print("conformation",number_rotatable_bonds,i,number_rotatable_bonds-i,conformation)
      #print("tmp_conf_number latest",tmp_conf_number,tmp_conf_number%3,120.0*(tmp_conf_number%3))
      float_conformation.append(180.0-120.0*(tmp_conf_number%3))
      #print("float_conformation",float_conformation)
      tmp_conf_number=int(tmp_conf_number/3)
    print("Conformation, Energy",energy_list[i],float_conformation)
  conf_num=lowest_energy_conformation
  tmp_conf_number=conf_num
  for i in range(0,number_rotatable_bonds):
    conformation[number_rotatable_bonds-i-1]=tmp_conf_number%3
    tmp_conf_number=int(tmp_conf_number/3)
  energy=find_all_conf_number()
  print()
  print("lowest energy",lowest_energy," at conformation",lowest_energy_conformation)


print()
printatoms()
print()
printxyzatoms()
#print(len(energy_list))
#print()
#printmolatoms()
#print("Number of conformations with accessible energies (< 200):",len(energy_list))
print()
end = time.process_time()
print("runningtime:", str(end-start))
