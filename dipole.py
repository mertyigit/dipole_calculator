from collections import defaultdict
import math

def dipole(coordinate, charge, atom_type1, atom_type2, first_frame, last_frame, step, min_x, max_x, min_y, max_y, min_z, max_z):
    """Getting charges"""
    try:
        file = open("xmolout","r+")
        number_of_atoms = int(file.readline().split()[0])
        line = file.readline().split()
        file.close()
        statistics = int(math.ceil((last_frame-first_frame)/step))
        stat=0
        permitivity = 0
        dipole = [[0 for x in range(3)] for y in range(statistics)]
        total_dipole = [0 for x in range(statistics)]
        atom_number = 0
        
        for k in range(first_frame, last_frame, step):
            dipole[stat][0]=0
            dipole[stat][1]=0
            dipole[stat][2]=0
            for j in range(0, number_of_atoms):
                if (coordinate[k][j][0] == atom_type1) or (coordinate[k][j][0] == atom_type2):
                    if float(coordinate[k][j][1]) > min_x:
                        if float(coordinate[k][j][1]) < max_x:
                             if float(coordinate[k][j][2]) > min_y:
                                 if float(coordinate[k][j][2]) < max_y:
                                     if float(coordinate[k][j][3]) > min_z:
                                         if float(coordinate[k][j][3]) < max_z:
                                             dipole[stat][0] = dipole[stat][0] + float(coordinate[k][j][1])*float(charge[k][j][1])*1.6e-29
                                             dipole[stat][1] = dipole[stat][1] + float(coordinate[k][j][2])*float(charge[k][j][1])*1.6e-29
                                             dipole[stat][2] = dipole[stat][2] + float(coordinate[k][j][3])*float(charge[k][j][1])*1.6e-29
                                             atom_number = atom_number + 1
                                             
            total_dipole[stat] = ((dipole[stat][0]*dipole[stat][0]+dipole[stat][1]*dipole[stat][1]+dipole[stat][2]*dipole[stat][2])**(0.5))
            print(str(dipole[stat][0])+' '+str(dipole[stat][1])+' '+str(dipole[stat][2]))           
            stat = stat + 1
            
        M2=0
        M=0
        for i in range(0, statistics):
            M2 = M2 + total_dipole[i]*total_dipole[i]           
            M = M + total_dipole[i]
            print(total_dipole[i])
            print(i)
        Volume = (max_x-min_x)*(max_y-min_y)*(max_z-min_z)*1e-30
        total_dipole_square_average = M2/statistics
        total_dipole_average_square = (M/statistics)*(M/statistics)
        
        permitivity = 4*(3.14)*(total_dipole_square_average-total_dipole_average_square)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        
        print ("Average dipole moment for one water molecule: "+str(M/309/statistics/3.34e-30))
        print ("Total dipole moment (Debye): "+str(M/statistics/3.34e-30))
        print(total_dipole_square_average-total_dipole_average_square)
        print((3*8.854187817e-12*1.38e-23*300*Volume))
        print(Volume)
        print(total_dipole_square_average)
        print(total_dipole_average_square)
        
        
        print('Calculation details:'+'\n')
        print('Number of atoms in the system: '+ str(number_of_atoms))
        print('Number of atoms contributed to dipole moment calculation: '+str(atom_number/statistics)+'\n')
        print('Averaging over '+str(statistics)+' frames'+'\n')
        print('Averaging over '+str(stat)+' frames'+'\n')
        
        
    except IOError:
        pass
    
    return permitivity
