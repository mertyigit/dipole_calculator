from collections import defaultdict
import math

def dipole_jw_vector_boundaries(coordinate, charge, atom_type1, atom_type2, first_frame, last_frame, step, min_y, max_y, min_an, max_an, a, c):
    """Getting charges"""
    try:
        file = open("xmolout","r+")
        number_of_atoms = int(file.readline().split()[0])
        line = file.readline().split()
        file.close()
        statistics = int(math.ceil((last_frame-first_frame)/step))
        stat=0
        permitivity = 0
        dipole = [0 for x in range(3)]
        dipol = 0
        total_dipole = [0 for x in range(statistics)]
        total_dipole_x = [0 for x in range(statistics)]
        total_dipole_y = [0 for x in range(statistics)]
        total_dipole_z = [0 for x in range(statistics)]
        
        atom = 0
        molecule = 0
                
        for k in range(first_frame, last_frame, step):
            print('Frame number: '+str(k))
            dipole[0] = 0
            dipole[1] = 0
            dipole[2] = 0
            for j in range(0, number_of_atoms):
                if (coordinate[k][j][0] == atom_type1) or (coordinate[k][j][0] == atom_type2):
                    if float(charge[k][j][0]) > min_an:
                        if float(charge[k][j][0]) < max_an:
                            dipole[0] = dipole[0] + float(coordinate[k][j][1])*float(charge[k][j][1])*1.6e-29
                            dipole[1] = dipole[1] + float(coordinate[k][j][2])*float(charge[k][j][1])*1.6e-29
                            dipole[2] = dipole[2] + float(coordinate[k][j][3])*float(charge[k][j][1])*1.6e-29
                            atom = atom + 1
                
            total_dipole[stat] = (dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2])**0.5
            total_dipole_x[stat] = dipole[0]
            total_dipole_y[stat] = dipole[1]
            total_dipole_z[stat] = dipole[2]
            
            stat = stat + 1
            
        M2=0
        M=0
        M2_x=0
        M_x=0
        M2_y=0
        M_y=0
        M2_z=0
        M_z=0
        
        for i in range(0, statistics):
            M2 = M2 + total_dipole[i]*total_dipole[i]           
            M = M + total_dipole[i]
            M2_x = M2_x + total_dipole_x[i]*total_dipole_x[i]           
            M_x = M_x + total_dipole_x[i]
            M2_y = M2_y + total_dipole_y[i]*total_dipole_y[i]           
            M_y = M_y + total_dipole_y[i]
            M2_z = M2_z + total_dipole_z[i]*total_dipole_z[i]           
            M_z = M_z + total_dipole_z[i]
        b = max_y - min_y
        
        Volume = a*b*c*1e-30
        total_dipole_square_average = M2/statistics
        total_dipole_average_square = (M/statistics)*(M/statistics)
        total_dipole_square_average_x = M2_x/statistics
        total_dipole_average_square_x = (M_x/statistics)*(M_x/statistics) 
        total_dipole_square_average_y = M2_y/statistics
        total_dipole_average_square_y = (M_y/statistics)*(M_y/statistics)
        total_dipole_square_average_z = M2_z/statistics
        total_dipole_average_square_z = (M_z/statistics)*(M_z/statistics) 
        
        permitivity_new = (4*3.14*((total_dipole_average_square_x+total_dipole_average_square_y+total_dipole_average_square_z) - ( (M2_x/statistics)**2 + (M2_y/statistics)**2 + (M2_z/statistics)**2 ))/(3*(8.854187817e-12)*Volume*1.38e-23*300))+1
        permitivity = (total_dipole_square_average-total_dipole_average_square)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_x = (total_dipole_square_average_x-total_dipole_average_square_x)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_y = (total_dipole_square_average_y-total_dipole_average_square_y)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_z = (total_dipole_square_average_z-total_dipole_average_square_z)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        
        molecule = atom/statistics/3
        print('\n')
        print('Calculation details:'+'\n')
   
        print('POLARIZATION: ') 
        print ("Average dipole moment for one water molecule: "+str(M/(molecule)/statistics/3.34e-30))
        print ("Total dipole moment (Debye): "+str(M/statistics/3.34e-30)) 
        print ("Total dipole moment (x direction)(Debye): "+str(M_x/statistics/3.34e-30)) 
        print ("Total dipole moment (y direction)(Debye): "+str(M_y/statistics/3.34e-30)) 
        print ("Total dipole moment (z direction)(Debye): "+str(M_z/statistics/3.34e-30)) 
        
        print('Number of atoms in the system: '+ str(number_of_atoms))
        print('Number of atoms contributed to dipole moment calculation: '+str(atom/statistics))
        print('Number of molecules contributed to dipole moment calculation: '+str(molecule))
        print('Averaging over '+str(statistics)+' frames')
        print('Averaging over '+str(stat)+' frames')
        
        print('\n')
        print('Static dielectric permitivity (new): '+str(permitivity_new))
        print('Static dielectric permitivity: '+str(permitivity))
        print('Static dielectric permitivity (x direction): '+str(permitivity_x))
        print('Static dielectric permitivity (y direction): '+str(permitivity_y))
        print('Static dielectric permitivity (z direction): '+str(permitivity_z))
        print('deneme')
        
    except IOError:
        pass
    
    return permitivity
