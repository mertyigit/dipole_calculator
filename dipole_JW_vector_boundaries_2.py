from collections import defaultdict
import math

def dipole_jw_vector_boundaries(coordinate, charge, atom_type1, atom_type2, first_frame, last_frame, step, min_y, max_y, a, c):
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
        atom = 0
        molecule = 0
                
        for k in range(first_frame, last_frame, step):
            dipole[0] = 0
            dipole[1] = 0
            dipole[2] = 0
            for j in range(0, number_of_atoms):
                if (coordinate[k][j][0] == atom_type1) or (coordinate[k][j][0] == atom_type2):
                    if float(coordinate[k][j][2]) > min_y:
                        if float(coordinate[k][j][2]) < max_y:
                            dipole[0] = dipole[0] + float(coordinate[k][j][1])*float(charge[k][j][1])*1.6e-29
                            dipole[1] = dipole[1] + float(coordinate[k][j][2])*float(charge[k][j][1])*1.6e-29
                            dipole[2] = dipole[2] + float(coordinate[k][j][3])*float(charge[k][j][1])*1.6e-29
                            atom = atom + 1
                
            total_dipole[stat] = (dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2])**0.5
            stat = stat + 1
            
        M2=0
        M=0
        for i in range(0, statistics):
            M2 = M2 + total_dipole[i]*total_dipole[i]           
            M = M + total_dipole[i]
        
        b = max_y - min_y
        
        Volume = a*b*c*1e-30
        total_dipole_square_average = M2/statistics
        total_dipole_average_square = (M/statistics)*(M/statistics)
        
        permitivity = (total_dipole_square_average-total_dipole_average_square)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
      
        molecule = atom/statistics/3
        print('\n')
        print('Calculation details:'+'\n')
   
        print('POLARIZATION: ') 
        print ("Average dipole moment for one water molecule: "+str(M/(molecule)/statistics/3.34e-30))
        print ("Total dipole moment (Debye): "+str(M/statistics/3.34e-30)) 
        print('Number of atoms in the system: '+ str(number_of_atoms))
        print('Number of atoms contributed to dipole moment calculation: '+str(atom/statistics))
        print('Number of molecules contributed to dipole moment calculation: '+str(molecule))
        print('Averaging over '+str(statistics)+' frames')
        print('Averaging over '+str(stat)+' frames')
        
        print('\n')
        print('Static dielectric permitivity: '+str(permitivity))
        
    except IOError:
        pass
    
    return permitivity
