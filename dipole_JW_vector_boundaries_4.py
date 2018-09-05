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
        total_dipole_x = [0 for x in range(statistics)]
        total_dipole_y = [0 for x in range(statistics)]
        total_dipole_z = [0 for x in range(statistics)]
        
        water = [0 for x in range(statistics)]
        #hydrogen = 0
        #oxygen = 0
        #a = 0
        atom = 0
        molecule = 0
        for k in range(first_frame, last_frame, step):
            print('Frame number: '+str(k))
            mol = 0
            dipole[0] = 0
            dipole[1] = 0
            dipole[2] = 0
            for j in range(0, number_of_atoms):
                if (coordinate[k][j][0] == 'O'):
                    if float(coordinate[k][j][2]) > min_y:
                        if float(coordinate[k][j][2]) < max_y:
                            #oxygen = oxygen + 1
                            
                            O_x = float(coordinate[k][j][1])
                            O_y = float(coordinate[k][j][2])
                            O_z = float(coordinate[k][j][3])
                            O_c = float(charge[k][j][1])
                            O_t = coordinate[k][j][0]
                            noh = 0
                            for i in range(0, number_of_atoms):
                                H_x = 0
                                H_y = 0
                                H_z = 0
                                if (coordinate[k][i][0] == 'H') and (O_x - 5 < float(coordinate[k][i][1]) < O_x + 5) and (O_y - 5 < float(coordinate[k][i][2]) < O_y + 5) and (O_z - 5 < float(coordinate[k][i][3]) < O_z + 5):
                                    H_x = float(coordinate[k][i][1])
                                    H_y = float(coordinate[k][i][2])
                                    H_z = float(coordinate[k][i][3])
                                    H_t = coordinate[k][i][0]
                                    p = 0
                                    R = (((H_x-O_x)**2)+((H_y-O_y)**2)+((H_z-O_z)**2))**0.5

                                    if (R < 1.2) and (R > 0) and (noh == 0):
                                        
                                        Hydrogen1_x = H_x
                                        Hydrogen1_y = H_y                                
                                        Hydrogen1_z = H_z
                                        Hydrogen1_c = float(charge[k][i][1])
                                        Hydrogen1_t = H_t
                                        noh = noh + 1
                                        p = 1
                                        
                                        #print(str(charge[k][i][0])+' '+str(Hydrogen1_t))
                                        #hydrogen = hydrogen + 1
                                        
                                    if (R < 1.2) and (R > 0) and (noh == 1) and (p==0):
                                    
                                        Hydrogen2_x = H_x
                                        Hydrogen2_y = H_y                                
                                        Hydrogen2_z = H_z
                                        Hydrogen2_c = float(charge[k][i][1])
                                        Hydrogen2_t = H_t
                                        noh = noh + 1 
                                        #print(str(charge[k][i][0])+' '+str(Hydrogen2_t))
                                        #hydrogen = hydrogen + 1
                            if (noh == 2):
                                Oxygen_x = O_x
                                Oxygen_y = O_y
                                Oxygen_z = O_z
                                Oxygen_c = O_c
                                Oxygen_t = O_t
                                
                            if (noh == 2):
                                
                                dipole_x = (Hydrogen1_x*Hydrogen1_c*1.6e-29) + (Hydrogen2_x*Hydrogen2_c*1.6e-29) + (Oxygen_x*Oxygen_c*1.6e-29)
                                dipole_y = (Hydrogen1_y*Hydrogen1_c*1.6e-29) + (Hydrogen2_y*Hydrogen2_c*1.6e-29) + (Oxygen_y*Oxygen_c*1.6e-29)
                                dipole_z = (Hydrogen1_z*Hydrogen1_c*1.6e-29) + (Hydrogen2_z*Hydrogen2_c*1.6e-29) + (Oxygen_z*Oxygen_c*1.6e-29)
                                
                                dipole[0] = dipole[0] + dipole_x
                                dipole[1] = dipole[1] + dipole_y
                                dipole[2] = dipole[2] + dipole_z
                                water_dipole = (dipole_x**2+dipole_y**2+dipole_z**2)**0.5

                                #print(float(water_dipole)/3.34e-30)
                                #print(str(Hydrogen1_t)+' '+str(Hydrogen1_x)+' '+str(Hydrogen1_y)+' '+str(Hydrogen1_z)+' '+str(Hydrogen1_c))
                                #print(str(Hydrogen2_t)+' '+str(Hydrogen2_x)+' '+str(Hydrogen2_y)+' '+str(Hydrogen2_z)+' '+str(Hydrogen2_c))
                                #print(str(Oxygen_t)+' '+str(Oxygen_x)+' '+str(Oxygen_y)+' '+str(Oxygen_z)+' '+str(Oxygen_c))
                                mol = mol + 1
                                molecule = molecule + 1
                            #if ( noh == 1 ):
                             #   print ("noh is 1")
                              #  a = a + 1
                            if (noh > 2):
                                print("Bond length criteria is not correct!")
            #print(a)    
            total_dipole[stat] = (dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2])**0.5
            #print(total_dipole[stat])
            #print(hydrogen)
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
        molecule = molecule / statistics
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
        permitivity = (4*3.14*(total_dipole_square_average-total_dipole_average_square))/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_x = 4*3.14*(total_dipole_square_average_x-total_dipole_average_square_x)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_y = 4*3.14*(total_dipole_square_average_y-total_dipole_average_square_y)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
        permitivity_z = 4*3.14*(total_dipole_square_average_z-total_dipole_average_square_z)/(3*(8.854187817e-12)*(1.38e-23)*300*Volume)+1
      
        print('\n')
        print('Calculation details:'+'\n')
   
        print('POLARIZATION: ') 
        print ("Average dipole moment for one water molecule: "+str(M/(molecule)/statistics/3.34e-30))
        print ("Total dipole moment (Debye): "+str(M/statistics/3.34e-30)) 
        print ("Total dipole moment (x direction)(Debye): "+str(M_x/statistics/3.34e-30)) 
        print ("Total dipole moment (y direction)(Debye): "+str(M_y/statistics/3.34e-30)) 
        print ("Total dipole moment (z direction)(Debye): "+str(M_z/statistics/3.34e-30)) 
        
        print('Number of atoms in the system: '+ str(number_of_atoms))
        print('Number of molecules contributed to dipole moment calculation: '+str(molecule))
        print('Averaging over '+str(statistics)+' frames')
        print('Averaging over '+str(stat)+' frames')
        
        print('\n')
        print('Static dielectric permitivity (new): '+str(permitivity_new))
        print('Static dielectric permitivity: '+str(permitivity))
        print('Static dielectric permitivity (x direction): '+str(permitivity_x))
        print('Static dielectric permitivity (y direction): '+str(permitivity_y))
        print('Static dielectric permitivity (z direction): '+str(permitivity_z))
        
    except IOError:
        pass
    
    return permitivity
