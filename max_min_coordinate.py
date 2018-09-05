from collections import defaultdict
import math

def max_min(coordinate, atom_type1, atom_type2, first_frame, last_frame, step):
    """Getting charges"""
    try:
        file = open("xmolout","r+")
        number_of_atoms = int(file.readline().split()[0])
        line = file.readline().split()
        file.close()
        statistics = int(math.ceil((last_frame-first_frame)/step))
        stat = 0
        maximum_x = [-999999999 for x in range(statistics)]
        maximum_y = [-999999999 for x in range(statistics)]
        maximum_z = [-999999999 for x in range(statistics)]
        minimum_x = [999999999 for x in range(statistics)]
        minimum_y = [999999999 for x in range(statistics)]
        minimum_z = [999999999 for x in range(statistics)]
        MAX_X = 0
        MAX_Y = 0
        MAX_Z = 0
        MIN_X = 0
        MIN_Y = 0
        MIN_Z = 0 
        
        
        for k in range(first_frame, last_frame, step):
            for j in range (0, number_of_atoms):
                if (coordinate[k][j][0] == atom_type1) or (coordinate[k][j][0] == atom_type2):
                    if (float(coordinate[k][j][1]) > maximum_x[stat]):
                        maximum_x[stat] = float(coordinate[k][j][1])
                    if (float(coordinate[k][j][1]) < minimum_x[stat]):
                        minimum_x[stat] = float(coordinate[k][j][1])
                    if (float(coordinate[k][j][2]) > maximum_y[stat]):
                        maximum_y[stat] = float(coordinate[k][j][2])
                    if (float(coordinate[k][j][2]) < minimum_y[stat]):
                        minimum_y[stat] = float(coordinate[k][j][2])
                    if (float(coordinate[k][j][3]) > maximum_z[stat]):
                        maximum_z[stat] = float(coordinate[k][j][3])
                    if (float(coordinate[k][j][3]) < minimum_z[stat]):
                        minimum_z[stat] = float(coordinate[k][j][3])      
                          
            MAX_X = MAX_X + maximum_x[stat]/statistics
            MAX_Y = MAX_Y + maximum_y[stat]/statistics
            MAX_Z = MAX_Z + maximum_z[stat]/statistics
            MIN_X = MIN_X + minimum_x[stat]/statistics
            MIN_Y = MIN_Y + minimum_y[stat]/statistics
            MIN_Z = MIN_Z + minimum_z[stat]/statistics
            stat = stat + 1
            
        print('Calculation details:'+'\n')
        print('COORDINATES: ')
        print('Maximum x coordinate: '+ str(MAX_X))
        print('Minimum x coordinate: '+ str(MIN_X))
        print('Maximum y coordinate: '+ str(MAX_Y))
        print('Minimum y coordinate: '+ str(MIN_Y))
        print('Maximum z coordinate: '+ str(MAX_Z))
        print('Minimum z coordinate: '+ str(MIN_Z))
        
        print('Average Volume: '+ str((MAX_X-MIN_X)*(MAX_Y-MIN_Y)*(MAX_Z-MIN_Z)))
       
        
        
    except IOError:
        pass
    
    return 0
