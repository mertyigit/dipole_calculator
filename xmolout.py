from collections import defaultdict
from charge_evaluater.filelength import filelen
def coordinates():
    """Getting charges"""
    try:
        file = open("xmolout","r+")
        number_of_atoms = int(file.readline().split()[0])
        file.close()
        
        frames = int(filelen("xmolout")/(number_of_atoms+2))
        
        coordinates = defaultdict(list)
        
        try:
            file = open("xmolout","r+")
            file.seek(0)

            
            for j in range(1,frames):
                
                null = file.readline().split()
                null = file.readline().split()
                
                for k in range(0,number_of_atoms):
                    
                    line = file.readline().split()
                    coordinates[j].append([line[0], line[1], line[2], line[3]])
                    
                    
        finally:
            file.close()
    except IOError:
        pass
    
    return coordinates
