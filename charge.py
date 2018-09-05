from collections import defaultdict
from charge_evaluater.filelength import filelen

def charge():
    """Getting charges"""
    try:
        file = open("fort.7","r+")
        number_of_atoms = int(file.readline().split()[0])
        file.close()

        frames = int(filelen("fort.7")/(number_of_atoms+3))
        charges = defaultdict(list)
        
        try:
            file = open("fort.7","r+")
            file.seek(0)
            
            
            for j in range(0,frames):
                
                null = file.readline().split()

                for k in range(0,number_of_atoms):
                    
                    line = file.readline().split()
                    charges[j+1].append([line[0], line[-1]])

                null = file.readline().split()
                null = file.readline().split()
                    
        finally:
            file.close()
    except IOError:
        pass
    
    return charges
