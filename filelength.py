from collections import defaultdict

def filelen(file_name):
    """Number of lines"""   
    try:
        
        with open(file_name) as file:
            for i, l in enumerate(file):
                pass
        return i+1
          
    finally:
            file.close()    
