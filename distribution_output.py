def distribution_output(distribution):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open("distribution","a+")        
        file.seek(0)
        
        try:
            for i in range(0, len(distribution)):
                for j in range(0, len(distribution[0])):
                    
                        file.write(str(round(distribution[i][j], 6)).rjust(9) + "  ") 
                file.write("\n")
            
        finally:
            file.close()
    except IOError:
        pass
