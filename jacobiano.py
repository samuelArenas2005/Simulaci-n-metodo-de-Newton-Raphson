import random
import numpy as np
import pandas as pd
from datetime import datetime

class Jacobiano:
    
    def __init__(self,n):
        self.matrixJacobiana = np.zeros((n,n))
    
    def cal_jacobiano(self,vec_0):
        
        v_ij= 0.2
        
        valores_derechos = [ecuacion * 81 - 1  for ecuacion in range(1, 10)]
        valores_bloque_A = [j * 81 + k for j in range(11,20) for k in range(0,3)]
        valores_bloque_B = [j * 81 + k for j in range(36,58) for k in range(5,8)]
        
        for ecuacion in range (729):
            if ((ecuacion % 81 == 0) or (ecuacion >= 0 and ecuacion<=81) or (ecuacion >= 81*8  and ecuacion<=728) or ecuacion in valores_derechos 
                    or ecuacion in valores_bloque_A or ecuacion in valores_bloque_B):
                
                self.matrixJacobiana[ecuacion,ecuacion] = 1
            
            else:
                #aqui van las derivadas parciales, en la posición [i,i] [i,i+1] [i,i-1].... Tal y como lo definimos
                i = ecuacion #Por comodidad visual 
                
                self.matrixJacobiana[ecuacion,i] = (-5/8)*(vec_0[i+1]-vec_0[i-1])
                
                self.matrixJacobiana[ecuacion,i+1] = (1/4) - (5/8)*vec_0[i]
                
                self.matrixJacobiana[ecuacion,i-1] = (1/4) + (5/8)*vec_0[i]
                
                self.matrixJacobiana[ecuacion,i+81] = (1/4) - (5/8)*v_ij
                
                self.matrixJacobiana[ecuacion,i-81] = (1/4) + (5/8)*v_ij
                
    def cal_inv_jacobiano(self):
        self.matrixJacobiana =  np.linalg.inv(self.matrixJacobiana)
                    
                    
    def showMatrixJacobiana(self, filename="jacobiano.xlsx", show=False):
        if show:
            try:
                mat = self.matrixJacobiana.toarray()
            except:
                mat = self.matrixJacobiana

        df = pd.DataFrame(mat)

        # timestamp para nombre único
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        new_filename = f"jacobiano_{timestamp}.xlsx"

        df.to_excel(new_filename, index=False, header=False)
        print(f"Jacobiano guardado en {new_filename}")

        
    