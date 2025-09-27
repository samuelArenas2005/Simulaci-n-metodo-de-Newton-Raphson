from vector import Vector
from jacobiano import Jacobiano

def main():
    
    xinit= Vector(729)
    xinit.vecInicial()
    xinit.showPlot()
    
    jacobiano = Jacobiano(729)
    jacobiano.cal_jacobiano(xinit.vec)
    
    #si quieres crear un nuevo archivo excel para ver la matriz pon show=True en los parametros xd
    jacobiano.showMatrixJacobiana() 
    
    jacobiano.cal_inv_jacobiano()
    jacobiano.showMatrixJacobiana()
    

if __name__ == "__main__":
    main()



