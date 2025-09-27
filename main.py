from vector import Vector
from jacobiano import Jacobiano

def main():
    
    xinit= Vector(729)
    xinit.vecInicial()
    xinit.showPlot()
    
    jacobiano = Jacobiano(729)
    jacobiano.cal_jacobiano(xinit.vec)
    jacobiano.showMatrixJacobiana()

    

if __name__ == "__main__":
    main()



