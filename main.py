from vector import Vector


def main():
    
    xinit= Vector(729)
    xinit.vecInicial()
    xinit.showPlot(True)
    xinit.cal_jacobiano()
    xinit.cal_inv_jacobiano()
    xinit.getfunction(1)
    xinit.showInConsoleFunction(True)
    #comentario prueba
    #comentario prueba


   
    
    #si quieres crear un nuevo archivo excel para ver la matriz pon show=True en los parametros xd
 

if __name__ == "__main__":
    main()



