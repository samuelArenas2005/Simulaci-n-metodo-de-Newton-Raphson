from vector import Vector


def main():
    xinit= Vector(729)
    xinit.vecInicial()
    
    for i in range(0, 8):
        xinit.cal_jacobiano()
        xinit.cal_inv_jacobiano()
        xinit.cal_function()
        xinit.newVector()
    
    xinit.showPlot(True)
        
    
   
    

if __name__ == "__main__":
    main()



