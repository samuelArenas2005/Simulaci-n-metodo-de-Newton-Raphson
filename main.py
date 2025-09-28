from vector import Vector


def main():
    xinit= Vector(729)
    xinit.vecInicial(20)
    for i in range(0, 2):
        xinit.cal_jacobiano()
        xinit.cal_inv_jacobiano()
        xinit.cal_function(20)
        xinit.newVector()
    xinit.showPlot(True)
    

if __name__ == "__main__":
    main()



