from vector import Vector
import numpy as np


def main():
    u = Vector(81, 9, [[10 ,20 ,0 ,2] , [58 ,80 ,5 ,9]], 1, 0.01, 5)
    u.vecInicial()
    
    max_iterations = 200
    tolerance = 1e-7
    
    print("Iniciando método de Newton-Raphson...")
    
    for i in range(max_iterations):
        print("Iteracion:", (i+1))
        u.cal_function()
        residuo_norm = np.linalg.norm(u.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        if residuo_norm < tolerance:
            print(f"¡Convergencia alcanzada en {i + 1} iteraciones!")
            break
        u.cal_jacobiano()   
        u.newVector()
    else:
        print(f"No convergió después de {max_iterations} iteraciones")
    
    u.showPlot(True)

main()


