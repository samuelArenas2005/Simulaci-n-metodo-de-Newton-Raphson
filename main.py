from vector import Vector
import numpy as np


def main():
    """ xinit = Vector(729)
    xinit.vecInicial()
    
    max_iterations = 200
    tolerance = 1e-7
    
    print("Iniciando método de Newton-Raphson...")
    
    for i in range(max_iterations):
        print("Iteracion:", (i+1))
        xinit.cal_function()
        residuo_norm = np.linalg.norm(xinit.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        if residuo_norm < tolerance:
            print(f"¡Convergencia alcanzada en {i + 1} iteraciones!")
            break
        xinit.cal_jacobiano()   
        xinit.newVector()
    else:
        print(f"No convergió después de {max_iterations} iteraciones")
    
    xinit.showPlot(True) """
    
    valores_derechos = [ecuacion * 81 - 1  for ecuacion in range(1, 10)]
    valores_bloque_A = [k * 81 + j for j in range(10,21) for k in range(0,3)]
    valores_bloque_B = [k * 81 + j for j in range(58,80) for k in range(5,9)]
    
    for ecuacion in range (729):
        if (ecuacion >= 81*8  and ecuacion<=728):
            print(ecuacion)


""" if __name__ == "__main__":
    main()

 """

