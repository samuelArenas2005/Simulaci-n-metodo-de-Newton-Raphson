from vector import Vector
import numpy as np


def main():
    xinit = Vector(729)
    xinit.vecInicialUniform()
    
    max_iterations = 50
    tolerance = 2
    
    print("Iniciando método de Newton-Raphson...")
    
    for i in range(max_iterations):
        print(f"Iteración {i + 1}:")
        
        # Calcular función y verificar convergencia
        xinit.cal_function()
        residuo_norm = np.linalg.norm(xinit.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        
        if residuo_norm < tolerance:
            print(f"¡Convergencia alcanzada en {i + 1} iteraciones!")
            break
            
        xinit.cal_jacobiano()
        xinit.cal_inv_jacobiano()
        xinit.newVector()
    else:
        print(f"No convergió después de {max_iterations} iteraciones")
    
    xinit.showPlot(True)


if __name__ == "__main__":
    main()



