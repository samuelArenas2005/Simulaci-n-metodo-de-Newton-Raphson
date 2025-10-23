from vector import Vector
import numpy as np


def main():
    u = Vector(Nx = 81, Ny = 9, blocks = [[10 ,20 ,0 ,2] , [58 ,80 ,5 ,8]], V0 = 1, v_ij = 0.01, h = 5)
    u.vecInicial()
    uOld =  Vector(Nx = 81, Ny = 9, blocks = [[10 ,20 ,0 ,2] , [58 ,80 ,5 ,8]], V0 = 1, v_ij = 0.01, h = 5)
    uOld.vecInicial()
    
    max_iterations = 200
    epsilon = 1e-14
    delta = 1e-14
    
    print("Iniciando método de Newton-Raphson...")
    print(f"Tolerancia residuo (epsilon): {epsilon}")
    print(f"Tolerancia cambio entre iteraciones (delta): {delta}")
    
    # Guardar el vector de la iteración anterior para comparar
    vec_anterior = None
    vec_anteriorOld = None
    
    for i in range(max_iterations):
        print("Iteracion:", (i+1))
        
        # Guardar vector actual antes de modificarlo
        if vec_anterior is not None:
            # Calcular diferencia entre iteraciones
            cambio_norm = np.linalg.norm(u.vec - vec_anterior)
            print(f"  Cambio entre iteraciones: {cambio_norm}")
            
            # Verificar condición de parada por delta
            if cambio_norm < delta:
                print(f"¡Convergencia por cambio mínimo alcanzada en {i + 1} iteraciones!")
                break
        
        # Calcular función y verificar convergencia por residuo
        u.cal_function()
        residuo_norm = np.linalg.norm(u.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        
        if residuo_norm < epsilon:
            print(f"¡Convergencia por residuo alcanzada en {i + 1} iteraciones!")
            break
        
        # Guardar vector actual antes de la actualización
        vec_anterior = u.vec.copy()
        
        u.cal_jacobiano()  
        u.newVector()
    else:
        print(f"No convergió después de {max_iterations} iteraciones")
    

    for i in range(max_iterations):
        print("Iteracion:", (i+1))
        
        # Guardar vector actual antes de modificarlo
        if vec_anteriorOld is not None:
            # Calcular diferencia entre iteraciones
            cambio_norm = np.linalg.norm(uOld.vec - vec_anteriorOld)
            print(f"  Cambio entre iteraciones: {cambio_norm}")
            
            # Verificar condición de parada por delta
            if cambio_norm < delta:
                print(f"¡Convergencia por cambio mínimo alcanzada en {i + 1} iteraciones!")
                break
        
        # Calcular función y verificar convergencia por residuo
        uOld.cal_function()
        residuo_norm = np.linalg.norm(uOld.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        
        if residuo_norm < epsilon:
            print(f"¡Convergencia por residuo alcanzada en {i + 1} iteraciones!")
            break
        
        # Guardar vector actual antes de la actualización
        vec_anteriorOld = uOld.vec.copy()
        
        uOld.cal_jacobiano()  
        uOld.newVectorInversa()
    else:
        print(f"No convergió después de {max_iterations} iteraciones")
    
    u.showPlot(True)

main()