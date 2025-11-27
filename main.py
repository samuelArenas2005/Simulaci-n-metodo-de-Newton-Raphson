from vector import NewtonRaphson
from plotVec import Plot
from interfaz import SimulationViewer, view_heatmap, view_spline_static, view_anim_history, view_anim_spline

import numpy as np

#Método de configuración que permite expandir el número de nodos 
def get_scaled_config(factor):
    """
    Recibe un factor de escalado (ej. 1, 3, 5.5).
    Retorna la configuración (Nx, Ny, h, blocks) escalada.
    """
    # 1. Configuración Base (Original)
    Nx_base = 81
    Ny_base = 9
    h_base = 5.0
    blocks_base = [[10, 20, 0, 2], [58, 80, 5, 8]]

    # 2. Calcular nuevas dimensiones
    Nx_new = int(Nx_base * factor)
    Ny_new = int(Ny_base * factor)
    
    # El paso h disminuye inversamente al factor para aumentar precisión
    h_new = h_base / factor

    # 3. Escalar Bloques con Corrección de Techo
    blocks_new = []
    
    for bloque in blocks_base:
        x1, x2, y1, y2 = bloque
        
        # Escalado simple inicial
        nx1 = int(x1 * factor)
        nx2 = int(x2 * factor)
        ny1 = int(y1 * factor)
        ny2 = int(y2 * factor)
        
        if y2 == (Ny_base - 1):
            ny2 = Ny_new - 1
            
        if x2 == (Nx_base - 1):
            nx2 = Nx_new - 1
            
        blocks_new.append([nx1, nx2, ny1, ny2])

    return Nx_new, Ny_new, h_new, blocks_new

#Metodo de iteraciones de Newton Raphson AQUI ESTA TODA LA LOGICA DEL ALGORITMO PROGRAMADA
def iteration(NewtonRaphson,callNewVector,saveHistory=False, proves=False):
    
    #Condiciones de parada
    max_iterations = 200
    epsilon = 1e-13
    delta = 1e-13
    
    if saveHistory : NewtonRaphson.vecInicial()
    historial = [] # Aquí se guardan las "fotos" de la simulación
    numeros_condicion = []
    # Guardamos el estado inicial (Iteración 0)
    historial.append(NewtonRaphson.vec.copy())
    
    print("\nIniciando método de Newton-Raphson con el método del", callNewVector )
    print(f"Tolerancia residuo (epsilon): {epsilon}")
    print(f"Tolerancia cambio entre iteraciones (delta): {delta}")
    
     # Guardar el vector de la iteración anterior para comparar
    vec_anterior = None
    
    for i in range(max_iterations):
        print("Iteracion:", (i+1))
        
        # Guardar vector actual antes de modificarlo
        if vec_anterior is not None:
            # Calcular diferencia entre iteraciones
            cambio_norm = np.linalg.norm(NewtonRaphson.vec - vec_anterior)
            print(f"  Cambio entre iteraciones: {cambio_norm}")
            
            # Verificar condición de parada por delta
            if cambio_norm < delta:
                print(f"¡Convergencia por cambio mínimo alcanzada en {i + 1} iteraciones!")
                break
        
        # Calcular función y verificar convergencia por residuo
        NewtonRaphson.cal_function()
        residuo_norm = np.linalg.norm(NewtonRaphson.vectFunction)
        print(f"  Norma del residuo: {residuo_norm}")
        
        if residuo_norm < epsilon:
            print(f"\n¡Convergencia por residuo alcanzada en {i + 1} iteraciones!")
            break
        
        # Guardar vector actual antes de la actualización
        vec_anterior = NewtonRaphson.vec.copy()
        NewtonRaphson.cal_jacobiano()  

        # Comprobar propiedades de la matriz
        if proves:
            # Calcular e imprimir el número de condición de la Jacobiana en esta iteración
            numero_cond = NewtonRaphson.cal_condition_number()
            numeros_condicion.append(numero_cond)
            print(f"  Número de condición del Jacobiano: {numero_cond}")
            print(f"  Cumple Richardson: {NewtonRaphson.is_Richardson()}")
            print(f"  Es simétrica: {NewtonRaphson.is_symmetrical()}")
            print(f"  Es diagonalmente dominante: {NewtonRaphson.is_d_dominant()}")
            print(f"  Es de rango completo: {NewtonRaphson.is_full_rank()}")

        (NewtonRaphson.newVector() if callNewVector == "gradiente Conjugado" else NewtonRaphson.newVectorInversa())
        historial.append(NewtonRaphson.vec.copy())
    else:
        print(f"No convergió después de {max_iterations} iteraciones")

    if proves and numeros_condicion:
        promedio_cond = np.mean(numeros_condicion)
        print(f"\nNúmero de condición promedio: {promedio_cond}")
    
    print("----------------------------------------------")
    
    return historial
    
#Metodo principal que se encarga de ejecutar las iteraciones de Newton Rhapson
def main():
    FACTORA = 1
    FACTORB = 3   #Hacer el factor más pequeño hace que la prueba B corra más rapido 
    nxA, nyA, h_valA, bloqA = get_scaled_config(FACTORA)
    nxB, nyB, h_valB, bloqB = get_scaled_config(FACTORB)
    
    #Instancia de un vector de pruebas para ver el vector inicial
    NewtonRaphsonItGCPrueba= NewtonRaphson( Nx=nxA,  Ny=nyA,  blocks=bloqA,  V0=1,  v_ij=0.01,  h=h_valA, v0ParedSuperior= 1)
    NewtonRaphsonItGCPrueba.vecInicial()
    
    
    
    # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 1 Y VELOCIDAD 1 EN LA PARED SUPERIOR
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCA = NewtonRaphson( Nx=nxA,  Ny=nyA,  blocks=bloqA,  V0=1,  v_ij=0.01,  h=h_valA, v0ParedSuperior= 1)
    NewtonRaphsonItGCA.vecInicial()

    #Instancia del vector que sera calculado por la inversa del jacobiano
    NewtonRapshonItJIA =  NewtonRaphson( Nx=nxA,  Ny=nyA,  blocks=bloqA,  V0=1,  v_ij=0.01,  h=h_valA, v0ParedSuperior= 1)
    NewtonRapshonItJIA.vecInicial()
    
    # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 1 Y VELOCIDAD 0 EN LA PARED SUPERIOR 
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCA2 = NewtonRaphson( Nx=nxA,  Ny=nyA,  blocks=bloqA,  V0=1,  v_ij=0.01,  h=h_valA, v0ParedSuperior= 0)
    NewtonRaphsonItGCA2.vecInicial()
    
    # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 1, VELOCIDAD 0 EN LA PARED SUPERIOR Y VELOCIDAD INICIAL 2
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCA3 = NewtonRaphson( Nx=nxA,  Ny=nyA,  blocks=bloqA,  V0=2,  v_ij=0.05,  h=h_valA, v0ParedSuperior= 0)
    NewtonRaphsonItGCA3.vecInicial()
    
    # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 3 Y VELOCIDAD 1 EN LA PARED SUPERIOR
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCB = NewtonRaphson( Nx=nxB,  Ny=nyB,  blocks=bloqB,  V0=1,  v_ij=0.01,  h=h_valB, v0ParedSuperior= 1)
    NewtonRaphsonItGCB.vecInicial()

    # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 3 Y VELOCIDAD 0 EN LA PARED SUPERIOR
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCB2 = NewtonRaphson( Nx=nxB,  Ny=nyB,  blocks=bloqB,  V0=1,  v_ij=0.01,  h=h_valB, v0ParedSuperior= 0)
    NewtonRaphsonItGCB2.vecInicial()
    
     # ---------------------------------------------------------------------------
    # INSTANCIAS SIMULACIÓN CON FACTOR: 3 Y VELOCIDAD 0 EN LA PARED SUPERIOR Y VELOCIDAD INICIAL 2
    # ---------------------------------------------------------------------------
    #Instancia del vector que sera calculado por el gradiente conjugado
    NewtonRaphsonItGCB3 = NewtonRaphson( Nx=nxB,  Ny=nyB,  blocks=bloqB,  V0=2,  v_ij=0.05,  h=h_valB, v0ParedSuperior= 0)
    NewtonRaphsonItGCB3.vecInicial()


    # ---------------------------------------------------------------------------
    # COMPARACIÓN USO DE JACOBIANO INVERSO Y GRADIENTE CONJUGADO
    # ---------------------------------------------------------------------------
    
    #Ciclo principal del método calculado por el gradiente conjugado
    iteration(NewtonRaphsonItGCA,"gradiente Conjugado")
    
    #Ciclo principal del método calculado por la inversa del jacobiano 
    iteration(NewtonRapshonItJIA,"jacobiano Inverso")
    
    print("\nComparacion de ambos métodos")
    
    residuo = np.linalg.norm(NewtonRaphsonItGCA.vec - NewtonRapshonItJIA.vec)
    
    print(residuo)
    
    
    # ---------------------------------------------------------------------------
    # VISUALIZADOR UNIFICADO
    # ---------------------------------------------------------------------------
    print("Iniciando Visor Interactivo...")
    
    viewer = SimulationViewer()
    
    # ---------------------------------------------------------------------------
    # GRAFICAS SIMULACION A
    # ---------------------------------------------------------------------------
    def generateSimulationA(show):
        if not show : return 
        iteration(NewtonRaphsonItGCA2,"gradiente Conjugado")
        iteration(NewtonRaphsonItGCA3,"gradiente Conjugado")
        
        historial_pasosA =  iteration(NewtonRaphsonItGCA, "gradiente Conjugado",True)
        
        #Instancias PLOT
        plotNewtonRaphsonItGCA = Plot(NewtonRaphsonItGCA)
        plotNewtonRaphsonItGCA2 = Plot(NewtonRaphsonItGCA2)
        plotNewtonRaphsonItGCA3 = Plot(NewtonRaphsonItGCA3)
        
        # --- GRUPO A1 ---
        viewer.add_slide("A1: Mapa de Calor", view_heatmap, plotNewtonRaphsonItGCA)
        viewer.add_slide("A1: Historial Convergencia", view_anim_history, plotNewtonRaphsonItGCA, historial_pasosA, 200)
        viewer.add_slide("A1: Spline Estático (x20)", view_spline_static, plotNewtonRaphsonItGCA, 20)
        viewer.add_slide("A1: Animación Resolución Spline", view_anim_spline, plotNewtonRaphsonItGCA, 10, 800)
        
        # --- GRUPO A2 ---
        viewer.add_slide("A2: Mapa de Calor", view_heatmap, plotNewtonRaphsonItGCA2)
        viewer.add_slide("A2: Spline Estático (x20)", view_spline_static, plotNewtonRaphsonItGCA2, 20)
        
        # --- GRUPO A3 ---
        viewer.add_slide("A3: Mapa de Calor V0 = 2", view_heatmap, plotNewtonRaphsonItGCA3)
        viewer.add_slide("A3: Spline Estático (x20) V0 = 2", view_spline_static, plotNewtonRaphsonItGCA3, 20)
        
        def showAllGraph():
            #GRAFICAS A1
                #Mapa de calor
            plotNewtonRaphsonItGCA.showPlot(True)
                #Animacion Mapa de calor
            plotNewtonRaphsonItGCA.animateHistoryLoop(historial_pasosA,800)
                #Spline
            plotNewtonRaphsonItGCA.showPlotSpline(20)
                #Animacion spline
            plotNewtonRaphsonItGCA.animateSpline(800)
            
            #GRAFICAS A2
            #Mapa de calor
            plotNewtonRaphsonItGCA2.showPlot(True)
                #Spline
            plotNewtonRaphsonItGCA2.showPlotSpline(20)
            
            #GRAFICAS A3
            #Mapa de calor
            plotNewtonRaphsonItGCA3.showPlot(True)
                #Spline
            plotNewtonRaphsonItGCA3.showPlotSpline(20)

        #showAllGraph()           #Descomentarlo hace que se vean las graficas de forma individual una por una
    # ---------------------------------------------------------------------------
    # GRAFICAS SIMULACION B
    # ---------------------------------------------------------------------------
    
    def generateSimulationB(show,metodo):
        if not show : return 
        iteration(NewtonRaphsonItGCB,metodo)
        iteration(NewtonRaphsonItGCB2,metodo)
        iteration(NewtonRaphsonItGCB3,metodo)
        
        historial_pasosB = iteration(NewtonRaphsonItGCB, metodo,True)
        
        #Instancias PLOT
        plotNewtonRaphsonItGCB = Plot(NewtonRaphsonItGCB)
        plotNewtonRaphsonItGCB2 = Plot(NewtonRaphsonItGCB2)
        plotNewtonRaphsonItGCB3 = Plot(NewtonRaphsonItGCB3)
        
        # --- GRUPO B1 (mayor escala)  ---
        viewer.add_slide("B1: Mapa de Calor", view_heatmap, plotNewtonRaphsonItGCB)
        viewer.add_slide("B1: Historial Convergencia", view_anim_history, plotNewtonRaphsonItGCB, historial_pasosB, 200)
        viewer.add_slide("B1: Spline Estático (x20)", view_spline_static, plotNewtonRaphsonItGCB, 20)
        viewer.add_slide("B1: Animación Resolución Spline", view_anim_spline, plotNewtonRaphsonItGCB, 10, 800)
        
        # --- GRUPO B2 ---
        viewer.add_slide("B2: Mapa de Calor", view_heatmap, plotNewtonRaphsonItGCB2)
        viewer.add_slide("B2: Spline Estático (x20)", view_spline_static, plotNewtonRaphsonItGCB2, 20)
        
        # --- GRUPO B3 ---
        viewer.add_slide("B3: Mapa de Calor V0 = 2", view_heatmap, plotNewtonRaphsonItGCB3)
        viewer.add_slide("B3: Spline Estático (x20) V0 = 2", view_spline_static, plotNewtonRaphsonItGCB3, 20)
        
        def showAllGraph():
            #GRAFICAS B1
                #Mapa de calor
            plotNewtonRaphsonItGCB.showPlot(True)
                #Animacion Mapa de calor
            plotNewtonRaphsonItGCB.animateHistoryLoop(historial_pasosB,800)
                #Spline
            plotNewtonRaphsonItGCB.showPlotSpline(20)
                #Animacion spline
            plotNewtonRaphsonItGCB.animateSpline(800)
            
            #GRAFICAS A2
            #Mapa de calor
            plotNewtonRaphsonItGCB2.showPlot(True)
                #Spline
            plotNewtonRaphsonItGCB2.showPlotSpline(20)
            
            #GRAFICAS A3
            #Mapa de calor
            plotNewtonRaphsonItGCB3.showPlot(True)
                #Spline
            plotNewtonRaphsonItGCB3.showPlotSpline(20)
            
        #showAllGraph()           #Descomentarlo hace que se vean las graficas de forma individual una por una    

    # --- GRAFICAS SIMULACION PRUEBA ---
    plotNewtonRaphsonItGCPrueba = Plot(NewtonRaphsonItGCPrueba)
    
    
    generateSimulationA(True)
    generateSimulationB(False,"JacobianoInversa") #PONER EN TRUE CORRE LA PRUEBA CON UN FACTOR DE 3, 
                                                 #LO CUAL PUEDE LLEVAR MUCHO TIEMPO EN MOSTRARSE
                                                 #Poner de segundo parametro "gradiente Conjugado" si se quiere ejecutar con ese metodo 
    
    
    viewer.add_slide("VALORES INICIALES", view_heatmap, plotNewtonRaphsonItGCPrueba) #Añade la grafica de valores iniciales
    
    viewer.run()
    
main()

