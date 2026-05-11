# Paso 1: librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from time import time

##################### apartado a #####################
# Paso 2: Parámetros generales, dominio y condiciones iniciales

# Dominio espacial
x_min = -1.0                        # límite inferior del dominio espacial
x_max = 1.0                         # límite superior del dominio
Nx = 101                            # número de puntos espaciales
x = np.linspace(x_min, x_max, Nx)   # malla (según el esquema debería ser el paso 3 pero lo prefiero aqui)
dx = x[1] - x[0]                    # tamaño del paso espacial

# Condiciones iniciales: T(x,0) = 100 exp(-20 x^2)
T0 = 100 * np.exp(-20 * x**2)

# Aseguramos Dirichlet T=0 en las fronteras desde el inicio
T0[0] = 0.0                         # borde izquierdo
T0[-1] = 0.0                        # borde derecho

# Tiempo de simulación total y paso temporal (ajustar para la estabilidad)
t_final = 1.0                       # tiempo total de simulación

# Paso 3: incializamos los arrays
D = 0.1                             # Coeficiente de difusión constante (apartado a, b, c)

# Paso 4: Paso temporal para método explícito FTCS cumple la condición de estabilidad:
dt_ftcs = dx**2 / (2*D) * 0.98      # 0.98 factor seguridad (por errores numéricos de python)
Nt_ftcs = int(t_final / dt_ftcs)    # Número de pasos temporales

# Paso 5: Inicialización de matriz temperatura y copia para evolución
T_ftcs = T0.copy()                  # hacemos copia para no modificar T0

# Paso 6: Función para aplicar condiciones de frontera de Dirichlet (T=0)
def aplicar_dirichlet(T):
    T[0] = 0.0                      # fijamos borde izquierdo
    T[-1] = 0.0                     # fijamos borde derecho
    return T                        # devolvemos T modificada

# función que realiza la evolución FTCS usando matrices sparse
def simulate_ftcs_matrix(D, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size             # número de puntos espaciales
    r = D * dt / dx**2              # coeficiente r = D*dt/dx^2

    # Paso 7: Matriz difusión (operador Laplaciano 1D)
    # Construir las diagonales del operador discreto L tal que L \cdot T ≈ (T_{j+1}-2T_j+T_{j-1})
    lower = np.ones(Nx-1)       # subdiagonal (son -1)
    main  = -2.0 * np.ones(Nx)  # diagonal principal
    upper = np.ones(Nx-1)       # superdiagonal (son +1)


    L = diags([lower, main, upper], offsets=[-1, 0, 1], shape=(Nx, Nx), format='csr')   # Construcción del operador discreto L

    # filas 0 y -1 deben ser 0 para no modificar frontera
    L = L.tolil()                # convertimos a formato LIL para modificar filas
    L[0,0] = 1                   # identidad en borde solo por estabilidad
    L[0, :] = 0.0                # fila 0 a cero
    L[-1, :] = 0.0               # fila final a cero
    L[-1,-1] = 1                 # identidad
    L = csr_matrix(L)            # volver a csr para velocidad

    T = T_initial.copy()         # inicializamos temperatura
    T = bc_func(T)               # aplicamos condiciones de contorno
    saved = []                   # creamos lista para ir guardando datos

    #Paso 10: Resolver
    start_time = time()
    for n in range(Nt):
        T_new = T + r * (L.dot(T))  # actualización explícita FTCS usando el operador L
        T_new = bc_func(T_new)      # aplicar condiciones de contorno Dirichlet
        T = T_new                   # asignar para siguiente paso
        if n in times_to_save:
            saved.append(T.copy())  # guardar si corresponde
    end_time = time()

    tiempo_total = end_time - start_time
    return saved, tiempo_total, T


times_to_save = [0, int(Nt_ftcs/6), int(Nt_ftcs/5), int(Nt_ftcs/4), int(Nt_ftcs/3), int(Nt_ftcs/2), Nt_ftcs-1]
saved_T_ftcs, tiempo_ftcs, T_final = simulate_ftcs_matrix(D=D, dx=dx, dt=dt_ftcs, Nt=Nt_ftcs, T_initial=T0, bc_func=aplicar_dirichlet, times_to_save=times_to_save)
print("Tiempo FTCS (matriz sparse):", tiempo_ftcs)

# Paso 11: Graficar 
plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_ftcs):
    t_real = times_to_save[i] * dt_ftcs
    plt.plot(x, T_snap, linewidth=2, alpha=0.9, label=f"t = {t_real:.3f}s")
plt.title("Difusión 1D — Método FTCS con Matrices", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("Temperatura", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()


##################### apartado c, vamos a repetir el a solo que ahora D depende de la posición (aunque no se verá afectado por la derivada) #####################
# Paso 6: Inicializamos D(x) como un array de Nx elementos (ejemplo: D aumenta con x^2)
alpha = 0.1                                # parámetro que controla la magnitud del coeficiente de difusión
D_x = alpha * (1 + np.exp(-x**2))          # definimos D(x) como una función dependiente de la posición

# Paso 7: Paso temporal para método explícito (FTCS)
dt_ftcs = dx**2 / (2*np.max(D_x)) * 0.98   # paso temporal estable usando el máximo D(x) (criterio FTCS)
Nt_ftcs = int(t_final / dt_ftcs)           # calculamos el número total de pasos de tiempo

T_ftcs = T0.copy()                         # copiamos la condición inicial para no modificar T0

# Paso 9: Función para condiciones Dirichlet
def aplicar_dirichlet(T):
    T[0] = 0.0                              # imponemos temperatura 0 en el borde izquierdo
    T[-1] = 0.0                             # imponemos temperatura 0 en el borde derecho
    return T                                # devolvemos el array modificado



# Paso 10: Función generalizada FTCS con D(x)
def simulate_ftcs_matrix_Dx(D_x, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size                     # número total de puntos espaciales
    T = T_initial.copy()                    # copiamos el campo inicial de temperatura
    T = bc_func(T)                          # aplicamos condiciones de borde iniciales
    saved = []                              # lista donde se almacenarán los valores

    
    # Inicializamos la matriz de Laplace como antes
    lower = np.ones(Nx-1)       # subdiagonal (coeficiente de T_{j-1})
    main  = -2.0 * np.ones(Nx)  # diagonal principal (coeficiente de T_j)
    upper = np.ones(Nx-1)       # superdiagonal (coeficiente de T_{j+1})
    
    L_in = diags([lower, main, upper], offsets=[-1, 0, 1], shape=(Nx, Nx), format='csr') # Construcción de la matriz sparse

    Dmat = diags(D_x, 0, shape=(Nx, Nx), format='csr')  #contruimos matriz cuyo elemento diagonal es D_x[j] (representamos D(x) en forma matricial)
    L = Dmat @ L_in                                     # hacemos el producto matricial
    
    # Fijar Dirichlet
    L = L.tolil()                            # convertimos a formato LIL para modificar filas
    L[0,:] = 0.0                             # anulamos la fila 0 para mantener la condición de frontera
    L[0,0] = 1.0                             # ponemos un 1 en la diagonal por estabilidad
    L[-1,:] = 0.0                            # anulamos la última fila
    L[-1,-1] = 1.0                           # ponemos un 1 en la diagonal en el borde derecho
    L = csr_matrix(L)                        # volvemos a CSR para cálculos rápidos


    #Paso 10: Resolver
    start_time = time()
    for n in range(Nt):
        T_new = T + dt/dx**2 * (L.dot(T))   # actualización explícita FTCS usando el operador L
        T_new = bc_func(T_new)              # aplicar condiciones de contorno Dirichlet
        T = T_new                           # asignar para siguiente paso
        if times_to_save is not None and n in times_to_save:
            saved.append(T.copy())          # guardar si corresponde
    end_time = time()
    tiempo_total = end_time - start_time
    return saved, tiempo_total, T


times_to_save = [0, int(Nt_ftcs/6), int(Nt_ftcs/5), int(Nt_ftcs/4), int(Nt_ftcs/3), int(Nt_ftcs/2), Nt_ftcs-1]
saved_T_ftcs, tiempo_ftcs, T_final_ftcs = simulate_ftcs_matrix_Dx(D_x, dx, dt_ftcs, Nt_ftcs, T0, aplicar_dirichlet, times_to_save)
print("Tiempo FTCS con D(x):", tiempo_ftcs)

plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_ftcs):
    t_real = times_to_save[i]*dt_ftcs
    plt.plot(x, T_snap, linewidth=2, alpha=0.8, label=f"t={t_real:.3f}")
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Difusión 1D con D(x) variable — FTCS")
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()









##################### apartado b #####################

# Parámetros y discretización temporal 
dt_cn = 0.005                           # paso temporal para Crank–Nicolson
Nt_cn = int(t_final / dt_cn)            # número de pasos temporales totales


# Paso 3: Función para aplicar condiciones de frontera mixtas
def aplicar_bc_cn(T,dx):
    T[0] = 0.0                          # Dirichlet en x=min: fijamos T en 0
    # Neumann modificada en x=1: dT/dx = -T
    # Aproximación de primer orden: (T_N - T_{N-1})/dx = -T_N => T_N*(1+dx) = T_{N-1} 
    T[-2] = T[-1]*(1 + dx)            # imponer la condición de frontera Neumann modificada
    return T                            # devolver el vector T con BC aplicadas



# Función Crank-Nicolson para D no dependiente de la posición
def simulate_cn_D(D, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size                  # número de puntos espaciales
    T = T_initial.copy()                 # copia del estado inicial para no modificar T_initial
    T = bc_func(T, dx)                   # aplicar condiciones de frontera iniciales
    saved = []                           # lista para almacenar la informacion


    # Construcción del Laplaciano estándar
    lower = np.ones(Nx-1)                # vector para la subdiagonal (coef. de T_{j-1})
    main  = -2*np.ones(Nx)               # vector para la diagonal principal (coef. de T_j)
    upper = np.ones(Nx-1)                # vector para la superdiagonal (coef. de T_{j+1})


    L = diags([lower, main, upper], [-1,0,1], shape=(Nx, Nx)) / dx**2 # dividimos por dx^2 para representar la derivada segunda
    I = diags([np.ones(Nx)], [0], shape=(Nx, Nx), format='csr')       # matriz identidad de tamaño Nx x Nx

    # A y B para Crank–Nicolson con D constante
    A = (I - 0.5 * dt * (D * L)).tolil()  # matriz izquierda: I - 0.5*dt*D*L, convertida a LIL para editar filas
    B = (I + 0.5 * dt * (D * L)).tolil()  # matriz derecha: I + 0.5*dt*D*L, convertida a LIL


    # Condiciones de contorno: frontera izquierda (Dirichlet)
    A[0, :] = 0                           # anular toda la fila 0 para imponer condición
    A[0, 0] = 1                           # poner 1 en la diagonal (T_0 = valor impuesto)
    B[0, :] = 0                           # anular fila 0 en B para consistencia
    B[0, 0] = 1                           # poner 1 en B[0,0] también

    # Condiciones de Neumann
    A[-1, :] = 0                          # anular toda la última fila antes de imponer la ecuación BC
    A[-1, -1] = 1 + dt                    # entrada A_{N,N} = 1 + dt  (en el informe se deriva la demo)
    A[-1, -2] = -1                        # entrada A_{N,N-1} = -1 para relacionar T_N con T_{N-1}
    B[-1, :] = 0                          # anulamos la última fila en B
    B[-1, -1] = 0                         # definimos B[-1,-1]=0 

    A = csr_matrix(A)                     # convertir A a formato CSR 
    B = csr_matrix(B)                     # convertir B a CSR


    # Bucle temporal de Crank–Nicolson
    start_time = time()                   # tiempo de inicio para medir cuanto tarda
    for n in range(Nt):                   # iterar Nt veces (cada iteración = un paso dt)
        b = B.dot(T)                      # construir vector lado derecho b = B * T^n
        T_new = spsolve(A, b)             # resolver el sistema A * T^{n+1} = b usando spsolve
        T_new = bc_func(T_new, dx)        # reimponer condiciones de frontera en la solución nueva
        T = T_new                         # actualizar T para el siguiente paso

        if times_to_save is not None and n in times_to_save:
            saved.append(T.copy())        # guardar snapshot si corresponde

    end_time = time()                     # tiempo final
    tiempo_total = end_time - start_time  # calcular tiempo de ejecución
    return saved, tiempo_total, T         # devolver snapshots, tiempo total y solución final

times_to_save_cn = [0, int(Nt_cn/6), int(Nt_cn/5), int(Nt_cn/4), int(Nt_cn/3), int(Nt_cn/2), Nt_cn-1]
saved_T_cn, tiempo_cn, T_final_cn = simulate_cn_D(D, dx, dt_cn, Nt_cn, T0, aplicar_bc_cn, times_to_save_cn)
print("Tiempo CN con D:", tiempo_cn)


plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_cn):
    t_real = times_to_save_cn[i]*dt_cn
    plt.plot(x, T_snap, '-', linewidth=2, alpha=0.8, label=f"t={t_real:.3f}")
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Difusión 1D con D no variable — Crank-Nicolson")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()



##################### apartado c, vamos a repetir el a solo que ahora D depende de la posición (aunque no se verá afectado por la derivada) #####################


# Función Crank-Nicolson para D(x) dependiente de la posición
def simulate_cn_Dx(D_x, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size                      # número de puntos espaciales
    T = T_initial.copy()                     # copia de la condición inicial
    T = bc_func(T, dx)                       # aplicar condiciones de frontera iniciales
    saved = []                               # lista donde se guardarán los datos

    # Construcción de la matriz como Laplace estándar multiplicada por D_x[j]
    A = lil_matrix((Nx, Nx))                 # inicializar matriz A vacía en formato LIL (para editar filas)
    B = lil_matrix((Nx, Nx))                 # inicializar matriz B vacía en formato LIL

    lower = np.ones(Nx-1)                    # subdiagonal del Laplaciano: coeficientes de T_{j-1}
    main  = -2*np.ones(Nx)                   # diagonal principal: coeficientes de T_j
    upper = np.ones(Nx-1)                    # superdiagonal: coeficientes de T_{j+1}


    L = diags([lower, main, upper], [-1,0,1], shape=(Nx,Nx))/dx**2  # Construimos el lapaciao
    Dmat = diags(D_x, 0, shape=(Nx, Nx), format='csr')              # Dmat es una matriz diagonal que contiene D(x) en cada punto de la malla. Permite multiplicar fácilmente L por D(x) usando producto matricial
    I = diags([np.ones(Nx)], [0], shape=(Nx, Nx), format='csr')     # Matriz identidad Nx x Nx en formato CSR (para construir A y B)

    # Crank–Nicolson: construcción de matrices implícita (A) y explícita (B)
    A = (I - 0.5 * dt * (Dmat @ L)).tolil()   # lado izquierdo: I - 0.5*dt*D*L, convertido a LIL para modificar filas
    B = (I + 0.5 * dt * (Dmat @ L)).tolil()   # lado derecho: I + 0.5*dt*D*L, convertido a LIL

    # Condiciones de contorno en la frontera izquierda (Dirichlet)
    A[0, :] = 0                           # anular toda la fila 0 para imponer condición
    A[0, 0] = 1                           # poner 1 en la diagonal (T_0 = valor impuesto)
    B[0, :] = 0                           # anular fila 0 en B para consistencia
    B[0, 0] = 1                           # poner 1 en B[0,0] también


    A[-1, :] = 0                          # anular toda la última fila antes de imponer la ecuación BC
    A[-1, -1] = 1 + dt                    # entrada A_{N,N} = 1 + dt  (en el informe se deriva la demo)
    A[-1, -2] = -1                        # entrada A_{N,N-1} = -1 para relacionar T_N con T_{N-1}
    B[-1, :] = 0                          # anulamos la última fila en B
    B[-1, -1] = 0                         # definimos B[-1,-1]=0 
     

    # Convertir a CSR para spsolve
    A = csr_matrix(A)
    B = csr_matrix(B)

    # Bucle temporal
    start_time = time()                     # tiempo inicial para medir duración de la simulación
    for n in range(Nt):
        b = B.dot(T)                        # vector lado derecho: b = B * T^n
        T_new = spsolve(A, b)               # resolver sistema A * T^{n+1} = b
        T_new = bc_func(T_new, dx)          # volver a imponer las condiciones de frontera
        T = T_new                           # actualizar T para el siguiente paso
        if times_to_save is not None and n in times_to_save:
            saved.append(T.copy())          # guardar si corresponde

    end_time = time()                        # tiempo final
    tiempo_total = end_time - start_time     # calcular tiempo de ejecución
    return saved, tiempo_total, T           # devolver snapshots, tiempo total y solución final


times_to_save_cn = [0, int(Nt_cn/6), int(Nt_cn/5), int(Nt_cn/4), int(Nt_cn/3), int(Nt_cn/2), Nt_cn-1]
saved_T_cn, tiempo_cn, T_final_cn = simulate_cn_Dx(D_x, dx, dt_cn, Nt_cn, T0, aplicar_bc_cn, times_to_save_cn)
print("Tiempo CN con D(x):", tiempo_cn)


plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_cn):
    t_real = times_to_save_cn[i]*dt_cn
    plt.plot(x, T_snap, '-', linewidth=2, alpha=0.8, label=f"t={t_real:.3f}")
plt.xlabel("x")
plt.ylabel("Temperatura")
plt.title("Difusión 1D con D(x) variable — Crank-Nicolson")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()


##################### apartado d, opcional ##################### 

alpha_values = [0.05, 0.1, 0.2, 0.5]          # Lista de valores de alpha (coeficiente de difusión) que vamos a estudiar
dt = 0.04                                     # Para ver cómo depende de dt vamos a graficar diferentes valores de dt (em diferentes ejecuciones)
for alpha in alpha_values:                    # Iteramos sobre cada valor de alpha
    print(f"Simulación con alpha = {alpha}")  # Mostramos el alpha actual
    
    D_x = alpha * (1 + np.exp(-x**2))         # Definimos el coeficiente de difusión D(x) como función de la posición x
    dt_ftcs = dx**2 / (2*np.max(D_x)) * 0.98  # Calculamos el paso temporal para FTCS usando el criterio de estabilidad
    Nt_ftcs = int(t_final / dt_ftcs)          # Calculamos el número de pasos temporales necesarios para alcanzar t_final
  
    saved_ftcs, tiempo_ftcs, T_ftcs_final = simulate_ftcs_matrix_Dx(D_x, dx, dt_ftcs, Nt_ftcs, T0, aplicar_dirichlet) # Ejecutamos la simulación FTCS
    saved_cn, tiempo_cn, T_cn_final = simulate_cn_Dx(D_x, dx, dt_cn, Nt_cn, T0, aplicar_bc_cn) # Ejecutamos la simulación Crank–Nicolson
    
    # Graficar comparación
    plt.figure(figsize=(10,5))
    plt.plot(x, T_ftcs_final, label='FTCS')
    plt.plot(x, T_cn_final, label='Crank–Nicolson')
    plt.xlabel('x')
    plt.ylabel('Temperatura')
    plt.title(f'Solución final t={t_final}s — alpha={alpha}')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    print(f"FTCS: dt={dt_ftcs:.4e}, pasos={Nt_ftcs}, tiempo_sim={tiempo_ftcs:.4f}s")
    print(f"Crank-Nicolson: dt={dt_cn}, pasos={Nt_cn}, tiempo_sim={tiempo_cn:.4f}s")



##################### extra, estudio de tiempo de simulación requerido por cada método #####################
tiempos_ftcs = []
tiempos_cn = []

for alpha in alpha_values:
    D_x = alpha * (1 + np.exp(-x**2))               # Definimos D(x) para cada alpha como función de la posición x
    
    dt_ftcs = dx**2 / (2*np.max(D_x)) * 0.98        # Calculamos el paso temporal FTCS usando la condición de estabilidad
    Nt_ftcs = int(t_final / dt_ftcs)                # Calculamos el número total de pasos temporales para FTCS
    
    t_ftcs = simulate_ftcs_matrix_Dx(D_x, dx, dt_ftcs, Nt_ftcs, T0, aplicar_dirichlet)[1]    # Ejecutamos la simulación FTCS y extraemos solo el tiempo de ejecución
    t_cn = simulate_cn_Dx(D_x, dx, dt_cn, Nt_cn, T0, aplicar_bc_cn)[1]                       # Ejecutamos la simulación Crank–Nicolson y extraemos solo el tiempo de ejecución
    
    tiempos_ftcs.append(t_ftcs)                     # Guardamos los tiempos de simulación para FTCS
    tiempos_cn.append(t_cn)                         # Guardamos los tiempos de simulación para Crank–Nicolson
    
plt.figure(figsize=(8,5))
plt.plot(alpha_values, tiempos_ftcs, 'o-', label='FTCS')
plt.plot(alpha_values, tiempos_cn, 's-', label='Crank–Nicolson')
plt.xlabel('Alpha (coeficiente de difusión)')
plt.ylabel('Tiempo de simulación [s]')
plt.title('Comparación de tiempo de simulación FTCS vs Crank–Nicolson')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.show()


# Fijamos alpha
alpha = 0.1
D_x = alpha * (1 + np.exp(-x**2))

# Lista de pasos temporales a estudiar
dt_values = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]

# Listas para guardar los tiempos de cada método por separado
tiempos_ftcs = []
tiempos_cn = []

for dt in dt_values:
    # FTCS
    Nt_ftcs = int(t_final / dt)                                                      # número de pasos
    t_ftcs = simulate_ftcs_matrix_Dx(D_x, dx, dt, Nt_ftcs, T0, aplicar_dirichlet)[1] #ejecutamos la función del método
    tiempos_ftcs.append(t_ftcs)                                                      # guardamos tiempos de simulación 
    
    # Crank–Nicolson
    Nt_cn = int(t_final / dt)                                       # número de pasos
    t_cn = simulate_cn_Dx(D_x, dx, dt, Nt_cn, T0, aplicar_bc_cn)[1] #ejecutamos la función del método
    tiempos_cn.append(t_cn)                                         # guardamos tiempos de simulación 

# Graficar comparación de tiempos
plt.figure(figsize=(8,5))
plt.plot(dt_values, tiempos_ftcs, 'o-', label='FTCS')
plt.plot(dt_values, tiempos_cn, 's-', label='Crank–Nicolson')
plt.xlabel('Paso temporal Δt')
plt.ylabel('Tiempo de simulación [s]')
plt.title('Comparación de tiempo de simulación FTCS vs Crank–Nicolson')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.show()

##################### extra, estudio de cómo aumenta el error al aumentar dt #####################
#para ello lo que vamos a ahcer es tomar un dt muy pequeño como solución y calcularemos
#la diferencia entre la señal y un dt mayor para ver cómo varía al aumentar el dt en cada método

def error_vs_dt_ftcs(D_x, dx, T0, aplicar_bc, dt_min, dt_max, num_dt, t_final): #Calcula el error numérico del método FTCS con D(x) al variar el paso temporal dt.
    # generamos una solución que tomaremos como referencia para comparar con el resto
    Nt_ref = int(t_final / dt_min)                 # número de pasos con dt muy pequeño
    print(f"Generando solución de referencia con Nt = {Nt_ref} ...")
    T_ref = simulate_ftcs_matrix_Dx(D_x, dx, dt_min, Nt_ref, T0, aplicar_bc)[2] # ejecutar FTCS con dt extremadamente pequeño 


    dts = np.linspace(dt_min, dt_max, num_dt)      # creamos una lista dt desde muy pequeño hasta valor grande

    errores_max = []                               # aquí guardamos error máximo para cada dt
    errores_med = []                               # aquí guardamos error medio (L1) para cada dt


    #simulamos y obtenemos el valor de cada dt
    for dt in dts:                                  # recorremos cada valor de dt
        Nt = int(t_final / dt)                      # número de pasos para ese dt
        T_num = simulate_ftcs_matrix_Dx(D_x, dx, dt, Nt, T0, aplicar_bc)[2] # simular FTCS con ese dt
        err = np.abs(T_num - T_ref)                 # calcular error absoluto punto a punto
        errores_max.append(np.max(err))             # guardamos error máximo del dominio
        errores_med.append(np.mean(err))            # guardamos error medio (promedio del error absoluto)


    plt.figure(figsize=(9,6))
    plt.loglog(dts, errores_max, 'o-', label="Error máximo")  # gráfica doble logaritmo
    plt.loglog(dts, errores_med, 's-', label="Error medio")
    plt.xlabel("dt")
    plt.ylabel("Error")
    plt.title("Error numérico FTCS con D(x) al variar dt")
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return dts, errores_max, errores_med


dt_min = 1e-5             # muy pequeño, lo tomaremos como solución casi exacta
dt_max = 5e-3             # grande FTCS empezará a volverse inestable
num_dt = 15               # número de puntos del estudio

dts, emax, emed = error_vs_dt_ftcs(D_x=D_x, dx=dx, T0=T0,aplicar_bc=aplicar_dirichlet,dt_min=dt_min,dt_max=dt_max,num_dt=num_dt,t_final=t_final)

def error_vs_dt_cn(D_x, dx, T0, aplicar_bc, dt_min, dt_max, num_dt, t_final):

    Nt_ref = max(10, int(t_final / dt_min))  # calcular número de pasos para la solución de referencia, evitar pasos excesivos que generen NaN
    print(f"Generando solución de referencia CN con Nt = {Nt_ref} ...")  # informar cuántos pasos se usarán para la referencia
    T_ref = simulate_cn_Dx(D_x, dx, dt_min, Nt_ref, T0, aplicar_bc)[2]  # ejecutar CN con dt muy pequeño y obtener la solución final como referencia

    dts = np.linspace(dt_min + dx**2 / (2*np.max(D_x))*0.001, dt_max, num_dt)  # crear lista de dt a estudiar, evitando dt extremadamente pequeño
    errores_max = []  # lista donde guardaremos el error máximo para cada dt
    errores_med = []  # lista donde guardaremos el error medio (promedio)

    for dt in dts:
        Nt = max(1, int(t_final / dt))  # calcular número de pasos para este dt, evitar Nt=0
        T_num = simulate_cn_Dx(D_x, dx, dt, Nt, T0, aplicar_bc)[2]  # ejecutar CN con este dt y obtener la solución final

        err = np.abs(T_num - T_ref)  # calcular error absoluto punto a punto
        errores_max.append(np.max(err))  # error máximo en todo el dominio
        errores_med.append(np.mean(err))  # error medio (promedio)

        print(f"dt={dt:.5e}, Error máximo={errores_max[-1]:.5e}, Error medio={errores_med[-1]:.5e}")  # mostrar errores por consola

    plt.figure(figsize=(9,6))  # crear figura para graficar errores
    plt.loglog(dts, errores_max, label="Error máximo")  # graficar error máximo en escala log-log
    plt.loglog(dts, errores_med, label="Error medio")  # graficar error medio en escala log-log
    plt.xlabel("dt")  # etiqueta eje x
    plt.ylabel("Error")  # etiqueta eje y
    plt.title("Error numérico Crank–Nicolson con D(x) al variar dt")  # título de la gráfica
    plt.grid(True, alpha=0.6)  # activar rejilla con transparencia
    plt.legend()  # mostrar leyenda
    plt.show()  # mostrar la gráfica

    return dts, errores_max, errores_med  # devolver dt estudiados y errores calculados

dt_min = dx**2 / (2*np.max(D_x))*0.1  # dt muy pequeño, referencia precisa
dt_max = dx**2 / (2*np.max(D_x)) * 0.99  # dt grande permitido por Crank–Nicolson
num_dt = 8  # número de dt a estudiar en el gráfico

dts_cn, emax_cn, emed_cn = error_vs_dt_cn(D_x=D_x, dx=dx, T0=T0, aplicar_bc=aplicar_bc_cn, dt_min=dt_min, dt_max=dt_max, num_dt=num_dt, t_final=t_final)  # ejecutar estudio de error



##################### extra, derivar correctamente D(x) #####################

def simulate_ftcs_matrix_Dx(D_x, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size                     # número de puntos espaciales
    T = T_initial.copy()                    # copia del estado inicial de temperatura
    T = bc_func(T)                          # aplicar condiciones de frontera iniciales
    saved = []                              # lista donde guardaremos snapshots
    

    ##### Construcción del operador Laplaciano para D(x) #####
    # Para D(x) variable se usa la aproximación discreta:
    #\frac{\partial}{\partial x} \left( D \frac{\partial T}{\partial x} \right) = (D_{j+1/2}*(T_{j+1}-T_j) - D_{j-1/2}*(T_j-T_{j-1}))/dx^2
    
    D_plus_half  = 0.5 * (D_x[1:] + D_x[:-1])   # D_{j+1/2} para j=0..Nx-2
    lower = -D_plus_half                 # subdiagonal (-D_{j-1/2})
    upper = -D_plus_half                    # superdiagonal (-D_{j+1/2})
    main  = np.zeros(Nx)                        # diagonal principal
    main[1:-1] = -(lower[:-1] + upper[1:])               # T_j * (-D_{j-1/2} - D_{j+1/2}) según discretización
    
    D_plus_half = 0.5 * (D_x[1:] + D_x[:-1])  # Nx-1 elementos


    # Construcción de la matriz sparse L
    L_in = diags([lower, main, upper], offsets=[-1, 0, 1], shape=(Nx, Nx), format='csr')
    
    # Fijar condiciones de Dirichlet en los bordes
    L = L_in.tolil()         # convertir a LIL para poder editar filas
    L[0, :] = 0.0            # anular fila izquierda
    L[0, 0] = 1.0            # poner 1 en diagonal izquierda
    L[-1, :] = 0.0           # anular fila derecha
    L[-1, -1] = 1.0          # poner 1 en diagonal derecha
    L = csr_matrix(L)        # volver a CSR para multiplicaciones rápidas
    
    start_time = time()
    for n in range(Nt):
        # Actualización FTCS: T^{n+1} = T^n + dt/dx^2 * L * T^n
        T_new = T + dt/dx**2 * L.dot(T)
        T_new = bc_func(T_new)              # reimponer condiciones de frontera
        T = T_new                           # actualizar temperatura para el siguiente paso
        
        if times_to_save is not None and n in times_to_save:
            saved.append(T.copy())          # guardar
    
    end_time = time()
    tiempo_total = end_time - start_time
    
    return saved, tiempo_total, T

times_to_save_ftcs = [0, int(Nt_ftcs/6), int(Nt_ftcs/5), int(Nt_ftcs/4), int(Nt_ftcs/3), int(Nt_ftcs/2), Nt_ftcs-1]
saved_T_ftcs, tiempo_ftcs, T_final_ftcs = simulate_ftcs_matrix_Dx(D_x=D_x,dx=dx,dt=dt_ftcs,Nt=Nt_ftcs,T_initial=T0,bc_func=aplicar_dirichlet,times_to_save=times_to_save_ftcs)
print("Tiempo FTCS (matriz sparse, D(x)):", tiempo_ftcs)
plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_ftcs):
    t_real = times_to_save_ftcs[i] * dt_ftcs
    plt.plot(x, T_snap, linewidth=2, alpha=0.9, label=f"t = {t_real:.3f}s")
plt.title("Difusión 1D — Método FTCS con D(x) variable", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("Temperatura", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()



def simulate_cn_matrix_Dx(D_x, dx, dt, Nt, T_initial, bc_func, times_to_save=None):
    Nx = T_initial.size           # número de puntos espaciales
    T = T_initial.copy()          # copia del estado inicial
    T = bc_func(T, dx)            # aplicar condiciones de frontera iniciales
    saved = []                    # lista para guardar snapshots

    #### Construcción del operador L de la derivada espacial:
    D_plus_half = 0.5 * (D_x[1:] + D_x[:-1]) # D_{j+1/2} = 0.5*(D_j + D_{j+1})

    # Subdiagonal y superdiagonal según la discretización central
    lower = -D_plus_half     # subdiagonal (-D_{j-1/2})
    upper = -D_plus_half     # superdiagonal (-D_{j+1/2})
    
    # Diagonal principal
    main = np.zeros(Nx)
    # Para puntos interiores: T_j * (-D_{j-1/2} - D_{j+1/2})
    main[1:-1] = -(lower[:-1] + upper[1:])

    # Construcción de la matriz sparse L tridiagonal
    L = diags([lower, main, upper], offsets=[-1, 0, 1], shape=(Nx, Nx), format='csr')

    # Construcción de las matrices Crank–Nicolson
    # A * T^{n+1} = B * T^n
    I = diags([np.ones(Nx)], [0], format='csr')    # matriz identidad

    # T^{n+1} - T^n = 0.5*dt/dx^2 * L*(T^{n+1} + T^n) (esto proviene de la discretización temporal CN)
    A = (I - 0.5*dt/dx**2 * L).tolil() # Matriz implícita A: I - 0.5*dt/dx^2*L
    B = (I + 0.5*dt/dx**2 * L).tolil() # Matriz explícita B: I + 0.5*dt/dx^2*L
    
    A[0, :] = 0                           # anular toda la fila 0 para imponer condición
    A[0, 0] = 1                           # poner 1 en la diagonal (T_0 = valor impuesto)
    B[0, :] = 0                           # anular fila 0 en B para consistencia
    B[0, 0] = 1                           # poner 1 en B[0,0] también


    A[-1, :] = 0                          # anular toda la última fila antes de imponer la ecuación BC
    A[-1, -1] = 1 + dt                    # entrada A_{N,N} = 1 + dt  (en el informe se deriva la demo)
    A[-1, -2] = -1                        # entrada A_{N,N-1} = -1 para relacionar T_N con T_{N-1}
    B[-1, :] = 0                          # anulamos la última fila en B
    B[-1, -1] = 0                         # definimos B[-1,-1]=0 
     
    # Convertir A y B a formato CSR para resolver el sistema rápidamente
    A = csr_matrix(A)
    B = csr_matrix(B)


    start_time = time()
    for n in range(Nt):
        b = B.dot(T)                  # lado derecho B*T^n
        T_new = spsolve(A, b)         # resolver A*T^{n+1} = b
        T_new = bc_func(T_new, dx)    # reimponer condiciones de frontera
        T = T_new                     # actualizar temperatura

        if times_to_save is not None and n in times_to_save:
            saved.append(T.copy())    # guardar snapshot si corresponde

    end_time = time()
    tiempo_total = end_time - start_time

    return saved, tiempo_total, T



#cuidado al cambiar dt_cn y alpha ya que la gráfica varía bastante y oscila, con dt_cn = 0.02 y alpha = 0.001 creo que la solucion es interesante
dt_cn = 0.02                          # paso temporal para Crank–Nicolson
alpha = 0.001
Nt_cn = int(t_final / dt_cn)            # número de pasos temporales totales
D_x = alpha * (1 + np.exp(-x**2))
alpha = 0.1
times_to_save_cn = [0, int(Nt_cn/6), int(Nt_cn/5), int(Nt_cn/4), int(Nt_cn/3), int(Nt_cn/2), Nt_cn-1]
saved_T_cn, tiempo_cn, T_final_cn = simulate_cn_matrix_Dx(D_x=D_x, dx=dx, dt=dt_cn, Nt=Nt_cn, T_initial=T0, bc_func=aplicar_bc_cn,times_to_save=times_to_save_cn)
print("Tiempo Crank–Nicolson (matriz sparse):", tiempo_cn)
plt.figure(figsize=(10,6))
for i, T_snap in enumerate(saved_T_cn):
    t_real = times_to_save_cn[i] * dt_cn
    plt.plot(x, T_snap, linewidth=2, alpha=0.9, label=f"t = {t_real:.3f}s")
plt.title("Difusión 1D — Método Crank–Nicolson con D(x) variable", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("Temperatura", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()




