# La documentacion de esta libreria se puede visitar en:


# LIBRERIAS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


''' Funcion Beta_1: 
    Nos perimte encontrar el valor de Beta 1 que es la proporcion de c que se toma para calcular a, en el diagrama de tensiones
    rectangulares de Whitney - ACI 318-19
    -Parametros de entrada:
        fc: resistencia a compresión del hormigon medida en psi
    -Valores de salida:
        Beta_1: es un valor flotante adimensional
'''
def Beta_1(fc):
  if fc<2500:
    print("ERROR, fc debe ser mayor a 2500psi")
  elif fc<4000:
    return 0.85
  elif fc<8000:
    return 0.85 - (fc - 4000) / 1000 * 0.05
  else:
    return 0.65

''' Funcion Phi_col: 
    Nos perimte encontrar el valor del factor de reduccion de resistencia Phi a emplear para la columna.
    En el ACI 318-19 este valor no es constante en columnas, porque depende de la ductilidad y de la falla, y de las
    incertidubres. En columnas armadas con espiral, la ductilidad es mayor a aquellas armadas con estribos, por 
    lo tanto el codigo permite valores mayores de Phi. 
    En general cuando la falla está dominada por traccion el Phi es 0.9 para todos los casos, pero cuando la falla esta
    dominada por compresión, segun si es zunchada o estribada el phi sera 0.75 o 0.65 respectivamente. Hay ademas una zona 
    de transicion entre falla dominada por traccion y falla dominada por compresion en donde el phi es una interpolacion lineal.
    
    Una salvedad que no pertenece al ACI 318 pero que esta aceptado, es que cuando se dimensionan pilotes en vez de columnas
    estos valores pueden ser distitnos segun la condicion de instalacion del pilote. 
    
'''
def Phi_col(zeta_t, zeta_y, zeta_c, tipo_armado_transversal):
  #zeta_t: deformacion unitaria de la barra más traccionada. Positivo = traccion
  #zeta_y: deformacion unitaria de fluencia del acero
  #zeta_c: deformacion unitaria maxima a compresión del hormigon
  #tipo_armado_transversal: indica si esta armada con estribos (=0), o con zunchos (=1)
  if (tipo_armado_transversal == False):
    if zeta_t < zeta_y:
      Phi_col = 0.65
    elif zeta_t < zeta_y + zeta_c:
      Phi_col = 0.65 + (0.9 - 0.65) * (zeta_t - zeta_y) / zeta_c
    else:
      Phi_col = 0.9
  else:
    if zeta_t < zeta_y:
      Phi_col = 0.75
    elif zeta_t < zeta_y + zeta_c:
      Phi_col = 0.75 + (0.9 - 0.75) * (zeta_t - zeta_y) / zeta_c
    else:
      Phi_col = 0.9
  return Phi_col


def Diagrama_Interaccion(fc, zeta_c, fy, zeta_y,Diametro, n_bp, A_bar_p, n_bc, A_bar_c, tipo_armado_transversal, rec_mec):
  # PARA QUE SIRVE ESTA FUNCION?
  #Esta funcion arma todos los valores del diagrama de interaccion a partir de los parametros del pilote
  #PARAMETROS DE ENTRADA:
  #zeta_c: deformación unitaria máxima del concreto, que en el ACI 318 es 0.003
  #zeta_y: deformacion unitaria de fluencia del acero

  #OUTPUT:
  #Array de 4 vectores columna: Mn, Pn, Phi*Mn, Phi*Pn

  ID=np.arange(0,n_bp+1) # Vector de identificador de las barras perimetrales, cada barra perimetral tiene su ID. El ID 0 le pertenece a la barra central (que puede no haber)
  alpha=np.zeros(n_bp) # Angulo entre la vertical y la barra i
  X=np.zeros(n_bp+1) # Coordenadas X de las barras
  Y=np.zeros(n_bp+1) # Coordenadas Y de las barras
  d_barra=np.zeros(n_bp+1) # Declarar y dimensionar el arreglo para las alturas d de cada barra
  zetas=np.zeros(n_bp+1) # Declarar y dimensionar el arreglo para las deformaciones zeta de cada barra
  Tensiones=np.zeros(n_bp+1) # Declarar y dimensionar el arreglo para las tensiones de cada barra
  Fuerzas=np.zeros(n_bp+1) # Declarar y dimensionar el arreglo para las fuerzas de cada barra
  Momento_barras=np.zeros(n_bp+1) # Declarar el vector para almacenar los momentos que genera cada barra
  max_d= 0 # Declarar la variable para el valor máximo de d
  max_d_ID= 0 # Declarar la variable para el ID de la barra con el valor máximo de d
  c = 0 # Declarar la variable para la distancia al eje neutro
  a = 0 # Declarar la variable para la variable "a"
  theta = 0 # Declarar la variable para el ángulo "theta"
  Phi = 0 # Factor de reduccion de capacidad para columna
  Beta_1_adopt = Beta_1(fc)
  As_total = n_bp * A_bar_p + n_bc * A_bar_c
  tamanio_salida = 5000

  Pn=np.zeros(tamanio_salida)
  Mn=np.zeros(tamanio_salida)
  phi_Pn=np.zeros(tamanio_salida)
  phi_Mn=np.zeros(tamanio_salida)
  zeta_barra_mas_traccionada=np.zeros(tamanio_salida)

  for contador in range(0, tamanio_salida):
    zeta_barra_mas_traccionada[contador] = -zeta_c + (zeta_c+50*zeta_y) * contador / tamanio_salida

  #Para encontrar la máxima traccion tengo que proponer una deformacion muy grande de la barra mas traccionada, es más practico
  #calcularla con la formula analitica y unir ese valor al final del vector de deformaciones

  #Ahora para cada valor del grafico debemos iterar (cada valor del grafico corresponde a una deformacion de la
  #barra mas traccionada

  #Calcular los ángulos alpha, así como las coordenadas X e Y y alturas d_barra para cada barra perimetral
  for i in range(0, n_bp):
    alpha[i] = i * 2 * np.pi / n_bp # Utilizamos la constante Pi como se definió anteriormente
    X[i] = np.sin(alpha[i]) * (Diametro / 2 - rec_mec)
    Y[i] = np.cos(alpha[i]) * (Diametro / 2 - rec_mec)
    # Calcular las alturas d_barra de cada barra
    d_barra[i] = Diametro / 2 - Y[i]

  X[n_bp] = 0 # La barra central
  Y[n_bp] = 0 # La barra central
  d_barra[n_bp] = Diametro / 2 # La barra central

  # Calcular el valor máximo de d entre todas las barras, que va a ser el d, y el ID de la barra correspondiente
  max_d = d_barra[0] # Inicializar max_d con el primer valor de d_barra
  for i in range(0, n_bp+1):
    # Calcular el valor máximo de d y el ID de la barra correspondiente
    if d_barra[i] > max_d:
      max_d = d_barra[i]

  #El primer valor del poligono que corresponde a presion pura lo puedo calcular de antemano
  Pn[0] = (0.85 * fc * (Diametro * Diametro * np.pi / 4 - As_total) + fy * 1000 * As_total) / 1000
  Mn[0] = 0
  Phi = Phi_col(0, zeta_y, zeta_c, tipo_armado_transversal)
  if tipo_armado_transversal == False:
    phi_Pn_max = 0.8 * Phi * Pn[0]  # Este es el truncamiento de phi Pn que impone el ACI 318
  else:
    phi_Pn_max = 0.85 * Phi * Pn[0] # Este es el truncamiento de phi Pn que impone el ACI 318
  phi_Pn[0] = phi_Pn_max
  phi_Mn[0] = 0


  for j in range(1, tamanio_salida - 2): # Los ultimos 2 valores los calculo a parte, el anteultimo le corresponde a traccion pura, y el ultimo es el primero para cerrar la poligonal
    Pn[j] = 0
    Mn[j] = 0
    # Calcular las deformaciones zetas, tensiones y fuerzas en las barras
    for i in range(0, n_bp): # Esto es para las barras perimetrales que tienen todas el mismo area
      # Calcular las deformaciones zetas
      zetas[i] = zeta_c - ((zeta_c + zeta_barra_mas_traccionada[j]) / max_d) * d_barra[i]
      # Calcular las tensiones en las barras en psi
      if zetas[i] > zeta_y:
        Tensiones[i] = fy * 1000
      elif zetas[i] < -zeta_y:
        Tensiones[i] = -fy * 1000
      else:
        Tensiones[i] = zetas[i] / zeta_y * fy * 1000
      # Calcular las fuerzas en las barras
      Fuerzas[i] = Tensiones[i] * A_bar_p / 1000 # Fuerza en kip
      # En esta fuerza falta restar la tension del hormigon en caso que la barra este en la zona de a, pero eso se hace mas abajo, cuando se determina a.
    #Ahora calculo para la barra central que puede ser distinta
    zetas[n_bp] = zeta_c - ((zeta_c + zeta_barra_mas_traccionada[j]) / max_d) * d_barra[n_bp]
    if zetas[n_bp] > zeta_y:
      Tensiones[n_bp] = fy * 1000
    elif zetas[n_bp] < -zeta_y:
      Tensiones[n_bp] = -fy * 1000
    else:
      Tensiones[n_bp] = zetas[n_bp] / zeta_y * fy * 1000
    Fuerzas[n_bp] = Tensiones[n_bp] * A_bar_c * n_bc / 1000 # Fuerza en kip
    # En esta fuerza falta restar la tension del hormigon en caso que la barra este en la zona de a, pero eso se hace mas abajo, cuando se determina a.
    #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    # Calcular la distancia al eje neutro "c" usando max_d en lugar de d
    #Este if es para evitar la indeterminacion:
    if zeta_barra_mas_traccionada[j] == -zeta_c:
      c = Diametro / Beta_1_adopt
    else:
      c = max_d / (zeta_barra_mas_traccionada[j] + zeta_c) * zeta_c

    #Calcular la variable "a" como a = c * Beta_1
    if c * Beta_1_adopt > Diametro:
      a = Diametro
    else:
      a = c * Beta_1_adopt

    # Ahora que conocemos el a, vemos que barras quedaron en la zona de compresion y a esas les restamos fc*Asbar, ya que sino estariamos sumando mas fuerza de compresion
    for i in range(0, n_bp): # Esto es para las barras perimetrales que tienen todas el mismo area
      # Calcular las fuerzas en las barras
      if d_barra[i] <= a: # ESTA BIEN ESTA CONDICION??????? MEDIO RARA
        Fuerzas[i] = Fuerzas[i] - 0.85 * fc * A_bar_p / 1000 # Fuerza en kip
        #Fuerzas[n_bp] =0 # PRUEBA A VER SI ASI DA COMO ETABS, NO SE ME OCURRE OTRA COSA
      Momento_barras[i] = Fuerzas[i] * Y[i] #kip.in
      Pn[j] = Pn[j] + Fuerzas[i]
      Mn[j] = Mn[j] + Momento_barras[i]
    #Ahora calculo para la barra central que puede ser distinta
    if d_barra[n_bp] <= a: # ESTA BIEN ESTA CONDICION??????? MEDIO RARA
      Fuerzas[n_bp] = Fuerzas[n_bp] - 0.85 * fc * A_bar_c * n_bc / 1000 # Fuerza en kip
      #Fuerzas[n_bp] =0 # PRUEBA A VER SI ASI DA COMO ETABS, NO SE ME OCURRE OTRA COSA
    Momento_barras[n_bp] = Fuerzas[n_bp] * Y[n_bp] #kip.in
    Pn[j] = Pn[j] + Fuerzas[n_bp]
    Mn[j] = Mn[j] + Momento_barras[n_bp]
    #############################################################################################

    # Calcular el ángulo "theta" en radianes

    theta = np.arccos((Diametro / 2 - a) / (Diametro / 2))

    # Calcular el área comprimida A_comprimida
    A_comprimida = Diametro * Diametro * (theta - np.sin(theta) * np.cos(theta)) / 4
    # Calcular la fuerza resultante de compresión (Cc)
    Cc = 0.85 * fc * A_comprimida / 1000 # fuerza de compresion en kip
    Pn[j] = Pn[j] + Cc
    # Calcular el baricentro del área comprimida (Y_raya)
    Y_raya = Diametro**3 * np.sin(theta)** 3 / (12 * A_comprimida)
    Mn[j] = (Mn[j] + Y_raya * Cc) / 12 # kip.ft

    #Calculo de los valores admisibles

    Phi = Phi_col(zeta_barra_mas_traccionada[j], zeta_y, zeta_c, tipo_armado_transversal)
    if Phi * Pn[j] < phi_Pn_max:
      phi_Pn[j] = Phi * Pn[j]
    else:
      phi_Pn[j] = phi_Pn_max
    phi_Mn[j] = Phi * Mn[j]

    #El anteultimo valor de la poligonal es la máxima capacidad a traccion que la calculamos analiticamente:
    Pn[tamanio_salida - 2] = -fy * As_total
    Mn[tamanio_salida - 2] = 0
    Phi = Phi_col(zeta_barra_mas_traccionada[tamanio_salida - 2], zeta_y, zeta_c, tipo_armado_transversal)
    phi_Pn[tamanio_salida - 2] = Phi * Pn[tamanio_salida - 2]
    phi_Mn[tamanio_salida - 2] = Phi * Mn[tamanio_salida - 2]

    #Para que la poligonal sea cerrada, el ultimo valor debe coincidir con el primero
    Pn[tamanio_salida - 1] = Pn[0]
    Mn[tamanio_salida - 1] = Mn[0]
    phi_Pn[tamanio_salida - 1] = phi_Pn[0]
    phi_Mn[tamanio_salida - 1] = phi_Mn[0]
  return pd.DataFrame({'ØMn':phi_Mn, 'ØPn':phi_Pn, 'Mn':Mn, 'Pn':Pn})

