from structurelab.concrete.concrete_piles import *

# Datos de entrada:
# USADOS PARA FLEXION (Y PARA CORTE ALGUNOS TAMBIEN)
fc=5000 #psi
zeta_c = 0.003
fy=60 # ksi
Ey=29000 # ksi
zeta_y= fy/Ey
Diametro = 20
n_bp=8
A_bar_p=0.79
n_bc=0
A_bar_c=0
tipo_armado_transversal=0
rec_mec=1.5+0.5+0.5

# DIAGRAMA DE ETABS
# SIN PHI
Pn=[1343.163,1343.163,1308.997,1091.644,827.368,522.97,309.48,94.765,-88.646,-294.928,-379.2]
Mn=[0,105.3584,195.1761,269.829,319.133, 339.9481,317.2957,257.5067,173.7797,59.0473,0]


# POR PHI
phi_Pn = [873.056,873.056,850.848,709.569,537.789,339.931,233.844,85.289,-79.782,-265.435,-341.28]
phi_Mn = [0,68.4829,126.8645,175.3888,207.4364,220.9663,239.7497,231.756,156.4018,53.1425,0]

df =Diagrama_Interaccion(fc, zeta_c, fy, zeta_y,Diametro, n_bp, A_bar_p, n_bc, A_bar_c, tipo_armado_transversal, rec_mec)
plt.plot(df.loc[:,"ØMn"], df.loc[:,"ØPn"], 'r') # plotting t, a separately
plt.plot(df.loc[:,"Mn"], df.loc[:,"Pn"], 'b') # plotting t, b separately
plt.plot(phi_Mn, phi_Pn, 'g') # ESTE SALE DE ETABS
plt.plot(Mn, Pn, 'g') # ESTE SALE DE ETABS
plt.show()

