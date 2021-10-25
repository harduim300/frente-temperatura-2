#---------------------------------------------------------
# Bibliotecas necessárias
#---------------------------------------------------------
import numpy as np
#---------------------------------------------------------
# pip install numpy
#---------------------------------------------------------
# Caso pip não funcione usar pip3
#---------------------------------------------------------


#---------------------------------------------------------
# constantes
#---------------------------------------------------------
alfa  = 0.1   # difusividade termica  
L     = 4.0   # comprimento maximo 
T1    = 1.0   # temperatura inicial     
T2    = 0.0   # temperatura final     
#---------------------------------------------------------


#---------------------------------------------------------
# Inicialização das variaveis
#---------------------------------------------------------
betha = 1.0/2.0   # valor de betha   
deltT = 0.04  # delta T
deltX = 0.4   # delta X
tmax  = 1.0   # tempo maxima 
maxEx = 100   # k da serie infinita
sigma = 0.0  
#---------------------------------------------------------
# valor de s
#---------------------------------------------------------
# se o valor é definido
# s = 1.0
# deltX = np.sqrt(alfa*(deltT)/s)
#---------------------------------------------------------
# se o valor não é definido
s     = alfa*(deltT)/(deltX**2)     
#---------------------------------------------------------

#---------------------------------------------------------
# Se o valor de C é defindo
#---------------------------------------------------------
# C     = 0.0
# u     = C*deltX/deltT # velocidade de advecção
#---------------------------------------------------------
# Se o valor de u é definido e C não
#---------------------------------------------------------
u = 0.25  # velocidade de advecção
C   = u*(deltT/deltX);  
#---------------------------------------------------------





#---------------------------------------------------------
# variaveis de iteração
#---------------------------------------------------------
jInterator = int(1 + 4/deltX)    # -2 <= x <= 2
nInterator = int(tmax/deltT)     #  0 <= t <= tmax 
#---------------------------------------------------------

#---------------------------------------------------------
# Como apenas dj eh variavel podemos defini-los como constantes
# no inicio do programa
#---------------------------------------------------------
aj =  (1+2*betha*(0.5*C*sigma +s))                
bj =  betha*(0.5*C*(sigma - 1) + s)                 
cj =  betha*(0.5*C*(1 + sigma) + s)    
#---------------------------------------------------------

#---------------------------------------------------------
# CLASSIFICACAO
#---------------------------------------------------------
classitip = ''
classadv = ''
if betha == 0 and sigma == 0:
    classitip = 'FTCS'
    classadv = 'DIFERENCAS CENTRADAS'
if betha == 0 and sigma == 1:
    classitip = 'EXPLICITO'
    classadv = 'UPWIND'
if betha == 0.5:
    classitip = 'CRANK-NICOLSON'
    if sigma == 0:
        classadv = 'DIFERENCAS CENTRADAS'
    else:
        classadv = 'UPWIND'
if betha == 1:
    classitip = 'IMPLICITO'
    if sigma == 0:
        classadv = 'DIFERENCAS CENTRADAS'
    else:
        classadv = 'UPWIND'
#---------------------------------------------------------

print('EQUACAO DE ADVECCAO-DIFUSAO:  ',classitip)
print('TERMO ADVECTIVO:  ',classadv)
print()