from init import initValue
import numpy as np

phiInit = []
def systemHeatEXE():
    #---------------------------------------------------------
    # Isso constroi o vetor com a condição de inicial
    #---------------------------------------------------------
    for i in range(initValue.jInterator):
        if i < int(initValue.jInterator/2):
            phiInit.append(initValue.T1)
        if i == int(initValue.jInterator/2):
            phiInit.append(0.5) 
            # consideramos a Temperatira 0,5, na discontinuidade
        if i > int(initValue.jInterator/2):
            phiInit.append(initValue.T2)
    #---------------------------------------------------------
    phiN = [phiInit]
    t = [0]
    #---------------------------------------------------------
    # Metodo de Thomas
    #--------------------------------------------------------- 
    for k in range(initValue.nInterator):
        t.append(t[k]+ initValue.deltT)
        phiJ = []
        Pitere = 0; Qitere = initValue.T1
        Pjth = [] ; Qjth = []
        #-----------------------------------------------------
        # um for para calcular todos os valores de P e Q
        #-----------------------------------------------------
        for i in range(initValue.jInterator):
            if i == 0:
                #-------------------------------------------
                # Definição de P1 e Q1
                #-------------------------------------------
                P1 = 0; Q1 = initValue.T1
                Pjth.append(P1) ; Qjth.append(Q1)
                #-------------------------------------------
            elif i > 0 and i < (initValue.jInterator-1):
                #-------------------------------------------
                # Definição de Pj e Qj de j=2 ah j = J-1 
                #-------------------------------------------
                dj = ((1-initValue.betha)*(0.5*initValue.C*(1+initValue.sigma)+
                initValue.s)*phiN[k][i-1]+(1-2*(1-initValue.betha)*(0.5*
                initValue.C*initValue.sigma +initValue.s ))*
                phiN[k][i]+(1-initValue.betha)*(0.5*initValue.C*
                (initValue.sigma-1)+initValue.s)*phiN[k][i+1])

                Qitere= ((dj + initValue.cj*Qitere)/
                (initValue.aj -initValue.cj*Pitere))

                Pitere = initValue.bj/(initValue.aj -initValue.cj*Pitere)

                Pjth.append(Pitere) ; Qjth.append(Qitere)
                #-------------------------------------------
            else:
                #-------------------------------------------
                # Definição de PJ e QJ
                #-------------------------------------------
                PJ = 0; QJ = initValue.T2
                Pjth.append(PJ) ; Qjth.append(QJ)
        #-----------------------------------------------------
        

        #-----------------------------------------------------
        # um for para calcular phi no n+1
        #-----------------------------------------------------
        for i in range(initValue.jInterator-1,-1,-1):
            if i == initValue.jInterator-1:
                #-------------------------------------------
                # phi de n+1 em J
                #-------------------------------------------
                phiJ.append(Pjth[i]*phiN[k][i] + Qjth[i])
                #-------------------------------------------
            elif i > 0 and i < (initValue.jInterator-1):
                #-------------------------------------------
                # phi n+1 de j = J-1 ah j = 2 
                #-------------------------------------------
                phiJ.insert(0,Pjth[i]*phiJ[0] + Qjth[i])
                #-------------------------------------------
            else:
                #-------------------------------------------
                # phi de n+1 em 1
                #-------------------------------------------
                phiJ.insert(0,Pjth[i]*phiN[k][i] + Qjth[i])
                phiN.append(phiJ)
        #-----------------------------------------------------

    #---------------------------------------------------------
    xi = -2
    xf = 2
    T = []
    x = np.arange(xi, xf+initValue.deltX, initValue.deltX)
    #---------------------------------------------------------
    # expressao analitica
    #---------------------------------------------------------
    for i in range(initValue.jInterator):
        termSom = []
        for k in range(1,initValue.maxEx+1):
            termSom.append((1/(2*k-1))*(np.e**((-initValue.alfa*((2*k-1)**2)*
            (np.pi**2)*initValue.tmax)/(initValue.L**2)))*
            np.sin((2*k-1)*((np.pi*(x[i]-initValue.u*
            initValue.tmax))/(initValue.L))))
        somatorio = sum(termSom)
        T.append(0.5 - (2/np.pi)*somatorio)
    #---------------------------------------------------------

    quad = 0.0
    #---------------------------------------------------------
    # Erro Medio Quadratico
    #---------------------------------------------------------
    for i in range(initValue.jInterator):
        quad = (quad + ((phiN[initValue.nInterator][i]-T[i])**2))
    EMQ = np.sqrt(quad/initValue.jInterator)
    #---------------------------------------------------------


    print('JMAX=\t',initValue.jInterator,' MAXEX=\t',initValue.maxEx)   
    print('TMAX=\t',initValue.tmax,' T1=\t',initValue.T1,'T2=\t',initValue.T2)            
    print('BETA=\t',initValue.betha,' SIGMA=\t',initValue.sigma)            
    print('DELT=\t',initValue.deltT,' DELTX=\t',initValue.deltX)            
    print('S=\t','%.3f' % initValue.s,' ALPH=\t',initValue.alfa)            
    print('U=\t','%.3f' % initValue.u,' C= ',initValue.C)            
    print()
    for i in range(initValue.nInterator+1):
        var1 = "\t".join(map(str,np.around(phiN[i],3)))
        print('t= ','%.3f' % t[i],' TN=',var1)
    print()
    var2 = "\t".join(map(str,np.around(T,3)))
    print('t= ','%.3f' % initValue.tmax,' TE=',var2)
    print()
    print('EMQ= ', '{:0.3E}'.format(EMQ))
    #---------------------------------------------------------
    # FIM DO PROGRAMA
    #---------------------------------------------------------