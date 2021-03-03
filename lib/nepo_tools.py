import numpy as np
import cmath

def dPhasor(r, phi, polar=True):
    z = cmath.rect(r,phi*cmath.pi/180)
    return Phasor(z.real,z.imag)

class Phasor(complex):
    def __init__(self, real, imag):
        self.r, self.th_r = cmath.polar(self.real+1j*self.imag)
        self.phi = self.th_r*180/cmath.pi #Degrees
        self.rect = self.real + 1j*self.imag
        self.polar = (self.r, self.phi)
    
    def __repr__(self):
        r = round(self.r,4)
        phi = round(self.phi,2)
        return "{}/_{}{}".format(r, phi, u'\xb0')
    
    def __round__(self, n=2):
        return round(self.r,n), round(self.phi,n)

    @staticmethod
    def to_polar(z):
        z = cmath.polar(z.real+1j*z.imag)
        return (z[0], z[1]*180/cmath.pi)
    @staticmethod
    def rect2polarstr(z):
        z = Phasor.to_polar(z)
        r = round(z[0],4)
        phi = round(z[1],2)
        return "{}/_{}{}".format(r, phi,u'\xb0')
    @staticmethod
    def polar_matrix_str(M):
        '''
        Takes a complex matrix and converts it to polar form (only string)
        '''
        M = np.array(M, dtype='complex_')
        Mnew = ''

        for i in range(M.shape[0]):
            Mnew += '\n'
            for j in range(M.shape[1]):
                Mnew += Phasor.rect2polarstr(M[i][j]) + '  '
        return Mnew
    @staticmethod
    def rect_matrix_str(M):
        '''
        Takes a complex matrix and converts it to rectangular form (only string)
        '''
        M = np.array(M, dtype='complex_')
        Mnew = ''

        for i in range(M.shape[0]):
            Mnew += '\n'
            for j in range(M.shape[1]):
                Mnew += str(round(M[i][j].real,4)) + '+'+ str(round(M[i][j].imag,4)) + 'j   '
        return Mnew

######### Matriz A #########
global a,A,Ainv

a = dPhasor(1,120)

A = np.array([
    [1, 1, 1],
    [1, a**2, a],
    [1, a, a**2]
])

Ainv = np.linalg.inv(A)

######### KRON REDUCTION #########

def kron_red(Y, p):
    '''
    Kron reduction algorithm on Ybus (or Zbus) for line p.
    It is expected p to be according to the numbers mathematically, not computationally,
    that means for a 3x3 matrix, p can be 1,2 or 3, and not 0,1, or 2.
    '''
    p = p-1 # So we can use 'p' in the programming convention, not mathematical
    Y = np.array(Y)
    Ynew = np.zeros(Y.shape, dtype='complex_')
    #TBD Check if it's square
    rows, columns = Y.shape
    
    for j in range(rows):
        for k in range(columns):
            if (j!=p) and (k!=p):
                Ynew[j][k] = Y[j][k]-(Y[j][p]*Y[p][k]/Y[p][p])
    # Remove p line and row
    Ynew = np.delete(Ynew,p,axis=0)
    Ynew = np.delete(Ynew,p,axis=1)
    return Ynew

def inv(Y):
    '''
    Ybus to Zbus or Zbus to Ybus
    Inverse matrix operation
    '''
    return np.linalg.inv(Y)

def get_impedance_matrix_from_Ybus(Ybus):
    Ybus = np.array(Ybus)
    Zmat = np.zeros(Ybus.shape)
    for i in range(Ybus.shape[0]):
        for j in range(Ybus.shape[1]):
            Zmat[i][j] = 1./Ybus[i][j]
    return Zmat

def get_impedance_matrix_from_Zbus(Zbus):
    Ybus = inv(Zbus)
    return get_impedance_matrix_from_Ybus(Ybus)
    

######### Zbus MODIFICATIONS #########

def Zbus_case1(Zbus, Zb, steps=True):
    '''
    Adicionando uma impedancia Zb de uma nova barra p para a referencia
    '''
    zb_orig = Zbus
    Zbus = np.array(Zbus, dtype=np.complex)
    rows, columns = Zbus.shape
    Zbus = np.append(Zbus,np.zeros((rows,1)),axis=1)
    Zbus = np.append(Zbus,np.zeros((1,columns+1)),axis=0)
    Zbus[-1][-1] = Zb
    if steps:
        msg = "CASO 1 - Adicionando uma impedancia Zb de uma nova barra p para a referencia\n"
        msg += "Colocar coluna e linha de zeros e Zb=%(Zb)s no final\n"
        msg += "Zbus original:" + Phasor.polar_matrix_str(zb_orig) + '\n\n'
        msg += "Zbus modificado:" + Phasor.polar_matrix_str(Zbus)
        print(msg % locals())
    return Zbus

def Zbus_case2(Zbus, Zb, k, steps=True, case3=False):
    '''
    Adicionando Zb de uma nova barra p para uma barra k ja existente
    '''
    kb = k
    zb_orig = Zbus
    # Convert mathematical index to programming index
    k = k-1
    #Zkk + Zb
    Zkkb = Zbus[k][k] + Zb
    Zbus = np.array(Zbus, dtype=np.complex)
    Zbus = np.vstack([Zbus, Zbus[k]])
    Zbus = np.append(Zbus, Zbus[:,k].reshape(-1,1), axis=1)
    Zbus[-1][-1] = Zkkb
    if steps:
        Zbuskk = Zbus[k][k]
        if case3 == False:
            msg = "CASO 2 - Adicionando Zb de uma nova barra p para uma barra k ja existente\n"
        else:
            msg = "CASO 3 - Adicionando Zb de uma existente k para a barra de referencia\n"
        msg += "Barra k = %(kb)s\n"
        msg += "Copiar coluna e linha k para coluna e linha p. Modificar o ultimo elemento, que sera:\n"
        msg += "Zkkb = Zbus[k][k] + Zb = %(Zbuskk)s + %(Zb)s = %(Zkkb)s, onde k=%(kb)s\n\n"
        msg += "Zbus original:" + Phasor.polar_matrix_str(zb_orig) + '\n\n'
        msg += "Zbus modificado:" + Phasor.polar_matrix_str(Zbus) +'\n\n'
        print(msg % locals())
    return Zbus

def Zbus_case3(Zbus, Zb, k, steps=True):
    '''
    Adicionando Zb de uma existente k para a barra de referencia
    '''
    Zbus = Zbus_case2(Zbus, Zb, k, steps, case3=True)
    p = Zbus.shape[0]
    Zbus = kron_red(Zbus,p)
    if steps:
        msg = "Fazendo reducao de Kron em Zbus na barra p=%(p)s\n\n"
        msg += "Zbus final:" + Phasor.polar_matrix_str(Zbus)
        print(msg % locals())
    return Zbus

def Zbus_case4(Zbus, Zb, j, k, steps=True):
    '''
    Adicionando Zb entre duas barras j e k existentes
    '''
    # Convert mathematical index to programming index
    zb_orig = Zbus
    assert(k>0 and j>0)
    j -= 1
    k -= 1
    rowj = Zbus[j]
    rowk = Zbus[k]
    subjk = np.subtract(rowj,rowk)
    Zthjk = Zbus[j][j] + Zbus[k][k] - 2*Zbus[j][k]
    Zbus = np.vstack([Zbus, subjk])
    subjk = np.r_[subjk,0]
    Zbus = np.append(Zbus, subjk.reshape(-1,1), axis=1)
    Zbus[-1][-1] = Zthjk + Zb
    p = Zbus.shape[0]
    Zbus = kron_red(Zbus,p)
    if steps:
        Zbusjj  = Phasor.rect2polarstr(Zbus[j][j])
        Zbuskk = Phasor.rect2polarstr(Zbus[k][k])
        Zbusjk_2 = Phasor.rect2polarstr(2*Zbus[j][k])
        Zthjk = Phasor.rect2polarstr(Zthjk)
        Zthb = Phasor.rect2polarstr(Zbus[-1][-1])
        jb = j+1
        kb = k+1
        msg = "CASO 4 - Adicionando Zb entre duas barras j e k existentes\n"
        msg += "Adicionando Zb=%(Zb)s entre barras %(jb)s e %(kb)s\n"
        msg += "Copiando linha j - linha k e coluna j - coluna k ao final da matriz\n"
        msg += "O ultimo elemento sera: Zthjk + Zb\nZthjk = Zbus[j][j] + Zbus[k][k] - 2*Zbus[j][k]\n"
        msg += "Zthjk = %(Zbusjj)s + %(Zbuskk)s + %(Zbusjk_2)s = %(Zthjk)s\n"
        msg += "Zthjk + Zb = %(Zthjk)s + %(Zb)s = %(Zthb)s\n\n"
        msg += "Zbus original:" + Phasor.polar_matrix_str(zb_orig) + '\n\n'
        msg += "Zbus final:" + Phasor.polar_matrix_str(Zbus) +"\n\n"
        print(msg % locals())
    return Zbus

###### FALTAS SIMETRICAS ######
def triphase_fault_zbus(Zbus, k_fault, Vf, steps=True):
    '''
    3-phase fault (Symmetrical)
    :input:
        Zbus
        k_fault - Which bus has a fault
        Vf - Pre-fault voltage
    :output:
    '''
    Zbus = np.array(Zbus, dtype='complex_')
    # math index -> prog index
    k = k_fault - 1
    If = Vf/Zbus[k][k]

    V = np.zeros(Zbus.shape[0],dtype='complex_')
    I = np.zeros(Zbus.shape,dtype='complex_')

    for j in range(len(V)):
        V[j] = Vf - Zbus[j][k]*If

    # Impedance Matrix - Not to be confused with Zbus!
    Zmat = get_impedance_matrix_from_Ybus(Zbus)
    for j in range(len(V)):
        for i in range(len(V)):
            I[i][j] = (V[i] - V[j])/Zmat[i][j]

    if steps:
        k = k_fault
        msg = "\nFalta trifasica na barra k=%(k)s com tensao pre-falta Vf de %(Vf)s pu\n"
        msg += "Zbus =\n%(Zbus)s \n\n"
        msg += "Calculando If:\n"
        msg += "If = Vf/Zbus[%(k)s][%(k)s] = %(Vf)s/" + str(round(Zbus[k][k],4)) +\
            " = " + str(Phasor.rect2polarstr(If)) + "\n\n"
        msg += "Calculando as tensoes\n"
        msg += "V[j] = Vf - Zbus[j][%(k)s]*If \n\n"
        msg += "Calculando as correntes:\n"
        msg += "I[i][j] = (V[i] - V[j])/Zb[i][j]\n\n"
        msg += "Tensoes finais:\n"
        msg += str(V.reshape(-1,1)) + "\n\n"
        msg += "Correntes finais: \n%(I)s\n"
        
        print(msg % locals())
    return V, I

###### FALTAS ASSIMETRICAS ######
def fault_calculate_voltage(Zkk1, Zkk2, Zkk0,vIfn, Vf=1, steps=True):
    Ifa0 = vIfn[0]
    Ifa1 = vIfn[1]
    Ifa2 = vIfn[2]

    Vka0 = -Zkk0*Ifa0
    Vka1 = Vf - Zkk1*Ifa1
    Vka2 = -Zkk2*Ifa2

    Vkan = np.array([Vka0, Vka1, Vka2]).reshape(-1,1)
    Vabc = (A @ Vkan) # Correcao de base
    # print("Vkan\n",Vkan)
    # print("A * Vkan\n",A @ Vkan)
    # print("Vabc\n",Vabc)
    # input()
    if steps:
        msg = "\n\nCalculando tensoes pos falta:\n"
        msg += "Vka0 = -Zkk0*Ifa0\nVka1 = Vf - Zkk1*Ifa1\nVka2 = -Zkk2*Ifa2\n\n"
        msg += "Vkan = " + Phasor.polar_matrix_str(Vkan) + '\n\n'
        msg += "Vabc = (A @ Vkan) / sqrt(3) = " + Phasor.polar_matrix_str(Vabc)
    else:
        msg = ''
    return (Vabc, Vkan, msg)

def fault_phase_gnd(Zkk1, Zkk2, Zkk0, Zf=0, Vf=1, steps=True):
    '''
    Phase-gnd fault
    Vka = 3*Zf*Ifa0
    Ifa0 = Ifa1 = Ifa2= Vf/(Zkk0+Zkk1+Zkk2+3*Zf)
    >>> nepo.fault_phase_gnd(0.2j,0.2j,0.05j,0.08j)
    '''
    Ifa0 = Vf/(Zkk0+Zkk1+Zkk2+3*Zf)
    Ifa = 3*Ifa0
    
    Ifa0, Ifa1, Ifa2 = Ifa0, Ifa0, Ifa0
    vIfn = np.array([Ifa0, Ifa1, Ifa2]).reshape(-1,1)
    # vIf = np.array(A @ vIfn).reshape(-1,1)
    vIf = np.matmul(A, vIfn)

    Vabc, Vkan, msgV = \
        fault_calculate_voltage(Zkk1, Zkk2, Zkk0, vIfn, Vf, steps=True)
    ret = vIf, vIfn, Vabc, Vkan
    if steps:
        Zff = 3*Zf
        Vf = Phasor.rect2polarstr(Vf)
        Ifa0 = Phasor.rect2polarstr(Ifa0)
        Ifa = Phasor.rect2polarstr(Ifa)
        msg = "\nFalta Fase-Terra com Vf=%(Vf)spu, Zkk1=%(Zkk1)s, Zkk2=%(Zkk2)s, Zkk0=%(Zkk0)s e Zf=%(Zf)s\n"
        msg += "\nCalculando Ifa0, Ifa1, Ifa2:\n"
        msg += "Ifa0 = Ifa1 = Ifa2 = Vf/(Zkk0+Zkk1+Zkk2+3*Zf)\n"
        msg += "Ifa0 = %(Vf)s/(%(Zkk0)s+%(Zkk1)s+%(Zkk2)s+%(Zff)s)\n"
        msg += "Ifa0 = %(Ifa0)s\n\n"
        msg += "Calculando Ifa:\n Ifa = Ifa0 + Ifa1 + Ifa2 = 3*Ifa0 \n Ifa = %(Ifa)s\n\n"
        msg += "vIfn =\n %(vIfn)s \n\n vIf =\n%(vIf)s \n\n"
        msg += msgV
        print(msg % locals())

    return ret

def fault_phase_phase(Zkk1, Zkk2, Zkk0=0, Zf=0, Vf=1, steps=True):
    '''
    Phase-Phase fault
    Vka1 - Vka2 = Ifa1*Zf
    Ifa = 0
    Ifa1=-Ifa2
    >>> nepo.fault_phase_phase(0.15j, 0.15j, 0, 0)
    '''
    #Zkk0 unused
    Ifa0 = 0
    Ifa1 = Vf/(Zkk1+Zkk2+Zf)
    Ifa2 = -1*Ifa1

    Ifa = Ifa0 + Ifa1 + Ifa2

    vIfn = np.array([Ifa0, Ifa1, Ifa2]).reshape(-1,1)
    vIf = np.array(A @ vIfn).reshape(-1,1)
    
    Vabc, Vkan, msgV = \
        fault_calculate_voltage(Zkk1, Zkk2, Zkk0,vIfn, Vf, steps=True)
    ret = vIf, vIfn, Vabc, Vkan

    if steps:
        Vf = Phasor.rect2polarstr(Vf)
        Ifa0 = Phasor.rect2polarstr(Ifa0)
        Ifa1 = Phasor.rect2polarstr(Ifa1)
        Ifa2 = Phasor.rect2polarstr(Ifa2)
        Ifa = Phasor.rect2polarstr(Ifa)
        
        msg = "\nFalta Fase-Fase com Vf=%(Vf)spu, Zkk1=%(Zkk1)s, Zkk2=%(Zkk2)s e Zf=%(Zf)s\n"
        msg += "\nCalculando Ifa0, Ifa1, Ifa2:\n"
        msg += "Ifa0 = 0\n"
        msg += "Ifa1 = -Ifa2 = Vf/(Zkk1+Zkk2+Zf)\n"
        msg += "Ifa1 = %(Vf)s/(%(Zkk1)s+%(Zkk2)s+%(Zf)s)\n"
        msg += "Ifa1 = %(Ifa1)s = -Ifa2\n"
        msg += "Ifa2 = %(Ifa2)s\n\n"
        msg += "Calculando Ifa:\n Ifa = Ifa0 + Ifa1 + Ifa2 \n Ifa = %(Ifa)s\n\n"
        msg += "vIf ="+Phasor.polar_matrix_str(vIf)+"\n\n"
        msg += "vIfn =" + Phasor.polar_matrix_str(vIfn)
        msg += msgV
        print(msg % locals())

    return ret

def fault_phase_phase_gnd(Zkk1, Zkk2, Zkk0, Zf=0, Vf=1, steps=True):
    '''
    Vka1 = Vka2 = Vka0 - 3*Zf*Ifa0
    Ifa0 + Ifa1 + Ifa2 = 0
    Ifa1 = Vf/(Zkk1 + (Zkk2*(Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf)))
    Ifa2 = -1*Ifa1*( (Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf) )
    Ifa0 = -1*Ifa1*( (Zkk2)/(Zkk2 + Zkk0 + 3*Zf) )
    '''
    Ifa1 = Vf/(Zkk1 + (Zkk2*(Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf)))
    Ifa2 = -1*Ifa1*( (Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf) )
    Ifa0 = -1*Ifa1*( (Zkk2)/(Zkk2 + Zkk0 + 3*Zf) )

    vIfn = np.array([Ifa0, Ifa1, Ifa2]).reshape(-1,1)
    vIf = np.array(A @ vIfn).reshape(-1,1)

    # Vka0 = -Zkk0*Ifa0
    # Vka1 = Vf - Zkk1*Ifa1
    # Vka2 = -Zkk2*Ifa2

    # Vkan = np.array([Vka0, Vka1, Vka2]).reshape(-1,1)
    # Vabc = (A @ Vkan)/np.sqrt(3) # Correcao de base
    Vabc, Vkan, msgV = \
        fault_calculate_voltage(Zkk1, Zkk2, Zkk0,vIfn, Vf, steps=True)
    ret = vIf, vIfn, Vabc, Vkan

    if steps:
        Ifa2 = Phasor.rect2polarstr(Ifa2)
        Ifa1 = Phasor.rect2polarstr(Ifa1)
        Ifa0 = Phasor.rect2polarstr(Ifa0)
        Zff = 3*Zf
        Vf = Phasor.rect2polarstr(Vf)
        msg = "\nFalta Fase-Fase-Terra com Vf=%(Vf)spu, Zkk1=%(Zkk1)s, Zkk2=%(Zkk2)s, Zkk0=%(Zkk0)s e Zf=%(Zf)s\n"
        msg += "\nCalculando Ifa0, Ifa1, Ifa2:\n"
        msg += "Ifa1 = Vf/(Zkk1 + (Zkk2*(Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf)))\n"
        msg += "Ifa1 = %(Vf)s/(%(Zkk1)s + (%(Zkk2)s*(%(Zkk0)s+%(Zff)s)/(%(Zkk2)s+%(Zkk0)s+%(Zff)s)))\n"
        msg += "Ifa1 = %(Ifa1)s\n\n"
        msg += "Ifa2 = -1*Ifa1*( (Zkk0+3*Zf)/(Zkk2+Zkk0+3*Zf) )\n"
        msg += "Ifa2 = -1*%(Ifa1)s*( (%(Zkk0)s+%(Zff)s)/(%(Zkk2)s+%(Zkk0)s+%(Zff)s) )\n"
        msg += "Ifa2 = %(Ifa2)s\n\n"
        msg += "Ifa0 = -1*Ifa1*( (Zkk2)/(Zkk2 + Zkk0 + 3*Zf) )\n"
        msg += "Ifa0 = -1*%(Ifa1)s*( (%(Zkk2)s)/(%(Zkk2)s + %(Zkk0)s + %(Zff)s) )\n"
        msg += "Ifa0 = %(Ifa0)s\n\n"
        msg += "Calculando Ifa, Ifb, Ifc:\n"
        msg += "Ifn =" + Phasor.polar_matrix_str(vIfn) +"\n\n[Ifa, Ifb, Ifc]' = A @ Ifn\n"
        msg += "[Ifa, Ifb, Ifc]' ="+ Phasor.polar_matrix_str(vIf)
        msg += msgV
        print(msg % locals())

    return ret

