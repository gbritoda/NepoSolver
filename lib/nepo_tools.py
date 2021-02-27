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
        phi = round(self.phi)
        return "{} /_{}°".format(r, phi)
    
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
        phi = round(z[1])
        return "{} /_{}°".format(r, phi)

def cmatrix2polar(M):
    '''
    Takes a complex matrix and converts it to polar form (only printout)
    >>> nepo.cmatrix2polar(zbus3)
    '''
    M = np.array(M, dtype='complex_')
    Mnew = ''

    for i in range(M.shape[0]):
        Mnew += '\n'
        for j in range(M.shape[1]):
            Mnew += Phasor.rect2polarstr(M[i][j]) + '  '
    print(Mnew)

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

######### Zbus MODIFICATIONS #########

def Zbus_case1(Zbus, Zb):
    '''
    Adicionando uma impedância Zb de uma nova barra p para a referencia
    '''
    Zbus = np.array(Zbus, dtype=complex)
    rows, columns = Zbus.shape
    Zbus = np.append(Zbus,np.zeros((rows,1)),axis=1)
    Zbus = np.append(Zbus,np.zeros((1,columns+1)),axis=0)
    Zbus[-1][-1] = Zb
    return Zbus

def Zbus_case2(Zbus, Zb, k):
    '''
    Adicionando Zb de uma nova barra p para uma barra k ja existente
    '''
    # Convert mathematical index to programming index
    k = k-1
    #Zkk + Zb
    Zkkb = Zbus[k][k] + Zb
    Zbus = np.array(Zbus, dtype=complex)
    Zbus = np.vstack([Zbus, Zbus[k]])
    Zbus = np.append(Zbus, Zbus[:,k].reshape(-1,1), axis=1)
    Zbus[-1][-1] = Zkkb
    return Zbus

def Zbus_case3(Zbus, Zb, k):
    '''
    Adicionando Zb de uma existente k para a barra de referencia
    '''
    Zbus = Zbus_case2(Zbus,Zb,k)
    p = Zbus.shape[0]
    Zbus = kron_red(Zbus,p)

    return Zbus

def Zbus_case4(Zbus, Zb, j, k):
    '''
    Adicionando Zb entre duas barras j e k existentes
    '''
    # Convert mathematical index to programming index
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

    for j in range(len(V)):
        for i in range(len(V)):
            I[i][j] = (V[i] - V[j])/Zbus[i][j]

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
        msg += "I[i][j] = (V[i] - V[j])/Zbus[i][j]\n\n"
        msg += "Tensoes finais:\n"
        msg += str(V.reshape(-1,1)) + "\n\n"
        msg += "Correntes finais: \n%(I)s\n"
        
        print(msg % locals())
    return V, I

###### FALTAS ASSIMETRICAS ######
def fault_phase_gnd(Zkk1, Zkk2, Zkk0, Zf, Vf=1, steps=True):
    '''
    Phase-gnd fault
    Ifa0 = Ifa1 = Ifa2= Vf/(Zkk0+Zkk1+Zkk2+3*Zf)
    >>> nepo.fault_phase_gnd(0.2j,0.2j,0.05j,0.08j)
    '''
    Ifa0 = Vf/(Zkk0+Zkk1+Zkk2+3*Zf)
    Ifa = 3*Ifa0
    
    Ifa0, Ifa1, Ifa2 = Ifa0, Ifa0, Ifa0
    vIfn = np.array([Ifa0, Ifa1, Ifa2]).reshape(-1,1)
    vIf = np.array(A @ vIfn).reshape(-1,1)

    ret = vIf, vIfn
    if steps:
        Zff = 3*Zf
        Vf = Phasor.rect2polarstr(Vf)
        Ifa0 = Phasor.rect2polarstr(Ifa0)
        Ifa = Phasor.rect2polarstr(Ifa)
        msg = "\nCalculando Ifa0, Ifa1, Ifa2:\n"
        msg += "Ifa0 = Ifa1 = Ifa2 = Vf/(Zkk0+Zkk1+Zkk2+3*Zf)\n"
        msg += "Ifa0 = %(Vf)s/(%(Zkk0)s+%(Zkk1)s+%(Zkk2)s+%(Zff)s)\n"
        msg += "Ifa0 = %(Ifa0)s\n\n"
        msg += "Calculando Ifa:\n Ifa = Ifa0 + Ifa1 + Ifa2 = 3*Ifa0 \n Ifa = %(Ifa)s\n\n"
        msg += "vIfn =\n %(vIfn)s \n\n vIf =\n%(vIf)s \n\n"
    print(msg % locals())

    return ret

def fault_phase_phase(Zkk1, Zkk2, Zf, Vf=1, steps=True):
    '''
    Phase-Phase fault
    Ifa = 0
    Ifa1=-Ifa2
    >>> nepo.fault_phase_phase(0.15j, 0.15j, 0)
    '''
    Ifa0 = 0
    Ifa1 = Vf/(Zkk1+Zkk2+Zf)
    Ifa2 = -1*Ifa1

    Ifa = Ifa0 + Ifa1 + Ifa2

    vIfn = np.array([Ifa0, Ifa1, Ifa2]).reshape(-1,1)
    vIf = np.array(A @ vIfn).reshape(-1,1)
    
    ret = vIf, vIfn
    if steps:
        Vf = Phasor.rect2polarstr(Vf)
        Ifa0 = Phasor.rect2polarstr(Ifa0)
        Ifa1 = Phasor.rect2polarstr(Ifa1)
        Ifa2 = Phasor.rect2polarstr(Ifa2)
        Ifa = Phasor.rect2polarstr(Ifa)

        msg = "\nCalculando Ifa0, Ifa1, Ifa2:\n"
        msg += "Ifa0 = 0\n"
        msg += "Ifa1 = -Ifa2 = Vf/(Zkk1+Zkk2+Zf)\n"
        msg += "Ifa1 = %(Vf)s/(%(Zkk1)s+%(Zkk2)s+%(Zf)s)\n"
        msg += "Ifa1 = %(Ifa1)s = -Ifa2\n"
        msg += "Ifa2 = %(Ifa2)s\n\n"
        msg += "Calculando Ifa:\n Ifa = Ifa0 + Ifa1 + Ifa2 \n Ifa = %(Ifa)s\n\n"
        msg += "vIf =\n %(vIf)s \n\n vIfn =\n%(vIfn)s \n\n"
    print(msg % locals())

    return ret

def fault_phase_phase_gnd():
    print('to be done')
    return 0
