from cmath import *
from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt

#frequência (Hz)
f = 60

#tensão eficaz (V)
V1 = 3

#frequência angular (rad/s)
w = 2*pi*f

#dispersão
k = 0.2

#RC
Rcarga = 10

#Rdc
Rdc = 0.5


def CalcularCorrentesEmSerie(Uf, Rdc, Rc, XC, XL, XM):
    Z=np.array([[Rdc+XL+XC, -XM],[-XM, XL+XC+Rdc+Rc]])
    V=np.array([Uf,0])
    i=np.dot(linalg.inv(Z),V)
    return i[0], i[1]

def CalcularCorrentesEmParalelo(Uf, Rdc, Rc, XC, XL, XM):
    Z=np.array([[Rdc+XL+XC, -XM],[-XM, (XL+Rdc)+((XC*Rc)/(XC+Rc))]])
    V=np.array([Uf,0])
    i=np.dot(linalg.inv(Z),V)
    return i[0], i[1]


lista_capacitancias = [150*10**(-9), 0.1*10**(-6), 0.47*10**(-6), 10**(-6), 4.7*10**(-6)]

# lista_eficiencias_paralelo = []
# lista_frequencia_paralelo = []
# lista_V2_paralelo = []
# lista_Rf_paralelo = []

# lista de frequências de ressonância(de 0 a 100kHz)
frequencias = np.linspace(0, 100000, 10000)
# Em série
for C in lista_capacitancias:
    print(f'Para a capacitância de {C} μF:')
    C1 = C
    XC = 1/(1j*w*C1)
    
    lista_eficiencia = []
    lista_V2 = []
    for f in frequencias:
        L = ((1/(2*pi*f))**2)*(1/C)
        XL = 1j*w*L
        M = k*L
        XM = 1j*w*M

        i1, i2 = CalcularCorrentesEmSerie(V1, Rcarga, XC, XL, XM)

        V2 = np.abs(i2*Rcarga)
        lista_V2.append(V2)

        Pot_saida = np.real(V2*i2.conjugate()*(1/2))
        Pot_entrada = np.real(V1*i1.conjugate()*(1/2))

        lista_eficiencia.append(Pot_saida/Pot_entrada)
    plt.plot(frequencias, lista_eficiencia, label=f'{C} μF')
    plt.xlabel('Frequência (Hz)')
    plt.ylabel('Eficiência')
    plt.legend()
plt.show()
    # print('Com resistor na saída de %.2f Ω:' %Rcarga)
    # print('i1 (eficaz) = %.3f A' %np.abs(i1))
    # print('i2 (eficaz) = %.3f A' %np.abs(i2))
    # print(f'V1 (eficaz): {V1} V')
    # print(f'V2 (eficaz): {V2} V')

    # print('Potência de saída: %.3f W' %Pot_saida)
    # print('Potência de entrada: %.3f W' %Pot_entrada)
    # print('Eficiência: %.3f' %Eficiencia)

    # Em paralelo
    # i1, i2 = CalcularCorrentesEmParalelo(Ufonte, Rcarga, XC1, XC2)
    # i1_conjugada = np.conj(i1)
    # i2_conjugada = np.conj(i2)
    # V1 = Ufonte
    # V2 = np.abs(i2*((XC2*Rcarga)/(XC2+Rcarga)))
    # P_saida = np.real(V2*i2.conjugate()*(1/2))
    # P_entrada = np.real(V1*i1.conjugate()*(1/2))
    # Eficiencia = Pot_saida/Pot_entrada
    # frequencia = 1 / (2 * pi * sqrt(L1 * C1))
    # Rf = Rcarga + (Rcarga*frequencia)/(10**5)

    # lista_eficiencias_paralelo.append(Eficiencia)
    # lista_frequencia_paralelo.append(frequencia)
    # lista_V2_paralelo.append(V2)
    # lista_Rf_paralelo.append(Rf)
    # print('Com resistor na saída de %.2f Ω:' %Rcarga)
    # print('i1 (eficaz) = %.3f A' %np.abs(i1))
    # print('i2 (eficaz) = %.3f A' %np.abs(i2))
    # print(f'V1 (eficaz): {V1} V')
    # print(f'V2 (eficaz): {V2} V')

    # print('Potência de saída: %.3f W' %Pot_saida)
    # print('Potência de entrada: %.3f W' %Pot_entrada)
    # print('Eficiência: %.3f' %Eficiencia)

    # print('--------------------------------------------------------')

    # print('Relação N1/N2: %.3f' %(np.abs(sqrt(L1/L2))))
    # print('Relação i2/i1: %.3f' %(np.abs(i2)/np.abs(i1)))
    # print('Relação V1/V2: %.3f' %(np.abs(Ufonte)/np.abs(i2*Rcarga)))
