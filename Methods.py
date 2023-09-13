from numpy import cos, pi, arcsin, sin, tan, exp, power,rad2deg, sqrt, log, arctan, real, imag
import pandas as pd


def N(n, k):
    '''
    Returns complex refractive index based on n and k

    Note:
        in some cases + and - are changed due to differences in other formulas
    Default:
        -
    '''
    return n + 1j*k

def eps(N = None, e1 = None, e2 = None):
    if N:
        return N ** 2
    elif e1 and e2:
        return e1 + 1j*e2


def beta(thickness, wave_length, N_j, theta_j):
        '''
        Returns beta of jth layer
        '''
        return ((2*pi*thickness)/wave_length) * N_j * cos(theta_j)

def complex_refractive_index(n, k):
    '''
    Returns complex refractive index based on n and k

    Note:
        in some cases + and - are changed due to differences in other formulas
    Default:
        -
    '''
    return n-1j*k

def theta_j(N_i, N_j, theta_i):
    '''
    Calculates theta angle of jth layer based on Schnells formula
    N0 sin(theta0) = N1 sin(theta1)
    '''
    
    return arcsin((N_i/N_j) * sin(theta_i))


def R_p(r_p):
    '''
    Returns p- reflectance but in a representative/comparative form rather than calculational
    '''
    return abs(r_p)**2
    
def R_s(r_s):
    '''
    Returns s- reflectance but in a representative/comparative form rather than calculational
    '''
    return abs(r_s)**2

def R_n(R_p, R_s):
    '''
    Returns mean value of p- and s- comparative reflectances
    
    Shown on refractiveindex.info so maybe somewhat useful
    '''
    return (R_p + R_s)/2


def rho(psi = None, delta = None, r_p = None, r_s = None):
        '''
        Probably the most important formula here calculating rho which is tan(psi) * exp(-1i*delta)
        '''
        if psi != None and delta != None:
            return tan(psi) * exp(-1j * delta)
        elif r_p != None and r_s != None:
            return r_p / r_s
        else: raise ValueError

def wavelength_to_eV(wavelength):
    return 1239.8/wavelength

def eV_to_wavelength(eV):
    return 1239.8/eV

def Cauchy_plot(A, B, C, wavelength):
    return A + B/wavelength**2 + C/wavelength**4

def merge_lists(list_of_lists):
    merged_list = []
    for list in list_of_lists:
        merged_list.extend(list)
    return merged_list

def k_FB(A, B, C, E, Eg):
    return (A*(E-Eg)**2)/(E**2-B*E+C)

def e1_TL(A, E0, C, Eg, E, einf):
    alpha_cp = 4 * E0 ** 2 - C ** 2
    gamma_cp = E0 ** 2 - (C ** 2) / 2
    if alpha_cp <= 0 or gamma_cp <= 0 or E0 <= 0 or C <= 0:
        return 0
    alpha = sqrt(alpha_cp)
    gamma = sqrt(gamma_cp)
    if (E0**2 + Eg**2 - alpha*Eg) < 0:
        return 0
    a_ln = (Eg**2 - E0**2) * E**2 + Eg**2 * C**2 - E0**2 * (E0**2 + 3*Eg**2)
    a_atan = (E**2 - E0**2) * (E0**2 + Eg**2) + Eg**2 * C**2
    dzeta4 = (E**2 - gamma**2)**2 + (alpha**2 * C**2)/4


    return einf + ((A*C)/(pi*dzeta4)) * (a_ln/(2*alpha*E0)) * log((E0**2 + Eg**2 + alpha*Eg)/(E0**2 + Eg**2 - alpha*Eg)) - (A/(pi*dzeta4)) * (a_atan/E0)*(pi - arctan((2*Eg+alpha)/C) + arctan((-2*Eg + alpha)/C)) + 2*((A*E0)/(pi*dzeta4*alpha)) * Eg*(E**2 - gamma**2)*(pi + 2*arctan(2*(gamma**2 - Eg**2)/(alpha*C))) - ((A*E0*C)/(pi*dzeta4)) * ((E**2 + Eg**2)/E) * log(abs(E-Eg)/(E+Eg)) + (2*(A*E0*C)/(pi*dzeta4)) *Eg * log((abs(E-Eg)*(E+Eg))/(sqrt((E0**2 - Eg**2)**2 + Eg**2 * C**2)))


def e2_TL(A, E0, C, Eg, E):
    if E > Eg:
        return ((A*E0*C*power(E-Eg, 2))/(power(power(E, 2) - power(E0, 2), 2) + power(C, 2)*power(E, 2))) * 1/E
    else:
        return 0
    
def n_TL(e1, e2):
    return sqrt((e1 + sqrt(e1**2 + e2**2))/2)

def k_TL(e1, e2):
    return sqrt((-e1 + sqrt(e1**2 + e2**2))/2)
    

def e2_to_k(e2, n):
    return e2/2*n

def EMA(ea, eb, q = 1/3, fa = 0.5):
    k = (1-q)/q
    try:
        return (ea*eb + k*ea*(fa*ea+(1-fa)*eb))/(k*ea+(fa*eb+(1-fa)*ea))
    except:
        return 1


def case_insensitive_pick(dictionary: dict, possible_picks: list):
    for key in dictionary.keys():
        if key.lower() in possible_picks:
            return key

def CompleteEASE_to_normal(filepath, sep):
    file = pd.read_csv(filepath, sep=sep)
    for delta in file[case_insensitive_pick(file, ['delta'])]:
        if delta < 0:
            file = file.replace(delta, rad2deg(2*pi)+delta)
    return file

def below_zero_range(base, subtractional, minimal_value):
    if (base - subtractional) < 0:
        return minimal_value