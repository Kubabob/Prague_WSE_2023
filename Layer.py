from numpy import deg2rad, arange, angle
from Methods import *
import matplotlib.pyplot as plt
import pandas as pd

class Layer:
    def __init__(self, complex_refractive_index: list, incidence_angle: float, wave_length: float, absorbing: bool, nk_file: pd.DataFrame = pd.DataFrame(), thickness: float = None, roughness: bool = False) -> None:
        '''
        incidence_angle
            degrees

        wave_length
            micrometers

        thickness
            nanometers
        '''
        self.complex_refractive_indexes = [N(1,0), N(complex_refractive_index[0], complex_refractive_index[1])]

        self.incidence_angles = [deg2rad(incidence_angle), theta_j(N_i=self.complex_refractive_indexes[0],
                                         N_j=self.complex_refractive_indexes[1],
                                         theta_i=deg2rad(incidence_angle))]
        



        self.thickness = thickness
        self.wave_length = wave_length

        self.nk_file = nk_file
        self.absorbing = absorbing
        self.roughness = roughness

    
    def r_p(self):
        '''
        Returns p- reflectance of layer
        Fresnells formula
        '''

        nominator = (self.complex_refractive_indexes[1] * cos(self.incidence_angles[0]) - self.complex_refractive_indexes[0] * cos(self.incidence_angles[1]))
        denominator = (self.complex_refractive_indexes[1] * cos(self.incidence_angles[0]) + self.complex_refractive_indexes[0] * cos(self.incidence_angles[1]))

        return nominator / denominator

    def r_s(self):
        '''
        Returns s- reflectance of layer
        Fresnells formula
        '''

        nominator = (self.complex_refractive_indexes[0] * cos(self.incidence_angles[0]) - self.complex_refractive_indexes[1] * cos(self.incidence_angles[1]))
        denominator = (self.complex_refractive_indexes[0] * cos(self.incidence_angles[0]) + self.complex_refractive_indexes[1] * cos(self.incidence_angles[1]))

        return nominator / denominator
    
    def t_p(self):
        '''
        Returns p- transmitance of layer
        '''

        nominator = 2*self.complex_refractive_indexes[0]*cos(self.incidence_angles[0])
        denominator = self.complex_refractive_indexes[1]*cos(self.incidence_angles[0]) + self.complex_refractive_indexes[0]*cos(self.incidence_angles[1])

        return nominator / denominator
    
    def t_s(self):
        '''
        Returns s- transmitance of layer
        '''

        nominator = 2*self.complex_refractive_indexes[0]*cos(self.incidence_angles[0])
        denominator = self.complex_refractive_indexes[0]*cos(self.incidence_angles[0]) + self.complex_refractive_indexes[1]*cos(self.incidence_angles[1])

        return nominator / denominator
    
    
    def reflectance_plot(self):
        '''
        Shows reflectance plot for layer
        '''
        R_p_plot = []
        R_s_plot = []
        R_n_plot = []
        aranged_space = arange(0,90,0.1)
        for x in aranged_space:
            self.incidence_angles = [deg2rad(x), theta_j(N_i=self.complex_refractive_indexes[0],
                                         N_j=self.complex_refractive_indexes[1],
                                         theta_i=self.incidence_angles[-1])]
            R_p_ = R_p(self.r_p())
            R_s_ = R_s(self.r_s())
            R_p_plot.append(R_p_)
            R_s_plot.append(R_s_)
            R_n_plot.append(R_n(R_p_,R_s_))
        else:
            plt.plot(aranged_space,R_p_plot)
            plt.plot(aranged_space,R_s_plot)
            plt.plot(aranged_space,R_n_plot)
            plt.grid()
            plt.show()

    def psi(self):
        '''
        Returns psi in radians
        '''

        return arctan(abs(self.r_p()) / abs(self.r_s())) 
    
    def delta(self):
        '''
        Returns delta in radians
        '''

        delta = angle(self.r_p()/self.r_s())
        if delta < 0:
            return 2*pi + delta
        else:
            return delta

    def psi_delta_plot(self, csv_filename: str, is_wavelength: bool = False):
        df = pd.read_csv(csv_filename)
        psi_plot = []
        delta_plot = []
        if is_wavelength:
            for idx, wavelength in enumerate(df['wl']):
                wavelength = float(wavelength)
                self.complex_refractive_indexes[1] = (N(df['n'][idx], df['k'][idx]))
                psi_plot.append(rad2deg(self.psi()))
                delta_plot.append(rad2deg(self.delta()))
            else:
                wavelengths = list(map(float, df['wl']))
                plt.plot(wavelengths, psi_plot)
                plt.plot(wavelengths, delta_plot)
                plt.grid()
                plt.show()
        else:
            eV_space = []
            for idx, wavelength in enumerate(df['wl']):
                eV = wavelength_to_eV(float(wavelength))
                eV_space.append(eV)
                

