from Layer import Layer
from Methods import *
from numpy import matrix, matmul, identity, arctan, angle, arange, rad2deg, deg2rad
import matplotlib.pyplot as plt
import pandas as pd

class Multilayer:
    def __init__(self, layers: list[Layer]) -> None:
        self.layers = layers

        self.complex_refractive_indexes = [N(1,0)]

        self.layers_nk = []

        for layer in self.layers:
            self.complex_refractive_indexes.append(layer.complex_refractive_indexes[1])
            if layer.nk_file.empty:
                self.layers_nk.append(None)
            else:
                self.layers_nk.append(layer.nk_file)

        
        self.incidence_angles = [layers[0].incidence_angles[0]]

        for i in range(len(self.complex_refractive_indexes)-1):
            self.incidence_angles.append(theta_j(N_i=self.complex_refractive_indexes[i],
                                                 N_j=self.complex_refractive_indexes[i+1],
                                                 theta_i=self.incidence_angles[-1]))
            
        self.thicknesses = [layer.thickness for layer in self.layers]

        self.wave_length = self.layers[0].wave_length
        
            
    def jth_r_s(self, number_of_layer):
        '''
        Returns s- reflectance of jth and j+1th layer
        '''
        nominator = (self.complex_refractive_indexes[number_of_layer] * cos(self.incidence_angles[number_of_layer]) - self.complex_refractive_indexes[number_of_layer+1] * cos(self.incidence_angles[number_of_layer+1]))
        denominator = (self.complex_refractive_indexes[number_of_layer] * cos(self.incidence_angles[number_of_layer]) + self.complex_refractive_indexes[number_of_layer+1] * cos(self.incidence_angles[number_of_layer+1]))

        return nominator / denominator
    
    def jth_r_p(self, number_of_layer):
        '''
        Returns p- reflectance of jth and j+1th layer
        '''
        nominator = (self.complex_refractive_indexes[number_of_layer+1] * cos(self.incidence_angles[number_of_layer]) - self.complex_refractive_indexes[number_of_layer] * cos(self.incidence_angles[number_of_layer+1]))
        denominator = (self.complex_refractive_indexes[number_of_layer+1] * cos(self.incidence_angles[number_of_layer]) + self.complex_refractive_indexes[number_of_layer] * cos(self.incidence_angles[number_of_layer+1]))

        return nominator / denominator
    
    def jth_t_s(self, number_of_layer):
        '''
        Returns s- transmitance of layer
        '''

        nominator = 2*self.complex_refractive_indexes[number_of_layer]*cos(self.incidence_angles[number_of_layer])
        denominator = self.complex_refractive_indexes[number_of_layer]*cos(self.incidence_angles[number_of_layer]) + self.complex_refractive_indexes[number_of_layer+1]*cos(self.incidence_angles[number_of_layer+1])

        return nominator / denominator
    
    def jth_t_p(self, number_of_layer):
        '''
        Returns p- transmitance of layer
        '''

        nominator = 2*self.complex_refractive_indexes[number_of_layer]*cos(self.incidence_angles[number_of_layer])
        denominator = self.complex_refractive_indexes[number_of_layer+1]*cos(self.incidence_angles[number_of_layer]) + self.complex_refractive_indexes[number_of_layer]*cos(self.incidence_angles[number_of_layer+1])

        return nominator / denominator
    

            
    def I_matrix(self, number_of_layer, s = None, p = None):
        '''
        Returns transfer matrix between jth and j+1th
        '''
        if s:
            r_s = self.jth_r_s(number_of_layer)
            return matrix([[1, r_s], [r_s, 1]]) / self.jth_t_s(number_of_layer)
        elif p:
            r_p = self.jth_r_p(number_of_layer)
            return matrix([[1, r_p], [r_p, 1]]) / self.jth_t_p(number_of_layer)
        else:
            raise ValueError

    
    def L_matrix(self, number_of_layer):
        return matrix([[exp(-1j*beta(self.thicknesses[number_of_layer-1], self.wave_length, self.complex_refractive_indexes[number_of_layer], self.incidence_angles[number_of_layer])), 0], [0, exp(1j*beta(self.thicknesses[number_of_layer-1], self.wave_length, self.complex_refractive_indexes[number_of_layer], self.incidence_angles[number_of_layer]))]])


    def S_matrix(self, s = None, p = None):
        S = identity(2)
        if s:
            S = matmul(S, self.I_matrix(0, s=True))
            for i in range(1, len(self.layers)):
                S = matmul(S, self.L_matrix(i))
                S = matmul(S, self.I_matrix(i, s=True))
            else:
                return S
        elif p:
            S = matmul(S, self.I_matrix(0, p=True))
            for i in range(1, len(self.layers)):
                S = matmul(S, self.L_matrix(i))
                S = matmul(S, self.I_matrix(i, p=True))
            else:
                return S
        else:
            raise ValueError
        
    def rs(self):
        Ss = self.S_matrix(s=True)

        rs = Ss[1,0]/Ss[0,0]
        return rs
    
    def rp(self):
        Sp = self.S_matrix(p=True)

        rp = Sp[1,0]/Sp[0,0]
        return rp
        
    def rho(self):
        '''
        Returns rho based on transfer matrix
        '''

        return self.rp()/self.rs()
    
    def psi(self, degrees = None, radians = None):
        '''
        Returns psi in radians
        '''


        if degrees:
            return rad2deg(arctan(abs(self.rp()) / abs(self.rs())))
        elif radians:
            return arctan(abs(self.rp()) / abs(self.rs()))
        else:
            raise ValueError
    
    def delta(self, degrees = None, radians = None):
        '''
        Returns delta in radians
        '''


        delta = -angle(self.rp()/self.rs())
        if degrees:
            if delta < 0:
                return rad2deg(2*pi + delta)
            else:
                return rad2deg(delta)
        elif radians:
            if delta < 0:
                return 2*pi + delta
            else:
                return delta
        else:
            raise ValueError
        
    def psi_delta_plot(self, reference_paths: list[str], eV = None, wave_length = None, degrees = None, radians = None):
        y_psi = []
        y_delta = []
        references = []

        for filepath in reference_paths:
            df = pd.read_csv(filepath, sep='\t')
            references.append(df)
        if eV:
            x_space = [wavelength*1000 for wavelength in references[0]['wvl']]
            for idx, wave_length in enumerate(x_space):

                self.complex_refractive_indexes = [N(1,0)]
                for i in range(len(references)):
                    self.complex_refractive_indexes.append(N(references[i]['n'][idx], references[i]['k'][idx]))

                self.wave_length = wave_length

                if degrees:
                    y_psi.append(self.psi(degrees=True))
                    y_delta.append(self.delta(degrees=True))
                elif radians:
                    y_psi.append(self.psi(radians=True))
                    y_delta.append(self.delta(radians=True))
                else:
                    raise ValueError
                
            else:
                plt.plot(x_space, y_psi, label='psi')
                plt.plot(x_space, y_delta, label='delta')
                plt.xlabel('photon energy[eV](more like wave_length in nm needs to be fixed)')
                plt.grid()
                plt.show()
        elif wave_length:

            x_space = [wavelength*1000 for wavelength in references[0]['wvl']] #need to be fixed, can't be hard coded
            for idx, wave_length in enumerate(x_space):

                self.complex_refractive_indexes = [N(1,0)]
                for i in range(len(references)):
                    self.complex_refractive_indexes.append(N(references[i]['n'][idx], references[i]['k'][idx]))

                self.wave_length = wave_length

                if degrees:
                    y_psi.append(self.psi(degrees=True))
                    y_delta.append(self.delta(degrees=True))
                elif radians:
                    y_psi.append(self.psi(radians=True))
                    y_delta.append(self.delta(radians=True))
                else:
                    raise ValueError
                
            else:
                plt.plot(x_space, y_psi, label='psi')
                plt.plot(x_space, y_delta, label='delta')
                plt.xlabel('wave length[um]')
                plt.grid()
                plt.show()

        else:
            raise ValueError

    def reflectance_plot(self):
        y_p = []
        y_s = []
        y_n = []
        x_space = arange(0,90,0.1)
        for angle in x_space:
            self.incidence_angles = [deg2rad(angle)]
            for i in range(len(self.complex_refractive_indexes)-1):
                self.incidence_angles.append(theta_j(N_i=self.complex_refractive_indexes[i],
                                                 N_j=self.complex_refractive_indexes[i+1],
                                                 theta_i=self.incidence_angles[-1]))
            y_p.append(R_p(self.rp()))
            y_s.append(R_s(self.rs()))
            y_n.append(R_n(y_p[-1], y_s[-1]))
        else:
            plt.plot(x_space,y_p)
            plt.plot(x_space,y_s)
            plt.plot(x_space,y_n)
            plt.grid()
            plt.show()

    

