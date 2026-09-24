# -*- coding: utf-8 -*-

# =============================================================================
# MIT License
# 
# The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
# 
# CSIRO grants you a licence to the Software on the terms of the MIT Licence.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# =============================================================================

'''
    File name: femticPy.py
    Author: Hoël Seillé 
    Date created: 01/09/2023
    Date last modified: 01/07/2024
    Python Version: 3.6
'''
__author__ = "Hoël Seillé"
__copyright__ = "Copyright 2024, CSIRO"
__credits__ = ["Hoël Seillé"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Hoël Seillé"
__email__ = "hoel.seille@csiro.au"
__status__ = "Beta"




import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import os
from copy import deepcopy

import EDItools as edi



"""
This script prepares the model and data inputs required for the inversion. 
It also read the results of the inversion

It consists of 3 parts:

    - dataGen creates the input data files for the inversion, using .edi files   
    
    - meshGen prepares the input data required to run the sequential shell 
    script that creates the input mesh and .vtk visualization files

    - InvResults

Files required: 

    - dataGen:
        - .edi files
        - data coordinates
        
    - meshGen:
        - data coordinates
        - coast_line 
        - topography 
        - bathymetry 
    
"""

class DataGen():

    def __init__(self, survey, outdir):
        self.survey = survey

        
        ## Size of computational domain (in km)
        self.analysis_domain = [[-500.0, 500.0],
                                [-500.0, 500.0],
                                [-500.0, 500.0]]
        
        ## Mesh center (in km)
        self.center = [0.0, 0.0, 0.0]
        
        ## Data to be inverted for
        self.invert_Z = True
        self.invert_VTF = False
        self.invert_PT = False
        self.error_floor_Z = [0.05, 0.05, 0.05, 0.05] # [Zxx, Zxy, Zyx, Zyy]
        self.error_floor_Tz = 0.05 #(absolute value)
        self.error_floor_PT = None 
        
        
        # MT sites data
        self.mt_coords_utm = None
        # self.mt_coords = None
        #self.mt_ids = []
        self.ids_Z   = []
        self.ids_VTF  = []
        self.ids_PT  = []
        self.data_Z  = []
        self.data_VTF = []
        self.data_PT = []
        #self.mt_data = None
        #self.mt_data_inversion = None
        self.nRx = 0
        self.nFreq = 0   
        self.freqs = []

        self.topography = False
        self.bathymetry = False
        self.coast_line = False
        
        self.outdir = outdir
        
    
        
    def _filter_dataframe_list(self, df_list, max_val=1e2):
        """
        Filters each dataframe in a list, removing rows where any numeric column
        except 'FREQ' exceeds max_val in absolute value.
        """
        new_list = []
        for i, df in enumerate(df_list):
            # skip empty or None
            if df is None or len(df) == 0:
                new_list.append(df)
                continue
    
            cols = [c for c in df.columns if c != "FREQ"]
            mask = (df[cols].abs() <= max_val).all(axis=1)
    
            removed = (~mask).sum()
            if removed > 0:
                print(f"[DataGen] Station {i}: removed {removed} rows")
    
            new_list.append(df.loc[mask].reset_index(drop=True))

        return new_list
    
    def read_MTdata(self, edi_path):

        site_ids=[]
        freq_list=[]
        
        for file in os.listdir(edi_path):
            if file.endswith(".edi"):
                site_ids.append(file[:-4])
        site_ids = np.sort(site_ids)
        
        for site_id in site_ids: 
                    self.nRx += 1
                    edi_file = '%s.edi'%site_id
                    edi_file_path = '%s/%s'%(edi_path, edi_file)
                    print ('    reading %s'%edi_file)
                    #read edi file
                    edi_data, site_id, coord = edi.read(edi_file_path)
                    

                    
                    if 'ZXXR' in edi_data.columns:

                        data_Z = edi_data[['FREQ',
                                 'ZXXR','ZXXI','ZXX.VAR',
                                 'ZXYR','ZXYI','ZXY.VAR',
                                 'ZYXR','ZYXI','ZYX.VAR',
                                 'ZYYR','ZYYI','ZYY.VAR']]
                        
                        # remove masked data, uncomment below if needed
                        # masked=np.unique(np.where(data_Z[['ZXXR','ZXXI','ZXX.VAR',
                        #                               'ZXYR','ZXYI','ZXY.VAR',
                        #                               'ZYXR','ZYXI','ZYX.VAR',
                        #                               'ZYYR','ZYYI','ZYY.VAR']]>= 1e10)[0])
                        # data_Z.drop(masked, inplace=True)
                        # data_Z = data_Z.reset_index(drop=True)
                        
                        self.data_Z.append(data_Z)
                        self.ids_Z.append(site_id)
                    
                    
                    if 'TXR.EXP' in edi_data.columns:

                        data_VTF = edi_data[['FREQ',
                                 'TXR.EXP','TXI.EXP','TXVAR.EXP',
                                 'TYR.EXP','TYI.EXP','TYVAR.EXP',]]
                        
                        # remove masked data, uncomment below if needed
                        # masked=np.unique(np.where(data_VTF[[
                        #          'TXR.EXP','TXI.EXP','TXVAR.EXP',
                        #          'TYR.EXP','TYI.EXP','TYVAR.EXP',]]>= 1e10)[1])
                        
                        # data_VTF.drop(masked, inplace=True)
                        # data_VTF = data_VTF.reset_index(drop=True)
                        
                        self.data_VTF.append(data_VTF)
                        self.ids_VTF.append(site_id)
                    
                    
                    # data.append([site_id, edi_data])
                    # site_ids.append(site_id)
                    # coords.append(coord)
                    # list all the existing frequencies in the dataset
                    for freq in range(len(edi_data['FREQ'].values)):
                        if not edi_data['FREQ'].values[freq] in freq_list:
                            if freq == 0:
                                freq_list.append(edi_data['FREQ'].values[freq])
                            else:
                                if min(np.abs(np.log10(freq_list) - np.log10(edi_data['FREQ'].values[freq]))) < 0.01:
                                    idx = (np.abs(freq_list - edi_data['FREQ'].values[freq])).argmin()
                                    edi_data['FREQ'][freq] = freq_list[idx]
                                else:
                                    freq_list.append(edi_data['FREQ'].values[freq])

        freq_list=np.sort(freq_list)[::-1]

        # self.mt_ids = site_ids
        # self.mt_data = data
        self.nRx_Z = len(self.data_Z)
        self.nRx_VTF = len(self.data_VTF)
        self.nFreq = len(freq_list)
        self.freqs = freq_list
        
        print('Read %d .edi files'%self.nRx)


    def read_MTdata_coordinates(self, coordinates_path, coordinates_file):
        
        # read coordinates_file
        coords = pd.read_csv('%s/%s'%(coordinates_path,coordinates_file), 
                             header=None, sep='\s+',
                             names = ['id','east','north','z'],
                             dtype={'id': 'string',
                                    'east': 'float64',
                                    'north': 'float64',
                                    'z': 'float64'})
        
        coords = coords.sort_values('id')
        coords = coords.reset_index(drop=True) 
        
        self.mt_coords_utm = coords
        self.mt_coords_utm['z'] = - self.mt_coords_utm['z']
        self.mt_coords = self.mt_coords_utm.copy(deep=True)
        
        

    
    def center_data(self):
        """
        This function defines the coordinates of the MT stations relative to
        the mesh center
        """    
        east_min = self.mt_coords_utm['east'].min()
        east_max = self.mt_coords_utm['east'].max()
        north_min = self.mt_coords_utm['north'].min()
        north_max = self.mt_coords_utm['north'].max()
        east_center = east_min + (east_max - east_min) / 2
        north_center = north_min + (north_max - north_min) / 2
        
        self.mt_coords['east'] -= east_center
        self.mt_coords['north'] -= north_center
        self.anchor = [north_center, east_center]



    def convert_units(self, data_orig):
        # Conversion to appropriate units: FEMTIC uses ohms
        #   1 ohm = 10000(4*pi) [mV/km/nT]
        C = 10000/(4*np.pi)
        data = data_orig.copy(deep=True)
        data.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] = data_orig.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] / C
        data.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] = data_orig.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] / (C**2)
        return data 
    
    
  

    def apply_error_floor(self, data, data_type = 'Z', error_type = 2):
        
        if data_type == 'Z':
            std = data[['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']]**0.5  
            
            if error_type == 0:
                #calculate error floor per component
                

                for freq in range(len(data)):
                    error_floor = []
                    comp = ['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']
                    for i in range(4):  
                        error_floor.append(
                            ((data[comp[2*i]][freq]**2+data[comp[2*i+1]][freq]**2)**0.5) *
                            np.array(self.error_floor_Z)[i])
                        
                    for i, comp in enumerate(['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']):
                        
                        if std[comp][freq] < error_floor[i]:
                            data[comp][freq] =  (error_floor[i])**2
                            
            elif error_type == 1 :      # orig 
                #calculate error original version
                for freq in range(len(data)):
                    
                    error_floor = ((((data['ZXYR'][freq]**2+data['ZXYI'][freq]**2)**0.5) * 
                                    ((data['ZYXR'][freq]**2+data['ZYXI'][freq]**2)**0.5))**0.5)  * np.array(self.error_floor_Z)
                    
                    for i, comp in enumerate(['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']):
                        if std[comp][freq] < error_floor[i]:
                            data[comp][freq] =  (error_floor[i])**2
                
            elif error_type == 2 : 
                
                for freq in range(len(data)):
                    error_floor = []
                    comp = ['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']
                    
                    i = 0
                    error_floor.append(
                            ((((data['ZXYR'][freq]**2+data['ZXYI'][freq]**2)**0.5) * 
                              ((data['ZXXR'][freq]**2+data['ZXXI'][freq]**2) ** 0.5))**0.5) *
                            np.array(self.error_floor_Z)[i])
                        
                    for i in [1, 2]:  
                        error_floor.append(
                            ((((data['ZXYR'][freq]**2+data['ZXYI'][freq]**2)**0.5) * 
                              ((data['ZYXR'][freq]**2+data['ZYXI'][freq]**2) ** 0.5))**0.5) *
                            np.array(self.error_floor_Z)[i])

                    
                    i = 3
                    error_floor.append(
                            ((((data['ZYXR'][freq]**2+data['ZYXI'][freq]**2)**0.5) * 
                              ((data['ZYYR'][freq]**2+data['ZYYI'][freq]**2) ** 0.5))**0.5) *
                            np.array(self.error_floor_Z)[i])
                    
                    
                    for i, comp in enumerate(['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']):
                        if std[comp][freq] < error_floor[i]:
                            data[comp][freq] =  (error_floor[i])**2
                            
            else:
                sys.exit('')
            
        elif data_type == 'VTF':
            std = data[['TXVAR.EXP','TYVAR.EXP']]**0.5   
            for freq in range(len(data)):
                for i, comp in enumerate(['TXVAR.EXP','TYVAR.EXP']):
                    if std[comp][freq] < self.error_floor_Tz:
                        data[comp][freq] =  (self.error_floor_Tz)**2

        return data
  

    def apply_defined_error(self, data, data_type='Z', comps=['ZXX.VAR', 'ZXY.VAR', 'ZYX.VAR', 'ZYY.VAR'],
                            err_val=[0.15, 0.05, 0.05, 0.15], error_type=2):
        """
        Applies fixed error values to specific components using signal-amplitude-based scaling.
    
        Parameters
        ----------
        data : list of DataFrames
            One DataFrame per frequency.
        data_type : str
            'Z' for impedance or 'VTF' for vertical transfer function.
        comps : list of str
            Target components to assign error to (e.g., 'ZXX.VAR').
        err_val : list of float
            Fixed error multipliers per component.
        error_type : int
            0, 1, or 2 — determines how the error is scaled.
        """
    
        if data_type == 'Z':
            for freq in range(len(data)):
                df = data[freq]  # Current frequency DataFrame
    
                error_dict = {}  # Will hold {component: computed_error}
    
                if error_type == 0:
                    for i, comp in enumerate(comps):
                        real = comp.replace('.VAR', 'R')
                        imag = comp.replace('.VAR', 'I')
                        amp = np.sqrt(df[real]**2 + df[imag]**2)
                        error_dict[comp] = amp * err_val[i]
    
                elif error_type == 1:
                    zxy_amp = np.sqrt(df['ZXYR']**2 + df['ZXYI']**2)
                    zyx_amp = np.sqrt(df['ZYXR']**2 + df['ZYXI']**2)
                    scalar = np.sqrt(zxy_amp * zyx_amp)
                    for i, comp in enumerate(comps):
                        error_dict[comp] = scalar * err_val[i]
    
                elif error_type == 2:
                    if 'ZXX.VAR' in comps:
                        i = comps.index('ZXX.VAR')
                        error_dict['ZXX.VAR'] = (
                            np.sqrt(np.sqrt(df['ZXYR']**2 + df['ZXYI']**2) *
                                    np.sqrt(df['ZXXR']**2 + df['ZXXI']**2)) * err_val[i])
    
                    if 'ZXY.VAR' in comps:
                        i = comps.index('ZXY.VAR')
                        error_dict['ZXY.VAR'] = (
                            np.sqrt(np.sqrt(df['ZXYR']**2 + df['ZXYI']**2) *
                                    np.sqrt(df['ZYXR']**2 + df['ZYXI']**2)) * err_val[i])
    
                    if 'ZYX.VAR' in comps:
                        i = comps.index('ZYX.VAR')
                        error_dict['ZYX.VAR'] = (
                            np.sqrt(np.sqrt(df['ZXYR']**2 + df['ZXYI']**2) *
                                    np.sqrt(df['ZYXR']**2 + df['ZYXI']**2)) * err_val[i])
    
                    if 'ZYY.VAR' in comps:
                        i = comps.index('ZYY.VAR')
                        error_dict['ZYY.VAR'] = (
                            np.sqrt(np.sqrt(df['ZYXR']**2 + df['ZYXI']**2) *
                                    np.sqrt(df['ZYYR']**2 + df['ZYYI']**2)) * err_val[i])
                else:
                    raise ValueError("Invalid error_type. Must be 0, 1, or 2.")
    
                # Apply errors
                for comp, err in error_dict.items():
                    data[freq][comp] = err ** 2
    
        elif data_type == 'VTF':
            for freq in range(len(data)):
                for i, comp in enumerate(comps):
                    data[freq][comp] = err_val[i] ** 2
    
        return data


    def add_gaussian_noise(self, data, data_type='Z', seed=None):
        """
        Add Gaussian noise to MT data components.
        
        Parameters
        ----------
        data : list of DataFrames
            List of DataFrames, one per frequency.
        data_type : {"Z", "VTF"}
            Whether we process impedance or tipper.
        seed : int or None
            If given, fixes the RNG seed for reproducibility.
            
        Returns
        -------
        noisy_data : deep copy of input data
        """
        import copy
        import numpy as np
    
        if seed is not None:
            np.random.seed(seed)
    
        noisy_data = copy.deepcopy(data)
    
        # ------------------------------------------------------------------
        # Define component mapping based on data_type
        # ------------------------------------------------------------------
        if data_type == 'Z':
            comp_pairs = [
                ('ZXXR','ZXXI','ZXX.VAR'),
                ('ZXYR','ZXYI','ZXY.VAR'),
                ('ZYXR','ZYXI','ZYX.VAR'),
                ('ZYYR','ZYYI','ZYY.VAR'),
            ]
    
        elif data_type == 'VTF':
            comp_pairs = [
                ('TXR.EXP','TXI.EXP','TXVAR.EXP'),
                ('TYR.EXP','TYI.EXP','TYVAR.EXP'),
            ]
    
        else:
            raise ValueError("data_type must be 'Z' or 'VTF'.")
    
        # ------------------------------------------------------------------
        # Loop on frequencies and add Gaussian noise
        # ------------------------------------------------------------------
        for fi, df in enumerate(noisy_data):
    
            # Compute std from VAR
            std_dict = {var: np.sqrt(df[var]) for _, _, var in comp_pairs}
    
            # Apply Gaussian noise to each component
            for real, imag, var in comp_pairs:
                sigma = std_dict[var]
    
                noisy_data[fi][real] = df[real] + np.random.normal(0, sigma)
                noisy_data[fi][imag] = df[imag] + np.random.normal(0, sigma)
    
        return noisy_data



    # def apply_defined_error(self, data, data_type='Z', comps=['ZXX.VAR', 'ZXY.VAR', 'ZYX.VAR', 'ZYY.VAR'],
    #                         err_val=[0.15, 0.05, 0.05, 0.15], error_type=2):
    #     """
    #     Applies fixed error values to specific components using signal-amplitude-based scaling.
    
    #     Parameters
    #     ----------
    #     data : DataFrame
    #         Magnetotelluric data.
    #     data_type : str
    #         'Z' for impedance or 'VTF' for vertical transfer function.
    #     comps : list of str
    #         Target components to assign error to (e.g., 'ZXX.VAR').
    #     err_val : list of float
    #         Fixed error multipliers per component.
    #     error_type : int
    #         0, 1, or 2 — determines how the error is scaled.
    #     """
    
    #     if data_type == 'Z':
    #         for freq in range(len(data)):
    #             error = []
    
    #             if error_type == 0:
    #                 # Independent error per component using norm(Zij)
    #                 for i, comp in enumerate(comps):
    #                     real = comp.replace('.VAR', 'R')
    #                     imag = comp.replace('.VAR', 'I')
    #                     amp = np.sqrt(data[freq][real]**2 + data[freq][imag]**2)
    #                     error.append(amp * err_val[i])
    
    #             elif error_type == 1:
    #                 # Common scalar using ZXY and ZYX
    #                 zxy_amp = np.sqrt(data[freq]['ZXYR']**2 + data[freq]['ZXYI']**2)
    #                 zyx_amp = np.sqrt(data[freq]['ZYXR']**2 + data[freq]['ZYXI']**2)
    #                 scalar = np.sqrt(zxy_amp * zyx_amp)
    #                 error = scalar * np.array(err_val)
    
    #             elif error_type == 2:
    #                 # Cross-term scaling per pair
    #                 # 0: ZXX vs ZXY, 1-2: ZYX vs ZXY, 3: ZYY vs ZYX
    #                 if 'ZXX.VAR' in comps:
    #                     i=comps.index('ZXX.VAR')
    #                     error.append(
    #                         np.sqrt(np.sqrt((data[freq]['ZXYR']**2 + data[freq]['ZXYI']**2)) *
    #                                 np.sqrt((data[freq]['ZXXR']**2 + data[freq]['ZXXI']**2))) * err_val[i])
        
    #                 if 'ZXY.VAR' in comps:
    #                     i=comps.index('ZXY.VAR')
    #                     error.append(
    #                         np.sqrt(np.sqrt((data[freq]['ZXYR']**2 + data[freq]['ZXYI']**2)) *
    #                                 np.sqrt((data[freq]['ZYXR']**2 + data[freq]['ZYXI']**2))) * err_val[i])
                    
    #                 if 'ZYX.VAR' in comps:
    #                     i=comps.index('ZYX.VAR')
    #                     error.append(
    #                         np.sqrt(np.sqrt((data[freq]['ZXYR']**2 + data[freq]['ZXYI']**2)) *
    #                                 np.sqrt((data[freq]['ZYXR']**2 + data[freq]['ZYXI']**2))) * err_val[i])
    #                 if 'ZYY.VAR' in comps:
    #                     i=comps.index('ZYY.VAR')
    #                     error.append(
    #                         np.sqrt(np.sqrt((data[freq]['ZYXR']**2 + data[freq]['ZYXI']**2)) *
    #                                 np.sqrt((data[freq]['ZYYR']**2 + data[freq]['ZYYI']**2))) * err_val[i])
    #             else:
    #                 raise ValueError("Invalid error_type. Must be 0, 1, or 2.")
    
    #             # Apply defined error (overwrite std dev with fixed error^2)
    #             for i, comp in enumerate(comps):
    #                 data[freq][comp] = error[i] ** 2
    
    #     elif data_type == 'VTF':
    #         for freq in range(len(data)):
    #             for i, comp in enumerate(comps):
    #                 data[freq][comp] = err_val[i] ** 2
    
    #     return data

        
    
    def conjugate(self, data, data_type = 'Z'):
        # FEMTIC convention takes the complex conjugate of Z (exp(+iwt) assumed)
        
        if data_type == 'Z':
            data['ZXXI'] = -data['ZXXI']
            data['ZXYI'] = -data['ZXYI']
            data['ZYXI'] = -data['ZYXI']
            data['ZYYI'] = -data['ZYYI']
            
        elif data_type == 'VTF':
            data['TXI.EXP'] = -data['TXI.EXP']
            data['TYI.EXP'] = -data['TYI.EXP']
            
        return data
    
    
            
    def prep_data(self, freq_bandwidth = None, subsampling = 1):
        
        # subsample the original frequency list:
        if isinstance(subsampling, int):
            # subsample every th frequency
            self.freqs_inversion = self.freqs[::subsampling]
        else:
            # subsambple at the indexes specified in the list provided as subsampling
            self.freqs_inversion = self.freqs[subsampling]
        self.freqs_inversion = self.freqs_inversion[np.where((self.freqs_inversion<=freq_bandwidth[1]) & (self.freqs_inversion>=freq_bandwidth[0]))]
        
        #self.mt_data_inversion = deepcopy(self.mt_data)
        # flag for specifying that data has already been prepared ?
        
        if self.invert_Z:
        
            for rx in range(self.nRx_Z):
                
                rx_data = self.data_Z[rx]
                rx_data = rx_data[rx_data['FREQ'].isin(self.freqs_inversion)]
                #rx_data = rx_data.iloc[::subsampling, :]
                rx_data.reset_index(inplace=True, drop=True)
                
                # we convert the units 
                rx_data_ = self.convert_units(rx_data)
                # we take the complex conjugate 
                rx_data_ = self.conjugate(rx_data_,data_type = 'Z')   
                # we apply the error floors
                rx_data_ = self.apply_error_floor(rx_data_,data_type = 'Z')
    
                self.data_Z[rx] = rx_data_


        if self.invert_VTF:
        
            for rx in range(self.nRx_VTF):
                                
                rx_data = self.data_VTF[rx]
                rx_data = rx_data[rx_data['FREQ'].isin(self.freqs_inversion)]
                #rx_data = rx_data.iloc[::subsampling, :]
                rx_data.reset_index(inplace=True, drop=True)

                # we take the complex conjugate
                rx_data_ = self.conjugate(rx_data,data_type = 'VTF')   
                # we apply the error floors
                rx_data_ = self.apply_error_floor(rx_data_,data_type = 'VTF')
                # if self.invert_Z:
                #     rx_data_ = self.apply_error_floor_Z(rx_data_)
                # if self.invert_VTF:
                #     rx_data_ = self.apply_error_floor_Tz(rx_data_)                
    
                self.data_VTF[rx] = rx_data_
                

    def filter_Z(self, max_val=1e2):
        """Filter impedance data (self.data_Z)."""
        self.data_Z = self._filter_dataframe_list(self.data_Z, max_val=max_val)
    
    def filter_VTF(self, max_val=1e2):
        """Filter tipper data (self.data_VTF)."""
        self.data_VTF = self._filter_dataframe_list(self.data_VTF, max_val=max_val)
    
    def filter_PT(self, max_val=1e2):
        """Filter phase tensor data if used."""
        self.data_PT = self._filter_dataframe_list(self.data_PT, max_val=max_val)

    def filter_all_data(self, max_val=1e2):
        """Filter all MT data: Z, VTF, PT."""
        self.data_Z   = self._filter_dataframe_list(self.data_Z,   max_val=max_val)
        self.data_VTF = self._filter_dataframe_list(self.data_VTF, max_val=max_val)
        self.data_PT  = self._filter_dataframe_list(self.data_PT,  max_val=max_val)     

    def read_observe(self, filepath):
        """
        Read Femtic observe.dat and optionally map StaIDs to station names
        using helper files (sites_imp.txt, sites_vtf.txt, sites_pt.txt).
    
        Behaviour matching EDI import:
        ----------------------------------
        • resp_Z / resp_VTF / resp_PT: StaID always a STRING
            -> station name if helper files exist
            -> else raw StaID integer converted to string
    
        • ids_Z / ids_VTF / ids_PT store the same strings.
    
        • mt_coords contains UNIQUE physical station IDs:
             - names if helper files exist
             - otherwise unified numeric IDs 1..N:
                   Z   : uid = raw_id
                   VTF : uid = raw_id - 1000
                   PT  : uid = raw_id - 2000
        """
    
        # --------------------------------------------------
        # Reset
        # --------------------------------------------------
        self.ids_Z = []
        self.ids_VTF = []
        self.ids_PT = []
    
        self.data_Z = []
        self.data_VTF = []
        self.data_PT = []
    
        self.mt_coords = None
        coord_raw = {}          # raw StaID → (north, east)
    
        # --------------------------------------------------
        # Load helper mapping files (raw_id → name)
        # --------------------------------------------------
        def load_map(fname):
            if not os.path.exists(fname):
                return None
            mapping = {}
            with open(fname) as f:
                for line in f:
                    if not line.strip():
                        continue
                    num, name = line.strip().split(maxsplit=1)
                    mapping[int(num)] = name.strip()
            return mapping
    
        map_Z   = load_map("sites_imp.txt")
        map_VTF = load_map("sites_vtf.txt")
        map_PT  = load_map("sites_pt.txt")
    
        station_map_raw = {}
        if map_Z:   station_map_raw.update(map_Z)
        if map_VTF: station_map_raw.update(map_VTF)
        if map_PT:  station_map_raw.update(map_PT)
    
        use_names = len(station_map_raw) > 0
        print("✓ Using station names" if use_names else "• No helper files → numeric IDs only")
    
        # --------------------------------------------------
        # Read observe.dat
        # --------------------------------------------------
        with open(filepath) as f:
            lines = [l.strip() for l in f]
    
        i = 0
        n = len(lines)
    
        def get_unified_id(raw):
            """Convert raw observe.dat ID → unique station ID."""
            if raw >= 2000:
                return raw - 2000
            elif raw >= 1000:
                return raw - 1000
            return raw
    
        # --------------------------------------------------
        # MAIN PARSER
        # --------------------------------------------------
        while i < n:
    
            line = lines[i]
    
            if line.startswith("END"):
                break
    
            # ===================== Z BLOCK =====================
            if line.startswith("MT"):
                _, nRx = line.split()
                nRx = int(nRx)
                self.nRx_Z = nRx
                i += 1
    
                for _ in range(nRx):
                    parts = lines[i].split()
                    raw_id = int(parts[0])
                    north  = float(parts[2])
                    east   = float(parts[3])
                    i += 1
    
                    coord_raw[raw_id] = (north, east)
    
                    sid = station_map_raw.get(raw_id, str(raw_id))
                    self.ids_Z.append(sid)
    
                    nFreq = int(lines[i]); i += 1
    
                    freq=[]; ZxxR=[]; ZxxI=[]; ZxyR=[]; ZxyI=[]
                    ZyxR=[]; ZyxI=[]; ZyyR=[]; ZyyI=[]
                    ZxxVAR=[]; ZxyVAR=[]; ZyxVAR=[]; ZyyVAR=[]
    
                    for __ in range(nFreq):
                        p = lines[i].split(); i += 1
                        freq.append(float(p[0]))
    
                        z = list(map(float, p[1:9]))
                        e = list(map(float, p[9:17]))
    
                        ZxxR.append(z[0]); ZxxI.append(z[1])
                        ZxyR.append(z[2]); ZxyI.append(z[3])
                        ZyxR.append(z[4]); ZyxI.append(z[5])
                        ZyyR.append(z[6]); ZyyI.append(z[7])
    
                        ZxxVAR.append(e[0]**2)
                        ZxyVAR.append(e[2]**2)
                        ZyxVAR.append(e[4]**2)
                        ZyyVAR.append(e[6]**2)
    
                    df = pd.DataFrame({
                        "StaID": [sid]*len(freq),
                        "Freq[Hz]": freq,
                        "Re(Zxx)_Obs": ZxxR, "Im(Zxx)_Obs": ZxxI,
                        "Re(Zxy)_Obs": ZxyR, "Im(Zxy)_Obs": ZxyI,
                        "Re(Zyx)_Obs": ZyxR, "Im(Zyx)_Obs": ZyxI,
                        "Re(Zyy)_Obs": ZyyR, "Im(Zyy)_Obs": ZyyI,
                        "Re(Zxx)_SD": np.sqrt(ZxxVAR),
                        "Re(Zxy)_SD": np.sqrt(ZxyVAR),
                        "Re(Zyx)_SD": np.sqrt(ZyxVAR),
                        "Re(Zyy)_SD": np.sqrt(ZyyVAR),
                    })
                    self.data_Z.append(df)
                continue
    
            # ===================== VTF BLOCK =====================
            if line.startswith("VTF"):
                _, nRx = line.split()
                nRx = int(nRx)
                self.nRx_VTF = nRx
                i += 1
    
                for _ in range(nRx):
                    parts = lines[i].split()
                    raw_id = int(parts[0])
                    north  = float(parts[2])
                    east   = float(parts[3])
                    i += 1
    
                    coord_raw[raw_id] = (north, east)
    
                    sid = station_map_raw.get(raw_id, str(raw_id))
                    self.ids_VTF.append(sid)
    
                    nFreq = int(lines[i]); i += 1
    
                    freq=[]; TzxR=[]; TzxI=[]; TzyR=[]; TzyI=[]
                    TzxVAR=[]; TzyVAR=[]
    
                    for __ in range(nFreq):
                        p = lines[i].split(); i += 1
                        freq.append(float(p[0]))
    
                        TzxR.append(float(p[1]));  TzxI.append(float(p[2]))
                        TzyR.append(float(p[3]));  TzyI.append(float(p[4]))
    
                        err = list(map(float, p[5:9]))
                        TzxVAR.append(err[0]**2)
                        TzyVAR.append(err[2]**2)
    
                    df = pd.DataFrame({
                        "StaID": [sid]*len(freq),
                        "Freq[Hz]": freq,
                        "Re(Tzx)_Obs": TzxR, "Im(Tzx)_Obs": TzxI,
                        "Re(Tzy)_Obs": TzyR, "Im(Tzy)_Obs": TzyI,
                        "Re(Tzx)_SD": np.sqrt(TzxVAR),
                        "Re(Tzy)_SD": np.sqrt(TzyVAR),
                    })
                    self.data_VTF.append(df)
                continue
    
            # ===================== PT BLOCK =====================
            if line.startswith("PT"):
                _, nRx = line.split()
                nRx = int(nRx)
                self.nRx_PT = nRx
                i += 1
    
                for _ in range(nRx):
                    parts = lines[i].split()
                    raw_id = int(parts[0])
                    north  = float(parts[2])
                    east   = float(parts[3])
                    i += 1
    
                    coord_raw[raw_id] = (north, east)
    
                    sid = station_map_raw.get(raw_id, str(raw_id))
                    self.ids_PT.append(sid)
    
                    nFreq = int(lines[i]); i += 1
    
                    freq=[]; TXXR=[]; TXXI=[]; TXYR=[]; TXYI=[]
                    TYXR=[]; TYXI=[]; TYYR=[]; TYYI=[]
                    TXXVAR=[]; TXYVAR=[]
    
                    for __ in range(nFreq):
                        p = lines[i].split(); i += 1
                        freq.append(float(p[0]))
    
                        vals = list(map(float, p[1:9]))
                        err  = list(map(float, p[9:13]))
    
                        TXXR.append(vals[0]); TXXI.append(vals[1])
                        TXYR.append(vals[2]); TXYI.append(vals[3])
                        TYXR.append(vals[4]); TYXI.append(vals[5])
                        TYYR.append(vals[6]); TYYI.append(vals[7])
    
                        TXXVAR.append(err[0]**2)
                        TXYVAR.append(err[2]**2)
    
                    df = pd.DataFrame({
                        "StaID": [sid]*len(freq),
                        "Freq[Hz]": freq,
                        "TXXR": TXXR, "TXXI": TXXI,
                        "TXYR": TXYR, "TXYI": TXYI,
                        "TYXR": TYXR, "TYXI": TYXI,
                        "TYYR": TYYR, "TYYI": TYYI,
                        "TXX.VAR": np.sqrt(TXXVAR),
                        "TXY.VAR": np.sqrt(TXYVAR),
                    })
                    self.data_PT.append(df)
                continue
    
            i += 1
    
        # ------------------------------------------------------
        # Build unified mt_coords
        # ------------------------------------------------------
        coord_unified = {}     # uid → dict(name, east, north, z)
    
        # build `station_map` (uid → name)
        station_map = {}
        if use_names:
            for raw_id, nm in station_map_raw.items():
                uid = get_unified_id(raw_id)
                station_map[uid] = nm
    
        for raw_id, (north, east) in coord_raw.items():
    
            uid = get_unified_id(raw_id)
    
            if use_names:
                name = station_map.get(uid, str(uid))
            else:
                name = str(uid)
    
            coord_unified[uid] = {
                "name": name,
                "east": east,
                "north": north,
                "z": 0.0
            }
    
        # convert to DataFrame
        rows = [
            (d["name"], d["east"], d["north"], d["z"])
            for uid, d in sorted(coord_unified.items())
        ]
    
        self.mt_coords = pd.DataFrame(rows, columns=["id", "east", "north", "z"])
    
        print("✓ Finished reading observe.dat")
        print(f"  Z stations:   {len(self.ids_Z)}")
        print(f"  VTF stations: {len(self.ids_VTF)}")
        print(f"  PT stations:  {len(self.ids_PT)}")




    def write_observe(self, write = True, freq_bandwidth = None, subsampling = 1, plot=True):
        
        self.prep_data(freq_bandwidth = freq_bandwidth, subsampling = subsampling)
        self.plot_inversion_periods()
        
        file_loc = os.path.join(self.outdir, 'observe.dat') 
        file = open(file_loc,'w')

        if self.invert_Z:
        
            file.write('MT    %d\n'%(self.nRx_Z))
            
            for rx in range(self.nRx_Z):
    
                ind = np.where(self.ids_Z[rx] == self.mt_coords['id'].values)[0][0]
                
                data = self.data_Z[rx]
                
                file.write('%d  %d  %2.3f  %2.3f  \n'%(rx+1, rx+1, 
                                                     self.mt_coords['north'][ind],
                                                     self.mt_coords['east'][ind]))                
                nFreq = len(data)
                file.write('%d\n'%nFreq)
                
                for freq in range(nFreq):
                    file.write('%.5f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n'%(
                               data['FREQ'][freq],
                               data['ZXXR'][freq],data['ZXXI'][freq],data['ZXYR'][freq],data['ZXYI'][freq],
                               data['ZYXR'][freq],data['ZYXI'][freq],data['ZYYR'][freq],data['ZYYI'][freq],
                               data['ZXX.VAR'][freq]**0.5,data['ZXX.VAR'][freq]**0.5,data['ZXY.VAR'][freq]**0.5,data['ZXY.VAR'][freq]**0.5,
                               data['ZYX.VAR'][freq]**0.5,data['ZYX.VAR'][freq]**0.5,data['ZYY.VAR'][freq]**0.5,data['ZYY.VAR'][freq]**0.5))
        
        
        if self.invert_VTF:
        
            file.write('VTF    %d\n'%(self.nRx_VTF))
            
            for rx in range(self.nRx_VTF):
    
                ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
                
                file.write('%d  %d  %2.3f  %2.3f  \n'%(rx+1001, rx+1001, 
                                                     self.mt_coords['north'][ind],
                                                     self.mt_coords['east'][ind]))
                data_tz = self.data_VTF[rx]
                
                # remove masked tipper points (>1+8) or inexistant points (<1e-8)
                
                for freq in range(len(data_tz)):
                    if (data_tz['TXR.EXP'][freq] > 1e+8 or
                        data_tz['TXR.EXP'][freq] > 1e+8 or 
                        data_tz['TXI.EXP'][freq] > 1e+8 or 
                        data_tz['TXI.EXP'][freq] > 1e+8 or
                        data_tz['TXR.EXP'][freq] < 1e+8 or
                        data_tz['TXR.EXP'][freq] < 1e+8 or 
                        data_tz['TXI.EXP'][freq] < 1e+8 or 
                        data_tz['TXI.EXP'][freq] < 1e+8) :
                            data_tz.drop(data_tz.index[freq])
                
                data_tz = data_tz.reset_index(drop=True)  
                
                nFreq = len(data_tz)
                file.write('%d\n'%nFreq)
                
                for freq in range(nFreq):
                    file.write('%.5f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n'%(
                               data_tz['FREQ'][freq],
                               data_tz['TXR.EXP'][freq],data_tz['TXI.EXP'][freq],data_tz['TYR.EXP'][freq],data_tz['TYI.EXP'][freq],
                               data_tz['TXVAR.EXP'][freq]**0.5,data_tz['TXVAR.EXP'][freq]**0.5,
                               data_tz['TYVAR.EXP'][freq]**0.5,data_tz['TYVAR.EXP'][freq]**0.5))
               
        
        
        file.write('END')
        file.close()
        

        
        print('Inversion data summary: \n')  
        print(f'  Inverting for {len(self.freqs_inversion)} frequencies\n\n')
        if self.invert_Z:        
            data_points_Z = np.hstack([self.data_Z[i]['FREQ'].values for i in range(self.nRx_Z)])
            unique, counts = np.unique(data_points_Z, return_counts=True)
            print('  Impedance Z: \n')
            print(f'  {self.nRx_Z} MT sites\n')
            print(f'  {len(data_points_Z)*8} Impedance tensor data points (full Z tensor, Real and Imag data) \n')
            print('  Data points per frequencies: \n')
            for i in range(len(unique)):
                #print(f'    T = {1/unique[i]} s  nData = {counts[i]} \n')
                print('    T = %.3Es  nData = %3d \n'%(1/unique[i],counts[i]))


        if self.invert_VTF:
            data_points_VTF = np.hstack([self.data_VTF[i]['FREQ'].values for i in range(self.nRx_VTF)])
            unique, counts = np.unique(data_points_VTF, return_counts=True)
            print('  Tipper VTF: \n')
            print(f'  {self.nRx_VTF} MT sites\n')
            print(f'  {len(data_points_VTF)*4} VTF full tensor (Real and Imag data) \n')
            print('  Data points per frequencies: \n')
            for i in range(len(unique)):
                #print(f'    T = {1/unique[i]} s  nData = {counts[i]} \n')
                print('    T = %.3Es  nData = %3d \n'%(1/unique[i],counts[i]))
            
                
    
            
    def plot_inversion_periods(self):
        
        fig, ax = plt.subplots(1,1, figsize=(10,3))
        
        if self.invert_Z:        
            data_points_Z = np.hstack([self.data_Z[i]['FREQ'].values for i in range(self.nRx_Z)])
            unique, counts = np.unique(data_points_Z, return_counts=True)
            ax.bar(np.log10(1/unique), counts, 0.1, color='k', alpha=0.5, label='Z') 

        if self.invert_VTF:
            data_points_VTF = np.hstack([self.data_VTF[i]['FREQ'].values for i in range(self.nRx_VTF)])
            unique, counts = np.unique(data_points_VTF, return_counts=True)
            ax.bar(np.log10(1/unique), counts, 0.1, color='r', alpha=0.5, label='VTF') 

                
        
        for i in range(self.nFreq):
            ax.axvline(x=np.log10(1/self.freqs[i]))
        for i in range(len(self.freqs_inversion)):
            ax.axvline(x=np.log10(1/self.freqs_inversion[i]), c= 'r', ls='--')   

        # ax.set_xscale('log')
        ax.set_title(f'{len(self.freqs_inversion)} inversion periods, {len(self.freqs)} original periods')
        
        ax.axvline(x=np.log10(1/self.freqs[0]), label='original periods')
        ax.axvline(x=np.log10(1/self.freqs_inversion[0]), c= 'r', ls='--', label='inversion periods')
        ax.legend()
        ax.set_xlabel('Log10 Periods (s)')
        ax.set_ylabel('nData')
            
            
    
    def write_inversion_control(self,
                                INV_METHOD = 1,
                                NUM_THREADS = 1,
                                DISTORTION = 0,
                                TRADE_OFF_PARAM = [3],
                                ITERATION = 10,
                                CONVERGE = 1.0,
                                ALPHA_WEIGHT = [1,1,1]
                                ):
        file_loc = os.path.join(self.outdir, 'control.dat') 
        file = open(file_loc,'w')
        file.write('NUM_THREADS\n')
        file.write('%d\n'%NUM_THREADS)
        file.write('MESH_TYPE\n')
        file.write('1\n')
        file.write('OUTPUT_PARAM_VTK\n')
        file.write('2\n')
        file.write('0 4\n')
        file.write('OFILE_TYPE\n')
        file.write('0\n')
        file.write('DISTORTION\n')
        file.write('%d\n'%DISTORTION)
        file.write('ALPHA_WEIGHT\n')
        for i in range(len(ALPHA_WEIGHT)):
            file.write('%.2f '%ALPHA_WEIGHT[i]) 
        file.write('\n')
        file.write('INV_METHOD\n')
        file.write('%d\n'%INV_METHOD)
        file.write('BOTTOM_RESISTIVITY\n')
        file.write('100.0\n')
        file.write('BOTTOM_ROUGHNING_FACTOR\n')
        file.write('1.0\n')
        file.write('OUTPUT_OPTION\n')
        file.write('0 0\n')
        file.write('TRADE_OFF_PARAM\n')
        for i in range(len(TRADE_OFF_PARAM)):
            file.write('%.2f '%TRADE_OFF_PARAM[i])      
        file.write('\n')    
        file.write('ITERATION\n')
        file.write('0 %d\n'%ITERATION)        
        file.write('RETRIAL\n')
        file.write('2\n') 
        file.write('CONVERGE\n')       
        file.write('%.1f\n'%CONVERGE)  
        file.write('STEP_LENGTH\n')
        file.write('0.5 0.1 1.0\n')           
        file.write('3\n')
        file.write('0.5 1.25\n')           
        file.write('END')
        file.close()    
                

    
    def plot_data_loc(self, plot_ids = False, zoom_core=False):
        fig, ax = plt.subplots(1,1,figsize=(8, 8))
        ax.plot(self.mt_coords['east'], self.mt_coords['north'], 'kv', ms=5, label = 'MT')
        if plot_ids:
            for mt_id in range(self.nRx_Z):
                ax.text(self.mt_coords['east'][mt_id], self.mt_coords['north'][mt_id]+0.1,
                         self.mt_coords['id'][mt_id])
            
        domain_boudary = np.array([[self.analysis_domain[0][0],
                          self.analysis_domain[0][0],
                          self.analysis_domain[0][1],
                          self.analysis_domain[0][1],
                          self.analysis_domain[0][0]],
                          [self.analysis_domain[1][0],
                           self.analysis_domain[1][1],
                           self.analysis_domain[1][1],
                           self.analysis_domain[1][0],
                           self.analysis_domain[1][0]]])
        
        ax.plot(domain_boudary[0], domain_boudary[1]  , 'k-', 
                 label = 'analysis_domain')
        
        if zoom_core:
            ax.set_xlim(self.mt_coords['east'].min()-10, self.mt_coords['east'].max()+10)
            ax.set_ylim(self.mt_coords['north'].min()-10, self.mt_coords['north'].max()+10)
        
        ax.legend()


    def plot_coast_line(self, ax=None):
   
        if ax is None:
            ax = plt.gca()
    
        nPolys = pd.read_csv(self.coast_line, nrows=1, names=['n'])['n'][0]
        coast_line_pd = pd.read_csv(self.coast_line, skiprows=1, sep='\s+', names=['north', 'east', 'i', 'j'])
        cb = np.where((coast_line_pd['i'] == 1) | (coast_line_pd['i'] == -1))[0]
    
        for i in range(nPolys):
            if i == 0:
                ax.plot(
                    coast_line_pd['east'][:cb[i]+1],
                    coast_line_pd['north'][:cb[i]+1],
                    'k-', lw=1
                )
            else:
                ax.plot(
                    coast_line_pd['east'][cb[i-1]+1:cb[i]+1],
                    coast_line_pd['north'][cb[i-1]+1:cb[i]+1],
                    'k-', lw=1
                )


    # def plot_coast_line(self):
    
    #     nPolys = pd.read_csv(self.coast_line, nrows=1, names=['n'])['n'][0]
    #     coast_line_pd = pd.read_csv(self.coast_line, skiprows=1, sep='\s+', names = ['north','east','i','j'])
    #     cb = np.where((coast_line_pd['i']==1) | (coast_line_pd['i']==-1))[0]
        
    #     for i in range(nPolys):
    #         if i == 0:
    #             plt.plot(coast_line_pd['east'][:cb[i]+1],coast_line_pd['north'][:cb[i]+1], 'k-',lw=1)
    #         else:
    #             plt.plot(coast_line_pd['east'][cb[i-1]+1:cb[i]+1],coast_line_pd['north'][cb[i-1]+1:cb[i]+1], 'k-',lw=1)
    
    
    def plot_topo_bathy(self, ax=None, vmin=None, vmax=None):

        if ax is None:
            ax = plt.gca()
    
        # Topography
        topo_pd = pd.read_csv(self.topography, sep='\s+', names=['north', 'east', 'z'])
        topo_pd = topo_pd[topo_pd['z'] > 0]
        ax.scatter(
            topo_pd['east'], topo_pd['north'], c=topo_pd['z'],
            cmap='gray', s=2, vmin=vmin, vmax=vmax
        )
    
        # Bathymetry
        bathy_pd = pd.read_csv(self.bathymetry, sep='\s+', names=['north', 'east', 'z'])
        bathy_pd = bathy_pd[bathy_pd['z'] > 0]
        ax.scatter(
            bathy_pd['east'], bathy_pd['north'], c=bathy_pd['z'],
            cmap='Blues', s=2, vmin=vmin, vmax=vmax
        )
    

    # def plot_topo_bathy(self, **kwargs):
    #     vmin = kwargs.get('vmin', None)
    #     vmax = kwargs.get('vmax', None)
    
    #     topo_pd = pd.read_csv(self.topography, sep='\s+', names=['north', 'east', 'z'])
    #     topo_pd = topo_pd[topo_pd['z'] > 0]
    #     plt.scatter(
    #         topo_pd['east'], topo_pd['north'], c=topo_pd['z'],
    #         cmap='gray', s=2, vmin=vmin, vmax=vmax
    #     )
    
    #     bathy_pd = pd.read_csv(self.bathymetry, sep='\s+', names=['north', 'east', 'z'])
    #     bathy_pd = bathy_pd[bathy_pd['z'] > 0]
    #     plt.scatter(
    #         bathy_pd['east'], bathy_pd['north'], c=bathy_pd['z'],
    #         cmap='Blues', s=2, vmin=vmin, vmax=vmax
    #     )           
        

    #def plot_input_data(self):
        



class MeshGen():
    """
    MeshGen prepares the input data required to run the sequential shell 
    script that creates the input model and .vtk visualization files
    
    Files required: 
        - data coordinates
        - coast_line 
        - topography 
        - bathymetry 
    
    Writes out:
        - analysis_domain.dat
        - control.dat
        - makeMtr.param
        - obs_site.dat
        - observing_site.dat
        - resistivity_attr.dat

    """
    
    def __init__(self, survey, analysis_domain, center, mt_coords, nRx, outdir):
        self.survey = survey
        
        self.outdir = outdir
                
        ## Domains to include
        self.land = True
        self.sea = False
        #self.coast_line = False
        
        
        ## Size of computational domain (in km)
        self.analysis_domain = analysis_domain
        
        ## Mesh 
        self.center = center  #center (in km)
        self.rotation = 0.0
        
        # MT sites data
        self.mt_coords = mt_coords
        self.median_elevation = - np.median(self.mt_coords['z'])
        self.center = [self.center[0], self.center[1], self.median_elevation]
        # self.mt_data = mt_data
        self.nRx = nRx

        
        # Information about the ellipsoids to  control edge lengths
        self.ellipsoids_control = [6,
                            [40.0,  1.0, 0.0, 0.5, 0.7],
                            [60.0,  5.0, 0.0, 0.3, 0.5],
                            [100.0, 10.0, 0.0, 0.1, 0.3],
                            [200.0, 20.0, 0.0, 0.0, 0.0],
                            [300.0, 30.0, 0.0, 0.0, 0.0],
                            [500.0, 50.0, 0.0, 0.0, 0.0]]

        
        # topography interpolation parameters
        self.interpolate_sr = 100
        self.interpolate_npts = 3
        self.interpolate_minval = 1e-6
        # altitude parameters:
        self.altitude_file = 'topography.dat'
        self.altitude_min = 0.0
        self.altitude_max = 1e20
        # sea_depth parameters:
        self.sea_depth_file = 'bathymetry.dat'
        self.sea_depth_min = 0.01
        self.sea_depth_max = 1e20     
        
        self.ellipsoids_mtr = [10,
                            [40.0,    1.0,  0.0,  0.7,  0.9],
                            [45.0,   1.5,  0.0,  0.5,  0.7],
                            [50.0,    3.0,  0.0,  0.4,  0.7],
                            [60.0,    5.0,  0.0,  0.3,  0.5],
                            [80.0,    8.0,  0.0,  0.1,  0.3],
                            [100.0,  10.0,  0.0,  0.0,  0.0],
                            [200.0,  20.0,  0.0,  0.0,  0.0],
                            [300.0,  30.0,  0.0,  0.0,  0.0],
                            [400.0,  40.0,  0.0,  0.0,  0.0],
                            [500.0,  50.0,  0.0,  0.0,  0.0]]
        
        self.ellipsoids_observing_sites =  [5, 
                                            0.1, 0.02,  
                                            0.3, 0.05,  
                                            1.0, 0.10,  
                                            3.0, 0.30,  
                                            5.0, 0.50]        

        self.ellipsoids_obs_sites = [6,
                            [0.5, 0.10, 0.3],
                            [1.0, 0.20, 0.3],
                            [1.5, 0.30, 0.3],
                            [2.0, 0.50, 0.3],
                            [3.0, 1.00, 0.3],
                            [5.0, 2.00, 0.3]]

        
        self.ellipsoids_resistivity_attr = [9,
                            [40.0,       2.0,  0.0,  0.7],
                            [45.0  ,     3.0,  0.0,  0.7],
                            [50.0,       5.0,  0.0,  0.7],
                            [60.0,      10.0,  0.0,  0.6],
                            [100.0,    100.0,  0.0,  0.5],
                            [200.0,    200.0,  0.0,  0.3],
                            [300.0,    300.0,  0.0,  0.2],
                            [500.0,    500.0,  0.0,  0.1],
                            [1000.0,  1000.0,  0.0,  0.0]] 
        
        
        self.ellipsoids_resistivity_attr_sites = [2,
                            [3.0, 2.0],
                            [5.0, 3.0]]
        
        self.resistivity_sea = 0.25
        self.resistivity_starting_model = 100

    def write_analysis_domain(self):
        file_loc = os.path.join(self.outdir, 'analysis_domain.dat') 
        file = open(file_loc,'w') 
        for i in range(3):
            file.write(' %4.1f %4.1f\n'%(self.analysis_domain[i][0], self.analysis_domain[i][1]))
        file.close() 
        
    
    def write_control(self):
        file_loc = os.path.join(self.outdir, 'control.dat') 
        file = open(file_loc,'w')      
        file.write('CENTER\n')
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('ROTATION\n')
        file.write('%4.1f\n'%(self.rotation))
        file.write('NUM_THREADS\n1\nSURF_MESH\nELLIPSOIDS\n')
        file.write('%d'%(self.ellipsoids_control[0]))
        for tt in self.ellipsoids_control[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.write('\n')
        file.write('INTERPOLATE\n')
        file.write('%4.1f\n'%self.interpolate_sr)
        file.write('%4.1f\n'%self.interpolate_npts)
        file.write('%4.7f\n'%self.interpolate_minval)

        file.write('ALTITUDE\n')
        file.write('topography.dat\n')
        file.write('%4.1f\n'%self.altitude_min)
        file.write('%4.1f\n'%self.altitude_max)          

        file.write('SEA_DEPTH\n')
        file.write('bathymetry.dat\n')
        file.write('%4.1f\n'%self.sea_depth_min)
        file.write('%4.1f\n'%self.sea_depth_max) 
        file.write('END')
        file.close() 


    def write_makeMtr(self):
        file_loc = os.path.join(self.outdir, 'makeMtr.param') 
        file = open(file_loc,'w')      
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('%4.1f\n'%(self.rotation))
        file.write('%d'%(self.ellipsoids_mtr[0]))
        for tt in self.ellipsoids_mtr[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.close() 
        

    def write_obs_site(self):
        file_loc = os.path.join(self.outdir, 'obs_site.dat') 
        file = open(file_loc,'w')   
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f %2.3f %2.3f\n'%(self.mt_coords['north'][rx],
                                              self.mt_coords['east'][rx],
                                              self.mt_coords['z'][rx]))
            file.write('%d'%(self.ellipsoids_obs_sites[0]))
            for tt in self.ellipsoids_obs_sites[1:]:
                file.write('\n')
                for ttt in tt:
                    file.write('%4.3f '%(ttt))
            file.write('\n')
        file.close() 



    def write_observing_site(self):
        file_loc = os.path.join(self.outdir, 'observing_site.dat') 
        file = open(file_loc,'w')    
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f  %2.3f  '%(self.mt_coords['north'][rx],
                                         self.mt_coords['east'][rx]))
            file.write('%d  '%(self.ellipsoids_observing_sites[0]))
            for tt in self.ellipsoids_observing_sites[1:]:
                file.write('%4.3f '%(tt))
            file.write('\n')
        file.close() 
        
        
    def write_dummy_coast_line(self):
        file_loc = os.path.join(self.outdir, 'coast_line.dat') 
        file = open(file_loc,'w')   
        # define a boundary which covers whole of the computational domain.
        domain_boudary = np.array([[self.analysis_domain[0][0]-1,
                          self.analysis_domain[0][0]-1,
                          self.analysis_domain[0][1]+1,
                          self.analysis_domain[0][1]+1,
                          self.analysis_domain[0][0]-1],
                          [self.analysis_domain[1][0]-1,
                            self.analysis_domain[1][1]+1,
                            self.analysis_domain[1][1]+1,
                            self.analysis_domain[1][0]-1,
                            self.analysis_domain[1][0]-1]])
        
        file.write('1\n')
        for corner in range(3):    
            file.write('%4.3f %4.3f 0 0\n'%(domain_boudary[0][corner],
                                            domain_boudary[1][corner]))
        file.write('%4.3f %4.3f 1 0\n'%(domain_boudary[0][3],
                                        domain_boudary[1][3]))
        file.close() 

        
        
    def write_dummy_topo_file(self):
        file_loc = os.path.join(self.outdir, 'topography.dat') 
        file = open(file_loc,'w') 
        for north in np.arange(self.analysis_domain[0][0]-11, self.analysis_domain[0][1]+11,10):
            for east in np.arange(self.analysis_domain[1][0]-11, self.analysis_domain[1][1]+11,10):
                file.write('%4.3f %4.3f %4.3f\n'%(north, east, 0.000))
        file.close() 



    def write_dummy_bathy_file(self):
        file_loc = os.path.join(self.outdir, 'bathymetry.dat') 
        file = open(file_loc,'w') 
        for north in np.arange(self.analysis_domain[0][0]-11, self.analysis_domain[0][1]+11,10):
            for east in np.arange(self.analysis_domain[1][0]-11, self.analysis_domain[1][1]+11,10):
                file.write('%4.3f %4.3f %4.3f\n'%(north, east, -0.01))
        file.close() 

        

    def write_resistivity_attr(self):
        file_loc = os.path.join(self.outdir, 'resistivity_attr.dat') 
        file = open(file_loc,'w') 
        self.region_attributes = np.sum([self.land, self.sea]) + 1
        file.write('%d\n'%self.region_attributes)
        file.write('10 1.0e+9 -1 1\n')
        if self.sea:
            file.write('20 %.3f  -1 1\n'%self.resistivity_sea)
            file.write('30 %.1f   9 0\n'%self.resistivity_starting_model)
        else:
            file.write('20 %.1f   9 0\n'%self.resistivity_starting_model)
            
        file.write('%4.1f %4.1f %4.1f\n'%(self.center[0],self.center[1],self.center[2]))
        file.write('%4.1f\n'%(self.rotation))
        
        file.write('%d'%(self.ellipsoids_resistivity_attr[0]))
        for tt in self.ellipsoids_resistivity_attr[1:]:
            file.write('\n')
            for ttt in tt:
                file.write('%4.1f '%(ttt))
        file.write('\n')
        file.write('%d\n'%self.nRx)
        for rx in range(self.nRx):
            file.write('%2.3f  %2.3f  %2.3f\n'%(self.mt_coords['north'][rx],
                                         self.mt_coords['east'][rx],
                                         self.mt_coords['z'][rx]))
            file.write('%d  '%(self.ellipsoids_resistivity_attr_sites[0]))
            for tt in self.ellipsoids_resistivity_attr_sites[1:]:
                file.write('\n')
                for ttt in tt:
                    file.write('%4.3f '%(ttt))
            file.write('\n')
        file.close() 
        
    
    def write_inputs(self):
        self.write_analysis_domain()
        self.write_makeMtr()
        self.write_obs_site()
        self.write_observing_site()
        self.write_resistivity_attr()
        self.write_control()
        


     
class InvResults():
    """
    InvResults plot the statistics of the inversion
    No Files required
    """
    
    # def __init__(self, survey, mt_coords, mt_data, mt_ids, results_directory):
    def __init__(self, survey, mt_coords, ids_Z, ids_VTF, results_directory):
        self.mt_coords = mt_coords
        self.ids_Z = ids_Z
        self.ids_VTF = ids_VTF
        #self.nRx = len(mt_data)
        self.dir = results_directory 
        self.plot_dir = None
 
        self.nIter = None
        
        self.rms_total = None
        self.rms_Z = None
        self.rms_VTF = None
        self.rms_breakDown_Z = None
        self.rms_breakDown_VTF = None
        self.log_inversion = None
        
        ## Data to be inverted for
        self.invert_Z = False
        self.invert_VTF = False
        self.invert_PT = False
        
        #self.obs = mt_data
        self.resp_Z = None  
        self.resp_VTF = None 
        self.resp_PT = None 
    
    

    def _get_background_image(self, inversion, extent, topo_vmin, topo_vmax, filename="background.png"):
    
        # If image already exists, return its path
        if os.path.exists(filename):
            print(f"Using cached background image: {filename}")
            return filename
    
        # Otherwise, generate and save the background
        print("Generating and saving background image...")
        fig, ax = plt.subplots()
        inversion.plot_topo_bathy(ax=ax, vmin=topo_vmin, vmax=topo_vmax)
        inversion.plot_coast_line(ax=ax)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect('equal')
        ax.axis('off')
        plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        return filename
    
    
    def create_res_dir(self):
        if self.plot_dir is None:
            self.plot_dir = './%s/plots/it%02d/'%(self.dir,self.nIter)
            if not os.path.isdir('./%s/plots/'%self.dir):
                os.mkdir('./%s/plots'%self.dir)
            if not os.path.isdir(self.plot_dir):
                os.mkdir(self.plot_dir)  
                

    def open_single_csv(self, csv_path):

        nMT  = 0
        nVTF = 0
        nPT  = 0
    
        with open(csv_path, 'r') as file:
            
            lines = file.readlines()
            
            for i, line in enumerate(lines):
                line = line.strip()
                if 'MT' in line:
                    data_type = 'MT'
                    line_MT = i
                    self.invert_Z = True            
                elif 'VTF' in line:
                    data_type = 'VTF'
                    line_VTF = i
                    self.invert_VTF = True
                elif 'PT' in line:
                    data_type = 'PT'
                    line_PT = i
                    self.invert_PT = True
                elif 'StaID' in line:
                    pass
                elif data_type:  
                    if data_type == 'MT':
                        nMT += 1
                    elif data_type == 'VTF':
                        nVTF += 1
                    elif data_type == 'PT':
                        nPT += 1
    
        #nDatasets = invert_Z + invert_VTF + invert_PT
        
        class output():
            pass
    
        if self.invert_Z:
            headers = [header.strip() for header in lines[line_MT+1].split(",")]
            data = [line.split(",") for line in lines[line_MT+2:line_MT+2+nMT]]
            df_Z = pd.DataFrame(data, columns=headers)
            df_Z.pop(df_Z.columns[-1])
            # cols = df_Z.columns.difference(['StaID'])
            # df_Z[cols] = df_Z[cols].astype(float).astype(int)
            # df_Z['StaID'] = df_Z['StaID'].apply(lambda x: x.strip())
            df_Z['StaID'] = df_Z['StaID'].astype(int) 
            df_Z = df_Z.astype({col: float for col in df_Z.columns[1:]})
            output.Z = df_Z
    
        if self.invert_VTF:
            headers = [header.strip() for header in lines[line_VTF+1].split(",")]
            data = [line.split(",") for line in lines[line_VTF+2:line_VTF+2+nVTF]]
            df_VTF = pd.DataFrame(data, columns=headers)
            df_VTF.pop(df_VTF.columns[-1])
            # cols = df_VTF.columns.difference(['StaID'])
            # df_VTF[cols] = df_VTF[cols].astype(float).astype(int)
            # df_VTF['StaID'] = df_VTF['StaID'].apply(lambda x: x.strip())
            df_VTF['StaID'] = df_VTF['StaID'].astype(int) 
            df_VTF = df_VTF.astype({col: float for col in df_VTF.columns[1:]})
            output.VTF = df_VTF
    
        if self.invert_PT:
            headers = [header.strip() for header in lines[line_VTF+1].split(",")]
            data = [line.split(",") for line in lines[line_PT+2:line_PT+2+nPT]]
            df_PT = pd.DataFrame(data, columns=headers)
            df_PT.pop(df_PT.columns[-1])
            # cols = df_PT.columns.difference(['StaID'])
            # df_PT[cols] = df_PT[cols].astype(float).astype(str)
            # df_PT['StaID'] = df_PT['StaID'].apply(lambda x: x.strip())
            df_PT['StaID'] = df_PT['StaID'].astype(int) 
            df_PT.astype({col: float for col in df_PT.columns[1:]})
            output.PT = df_PT
    
        return output
    
    
    def read_result_csv(self):
        
        self.create_res_dir()
    
        for root, dirs, files in os.walk(self.dir):
            for file in (files):
                if file.endswith('%d.csv'%self.nIter):
                    csv_path = '%s/%s'%(self.dir, file)
                    data_csv = self.open_single_csv(csv_path)
                    
                    if self.invert_Z:
                        try:
                            resp_Z_all
                        except NameError:
                            resp_Z_all = data_csv.Z
                        else:
                            resp_Z_all = pd.concat([resp_Z_all, data_csv.Z], ignore_index=True)
        
                            
                    if self.invert_VTF:
                        try:
                            resp_VTF_all
                        except NameError:
                            resp_VTF_all = data_csv.VTF
                        else:
                            resp_VTF_all = pd.concat([resp_VTF_all, data_csv.VTF], ignore_index=True)
                            
                            
                    if self.invert_PT:
                        try:
                            resp_PT_all
                        except NameError:
                            resp_PT_all = data_csv.PT
                        else:
                            resp_PT_all = pd.concat([resp_PT_all, data_csv.PT], ignore_index=True)
        
        
        if self.invert_Z:
            resp_Z_all.columns = resp_Z_all.columns.str.lstrip()
            resp_Z_all = resp_Z_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_Z_all = resp_Z_all.reset_index(drop=True)   
        
            # for site in range(len(self.ids_Z)):
            #     idx = np.where(resp_Z_all['StaID'] == site+1) 
            #     # print(idx[0], '----', self.mt_ids[site])
            #     resp_Z_all['StaID'][idx[0]] = self.ids_Z[site]
            
            #this attemps to fix the original code above which fails when station ID numbers aren''t sequential (if all stations don't have VTF or PT)
            unique_ids = resp_Z_all['StaID'].unique()
            if len(unique_ids) != len(self.ids_Z):
                print("Warning: Number of unique station IDs in CSV does not match self.ids_Z.")
                print("unique StaIDs found:", unique_ids)
                print("expected StaIDs:", self.ids_Z)

            for i, sta_id in enumerate(sorted(unique_ids)):
                if i < len(self.ids_Z):
                    idx = resp_Z_all['StaID'] == sta_id
                    resp_Z_all.loc[idx, 'StaID'] = self.ids_Z[i]
        
            self.resp_Z = resp_Z_all        
    

        if self.invert_VTF:
            resp_VTF_all.columns = resp_VTF_all.columns.str.lstrip()
            resp_VTF_all = resp_VTF_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_VTF_all = resp_VTF_all.reset_index(drop=True)   
        
            # for site in range(len(self.ids_VTF)):
            #     idx = np.where(resp_VTF_all['StaID'] == site+1001) 
            #     resp_VTF_all['StaID'][idx[0]] = self.ids_VTF[site]
            
            #this attemps to fix the original code above which fails when station ID numbers aren''t sequential (if all stations don't have VTF or PT)
            unique_ids = resp_VTF_all['StaID'].unique()
            if len(unique_ids) != len(self.ids_VTF):
                print("Warning: Number of unique station IDs in CSV does not match self.ids_VTF.")
                print("unique StaIDs found:", unique_ids)
                print("expected StaIDs:", self.ids_VTF)

            for i, sta_id in enumerate(sorted(unique_ids)):
                if i < len(self.ids_VTF):
                    idx = resp_VTF_all['StaID'] == sta_id
                    resp_VTF_all.loc[idx, 'StaID'] = self.ids_VTF[i]
            
            
            resp_VTF_all['Im(Tzx)_Cal'] = -resp_VTF_all['Im(Tzx)_Cal']
            resp_VTF_all['Im(Tzy)_Cal'] = -resp_VTF_all['Im(Tzy)_Cal']
            resp_VTF_all['Im(Tzx)_Obs'] = -resp_VTF_all['Im(Tzx)_Obs']
            resp_VTF_all['Im(Tzy)_Obs'] = -resp_VTF_all['Im(Tzy)_Obs']
        
            self.resp_VTF = resp_VTF_all    




        if self.invert_PT:
            resp_PT_all.columns = resp_PT_all.columns.str.lstrip()
            resp_PT_all = resp_PT_all.sort_values(['StaID', 'Freq[Hz]'],
                          ascending = [True, False])     
        
            resp_PT_all = resp_PT_all.reset_index(drop=True)   
        
            # for site in range(len(self.mt_ids)):
            #     idx = np.where(resp_PT_all['StaID'] == site+1) 
            #     resp_PT_all['StaID'][idx[0]] = self.mt_ids[site]

            #this attemps to fix the original code above which fails when station ID numbers aren''t sequential (if all stations don't have VTF or PT)
            unique_ids = resp_PT_all['StaID'].unique()
            if len(unique_ids) != len(self.ids_Z):
                print("Warning: Number of unique station IDs in CSV does not match self.ids_Z.")
                print("unique StaIDs found:", unique_ids)
                print("expected StaIDs:", self.ids_PT)

            for i, sta_id in enumerate(sorted(unique_ids)):
                if i < len(self.ids_PT):
                    idx = resp_PT_all['StaID'] == sta_id
                    resp_PT_all.loc[idx, 'StaID'] = self.ids_PT[i]

        
            self.resp_PT = resp_PT_all
    
    



    def plot_cnv(self, save_plot=True):
        # log_inv = pd.read_csv('%s/femtic.cnv'%self.dir)
        log_inv = pd.read_csv('%s/femtic.cnv'%self.dir, sep='\s+')
        
        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('Iter#')
        ax1.set_ylabel('RMS', color=color)
        ax1.plot(log_inv['Iter#'], log_inv['RMS'], 'o-',color=color ,ms=3)
        ax1.tick_params(axis='y', labelcolor=color)        
        
        ax1.set_xlim(0,log_inv['Iter#'].max() +1)
        ax1.set_yticks(np.arange(0,np.ceil(log_inv['RMS'].max()) +1 ))
        # ax.set_yticks(np.arange(ylim[0], ylim[1]))
        ax1.grid(lw=0.3)

        ax2 = ax1.twinx() 

        color = 'tab:blue'
        ax2.set_ylabel('Roughness', color=color)  
        ax2.plot(log_inv['Iter#'], log_inv['Roughness'], 'o-',color=color ,ms=3)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout() 
        # plt.show
        
        if save_plot:
            plt.savefig('%s/log_rms.pdf'%(self.plot_dir),dpi=300,format='pdf', bbox_inches='tight')
            # plt.close('all')


    
        
        
    def z2rhophy(self, FREQ,ZR,ZI,dZ):
        
        ZR *= 10000/(4*np.pi)
        ZI *= - 10000/(4*np.pi)
        dZ *= 10000/(4*np.pi)
        
        # # calcul of apparent resistivity and phases
        rho = ((ZR**2+ZI**2)*0.2/(FREQ))
        phy = np.degrees(np.arctan2(ZI,ZR))
        # # calcul of errors
        drho = 2*rho*dZ / (((ZR**2+ZI**2)**0.5))
        dphy = np.degrees(0.5 * (drho/rho))
        log10_drho = (0.3772 * (dZ**2/(ZR**2+ZI**2)))**0.5

        return (rho, phy, drho, dphy,log10_drho)  
    
    
    
    def compute_rms(self, residuals):
        lists_res = [np.array(residual.values.tolist(), dtype=float).flatten() for residual in residuals]
        combined = np.array([item for sublist in lists_res for item in sublist], dtype=float)
        rms = (sum(combined**2) / (len(combined)))**0.5
        return rms




    def compute_rms_breakDown(self):
        
        # total RMS
        residuals = []
        if self.invert_Z:
            filter_col = [col for col in self.resp_Z if col.endswith(('Res'))]
            residuals_Z = self.resp_Z[filter_col]
            self.rms_Z  = self.compute_rms([residuals_Z])
            residuals.append(residuals_Z)
        if self.invert_VTF:
            filter_col = [col for col in self.resp_VTF if col.endswith(('Res'))]
            residuals_VTF = self.resp_VTF[filter_col]
            self.rms_VTF  = self.compute_rms([residuals_VTF])
            residuals.append(residuals_VTF)            
        if self.invert_PT:
           filter_col = [col for col in self.resp_PT if col.endswith(('Res'))]
           residuals_PT = self.resp_PT[filter_col]
           residuals.append(residuals_PT)

        self.rms_total  = self.compute_rms(residuals)
        
        self.rms_total_by_site = {}

        if self.invert_Z or self.invert_VTF or self.invert_PT:
            all_residuals = []
        
            if self.invert_Z:
                all_residuals.append(self.resp_Z)
        
            if self.invert_VTF:
                all_residuals.append(self.resp_VTF)
            
            if self.invert_PT:
                all_residuals.append(self.resp_PT)
        
            combined_df = pd.concat(all_residuals, ignore_index=True)
            
            for site in self.ids_Z:
                z_id   = site
                vtf_id = str(int(site) + 1000)
            
                df_site = combined_df[(combined_df['StaID'] == z_id) |
                                      (combined_df['StaID'] == vtf_id)]
            
                filter_cols = [col for col in df_site.columns if col.endswith('Res')]
            
                if not df_site.empty and filter_cols:
                    residuals_site = df_site[filter_cols].values.flatten()
                    residuals_site = residuals_site[~np.isnan(residuals_site)]
            
                    rms_site = np.sqrt(np.mean(residuals_site ** 2))
                    self.rms_total_by_site[site] = rms_site
                else:
                    self.rms_total_by_site[site] = np.nan
            
            # for site in self.ids_Z:
            #     df_site = combined_df[combined_df['StaID'] == site]
            #     filter_cols = [col for col in df_site.columns if col.endswith('Res')]
            #     if not df_site.empty and filter_cols:
            #         residuals_site = df_site[filter_cols].values.flatten()
            #         residuals_site = residuals_site[~np.isnan(residuals_site)]
            #         rms_site = np.sqrt(np.mean(residuals_site ** 2))
            #         self.rms_total_by_site[site] = rms_site
        
        if self.invert_Z:
        
            self.rms_breakDown_Z = pd.DataFrame(columns=['StaID','Total','Zxx','Zxy','Zyx','Zyy'])
        
            # RMS / site
            rms_sites = []
            rms_sites_id = []
            filter_col = [col for col in self.resp_Z if col.endswith(('StaID','Res'))]
            residuals = self.resp_Z[filter_col]   
            
            for site in range(len(self.ids_Z)):
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_Z[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                rms_site = self.compute_rms([residuals_site])
                rms_sites_id.append(self.ids_Z[site])
                rms_sites.append(rms_site)
            rms_sites = np.array(rms_sites)
            
            
            # RMS / site / component
            rms_sites_comps = []
            for site in range(len(self.ids_Z)):
                rms_comps = []
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_Z[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                comp = ['Zxx','Zxy','Zyx','Zyy']
                for i in range(4):
                    filter_col = [col for col in residuals_site if comp[i] in col]
                    residuals_site_comp = residuals_site[filter_col]
                    rms_comp = self.compute_rms([residuals_site_comp])
                    rms_comps.append(rms_comp)
                rms_sites_comps.append(rms_comps)
            rms_sites_comps  = np.array(rms_sites_comps)
    
            self.rms_breakDown_Z['StaID'] = rms_sites_id
            self.rms_breakDown_Z['Total'] = rms_sites
            self.rms_breakDown_Z[['Zxx','Zxy','Zyx','Zyy']] = rms_sites_comps


        if self.invert_VTF:
        
            self.rms_breakDown_VTF = pd.DataFrame(columns=['StaID','Total','Tzx','Tzy'])
        
            # RMS / site
            rms_sites = []
            rms_sites_id = []
            filter_col = [col for col in self.resp_VTF if col.endswith(('StaID','Res'))]
            residuals = self.resp_VTF[filter_col]   
            
            for site in range(len(self.ids_VTF)):
                print(f"Processing site {site} / {len(self.ids_VTF)}")
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_VTF[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                print(f"Residuals shape: {residuals_site.shape}")
                rms_site = self.compute_rms([residuals_site])
                rms_sites_id.append(self.ids_VTF[site])
                rms_sites.append(rms_site)
            rms_sites = np.array(rms_sites)
            
            
            # RMS / site / component
            rms_sites_comps = []
            for site in range(len(self.ids_VTF)):
                rms_comps = []
                residuals_site = residuals.loc[np.where(residuals['StaID'] == self.ids_VTF[site])[0]]
                residuals_site = residuals_site.drop(['StaID'], axis=1)
                comp = ['Tzx','Tzy']
                for i in range(2):
                    filter_col = [col for col in residuals_site if comp[i] in col]
                    residuals_site_comp = residuals_site[filter_col]
                    rms_comp = self.compute_rms([residuals_site_comp])
                    rms_comps.append(rms_comp)
                rms_sites_comps.append(rms_comps)
            rms_sites_comps  = np.array(rms_sites_comps)
    
            self.rms_breakDown_VTF['StaID'] = rms_sites_id
            self.rms_breakDown_VTF['Total'] = rms_sites
            self.rms_breakDown_VTF[['Tzx','Tzy']] = rms_sites_comps



    def read_distorsion(self):
        for root, dirs, files in os.walk(self.dir):
            for file in (files):
                if file.endswith('distortion_iter%d.dat'%self.nIter):
                    file_path = '%s/%s'%(self.dir, file)

        self.gd_data = pd.read_csv(file_path, skiprows=1 ,sep='\s+',names=['Cxx','Cxy','Cyx','Cyy','act'])


    def calc_gd_strength(self):
        # following Adveeda et al., 2015
        self.gd_data['gd'] = 0
        for i in range(len(self.gd_data)):
            gd_matrix = np.array([[self.gd_data['Cxx'].iloc[i], self.gd_data['Cxy'].iloc[i]], [self.gd_data['Cyx'].iloc[i], self.gd_data['Cyy'].iloc[i]]])
            self.gd_data['gd'].iloc[i] = np.linalg.norm(gd_matrix - np.identity(2))
    
    def plot_distorsion_map(self,
                            xlim=[-100, 100], ylim=[-100, 100], 
                            vmin=0, vmax=3,
                            save_plot=False,
                            label_sites=False,
                            topo_vmin=0,
                            topo_vmax=1.2,
                            filename="gd_map.pdf",
                            inversion=None):
        """
        Plots Galvanic Distortion strength as a colored scatter map, with optional site labels.
        """
    
        self.read_distorsion()
        self.calc_gd_strength()
    
        extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
    
        fig, ax = plt.subplots(figsize=(6, 6))
    
        # Background from cached image
        background_path = None
        if inversion is not None:
            background_path = self._get_background_image(
                inversion,
                extent=extent,
                topo_vmin=topo_vmin,
                topo_vmax=topo_vmax,
                filename="./background.png"
            )
            if background_path:
                img = plt.imread(background_path)
                ax.imshow(img, extent=extent, origin='upper', aspect='equal', zorder=0)
    
        print("Plotting Galvanic Distortion strength map...")
    
        # Color-coded scatter
        sc = ax.scatter(self.mt_coords['east'], self.mt_coords['north'],
                        c=self.gd_data['gd'],
                        cmap='jet',
                        vmin=vmin,
                        vmax=vmax,
                        s=40,
                        edgecolors='k',
                        linewidths=0.5)
    
        # Add optional labels
        if label_sites:
            for i, row in self.mt_coords.iterrows():
                ax.annotate(row['id'], (row['east'], row['north']),
                            textcoords="offset points", xytext=(3, 3),
                            ha='left', fontsize=6, color='black')
    
        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
        cbar.set_label("Galvanic Distortion Strength")
    
        # Formatting
        ax.set_title("Galvanic Distortion Strength")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')
        ax.grid(lw=0.1)
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")
    
        if save_plot:
            os.makedirs(self.plot_dir, exist_ok=True)
            path = os.path.join(self.plot_dir, filename)
            plt.savefig(path, dpi=300, format='pdf', bbox_inches="tight")
            print(f"Saved: {path}")
            plt.close(fig)
        else:
            plt.show()



    # def plot_distorsion_map(self,
    #                  ax, 
    #                  xlim = [-100, 100],ylim = [-100, 100], 
    #                  vmin=1,vmax=3,
    #                  save_plot = False):

    #     self.read_distorsion()
    #     self.calc_gd_strength()

    #     print('Plotting nRMSE map for the impedance Z...')
    #     pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c='w',
    #             #vmin=vmin,vmax=vmax,cmap=cmap,
    #             marker='o', s=1+self.gd_data['gd'], linewidths=0.5,edgecolors='k',)
    #     ax.scatter(10,-10,s=1, linewidths=0.5,edgecolors='k',)
    #     ax.scatter(10,-10,s=10, linewidths=0.5,edgecolors='k',)
    #     ax.set_title('Galvanic Distortion strength')
    
    #     #ax.axis('equal')
    #     ax.set_xlim(xlim)
    #     ax.set_ylim(ylim)
    #     ax.grid(lw=0.1)
    #     ax.set_xlabel('East (km)')
    #     ax.set_ylabel('North (km)')
    #     #plt.tight_layout()
        
    #     if save_plot:
    #             plt.savefig('%s/gd_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')



    
    def plot_rms_map(self,
                     Z=True,
                     VTF=False,
                     PT=False,
                     xlim=[-100, 100], ylim=[-100, 100], 
                     vmin=1, vmax=3,
                     ms_size=40,
                     save_plot=True,
                     filename=None,
                     topo_vmin=0,
                     topo_vmax=1.2,
                     inversion=None):
           
        cmap = cm.get_cmap('jet', 12)
        plot_types = []
    
        if Z: plot_types.append("Z")
        if VTF: plot_types.append("VTF")
        if PT: plot_types.append("PT")
    
        if not plot_types:
            print("Nothing selected to plot.")
            return
    
        n_plots = len(plot_types)
        fig, axs = plt.subplots(1, n_plots, figsize=(6 * n_plots, 6))
        axs = np.atleast_1d(axs)
    
        # === Prepare reusable background (topo + coast) ===
        background_artists = []
        # Define extent once
        extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
        
        # Generate or retrieve background image file
        background_path = None
        if inversion is not None:
            background_path = self._get_background_image(
                inversion,
                extent=extent,
                topo_vmin=topo_vmin,
                topo_vmax=topo_vmax,
                filename="./background.png"
                )
    
        for i, plot_type in enumerate(plot_types):
            ax = axs[i]
        
            # Use imshow to show the cached background image
            if background_path:
                img = plt.imread(background_path)
                ax.imshow(img, extent=extent, origin='upper', aspect='equal', zorder=0)
    
            # === Plot RMS ===
            if plot_type == "Z":
                print("Plotting nRMSE for Z...")
                pc = ax.scatter(
                    self.mt_coords['east'], self.mt_coords['north'],
                    c=self.rms_breakDown_Z['Total'],
                    vmin=vmin, vmax=vmax, cmap=cmap,
                    marker='o', s=ms_size, linewidths=0.5, edgecolors='k',
                    label='Z nRMSE'
                )
                ax.set_title(f'nRMSE Z (Total = {self.rms_Z:.2f})')
    
            elif plot_type == "VTF":
                print("Plotting nRMSE for VTF...")
                self.mt_coords_VTF = pd.DataFrame(columns=self.mt_coords.columns)
                for rx in range(len(self.ids_VTF)):
                    ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
                    self.mt_coords_VTF = pd.concat([self.mt_coords_VTF, pd.DataFrame([self.mt_coords.loc[ind]])], ignore_index=True)
    
                pc = ax.scatter(
                    self.mt_coords_VTF['east'], self.mt_coords_VTF['north'],
                    c=self.rms_breakDown_VTF['Total'],
                    vmin=vmin, vmax=vmax, cmap=cmap,
                    marker='s', s=ms_size, linewidths=0.5, edgecolors='k',
                    label='VTF nRMSE'
                )
                ax.set_title(f'nRMSE VTF (Total = {self.rms_VTF:.2f})')
    
            elif plot_type == "PT":
                print("Plotting nRMSE for PT...")
                pc = ax.scatter(
                    self.mt_coords['east'], self.mt_coords['north'],
                    c=self.rms_breakDown['Total'],
                    vmin=vmin, vmax=vmax, cmap=cmap,
                    marker='^', s=ms_size, linewidths=0.5, edgecolors='k',
                    label='PT nRMSE'
                )
                ax.set_title("nRMSE PT")
    
            # === Formatting ===
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect('equal')
            ax.grid(lw=0.1)
            ax.set_xlabel("East (km)")
            ax.set_ylabel("North (km)")
            plt.colorbar(pc, ax=ax)
            ax.legend()
    
        plt.tight_layout()
        plt.subplots_adjust(wspace=0)  # reduce horizontal space between subplots
        
        if save_plot:
            if filename is None:
                parts = []
                if Z: parts.append("Z")
                if VTF: parts.append("VTF")
                if PT: parts.append("PT")
                filename = f"rms{'_'.join(parts)}_map.pdf"
    
            os.makedirs(self.plot_dir, exist_ok=True)
            filepath = f"{self.plot_dir}/{filename}"
            plt.savefig(filepath, dpi=300, format='pdf', bbox_inches="tight")
            print(f"Saved: {filepath}")
            plt.close(fig)


    # def plot_rms_map(self,
    #                  ax, 
    #                  Z  = True,
    #                  VTF = False,
    #                  PT = False,
    #                  xlim = [-100, 100],ylim = [-100, 100], 
    #                  vmin=1,vmax=3,
    #                  ms_size = 40,
    #                  save_plot = True):

        
    #     #plt.figure(1,figsize=figsize)
    #     # cmap = cm.get_cmap('bone_r', 12)
    #     cmap = cm.get_cmap('jet',12)
        
    #     if Z:
    #         print('Plotting nRMSE map for the impedance Z...')
    #         pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown_Z['Total'],
    #                 vmin=vmin,vmax=vmax,cmap=cmap,
    #                 marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
    #         ax.set_title('nRMSE Map    Total nRMSE Z = %.2f'%self.rms_Z)
            
    #     if VTF:
    #         print('Plotting nRMSE map for the vertical transfer function VTF...')
            
    #         self.mt_coords_VTF = pd.DataFrame(columns=self.mt_coords.columns)
    #         print(self.mt_coords_VTF)
    #         for rx in range(len(self.ids_VTF)):
    #             ind = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0][0]
    #             #self.mt_coords_VTF = self.mt_coords_VTF.append(self.mt_coords.loc[ind])
    #             self.mt_coords_VTF  = pd.concat([self.mt_coords_VTF, pd.DataFrame([self.mt_coords.loc[ind]])], ignore_index=True)
            
    #         pc=ax.scatter(self.mt_coords_VTF['east'],self.mt_coords_VTF['north'],c=self.rms_breakDown_VTF['Total'],
    #                 vmin=vmin,vmax=vmax,cmap=cmap,
    #                 marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
    #         ax.set_title('nRMSE Map    Total nRMSE VTF = %.2f'%self.rms_VTF)
            
    #     if PT:
    #         pc=ax.scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown['Total'],
    #                 vmin=vmin,vmax=vmax,cmap=cmap,
    #                 marker='o', s=ms_size, linewidths=0.5,edgecolors='k',)
    #     plt.colorbar(pc, ax=ax)
        
    #     #ax.axis('equal')
    #     ax.set_xlim(xlim)
    #     ax.set_ylim(ylim)
    #     ax.grid(lw=0.1)
    #     ax.set_xlabel('East (km)')
    #     ax.set_ylabel('North (km)')
    #     #plt.tight_layout()
        
    #     if save_plot:
    #         if Z:
    #             plt.savefig('%s/rmsZ_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
    #         if VTF:
    #             plt.savefig('%s/rmsVTF_map.png'%(self.plot_dir),dpi=300,bbox_inches='tight')
    

    
    def plot_rms_map_components(self,figsize=(16,3.5), 
                                xlim = [-100, 100],ylim = [-100, 100], 
                                vmin=1,vmax=3,
                                ms_size = 40,
                                save_plot = True):

        print('Plotting nRMSE map ...')
        fig, axs = plt.subplots(1,4,sharey=True, figsize=figsize)
        # cmap = cm.get_cmap('bone_r', 12)
        cmap = cm.get_cmap('jet', 12)

        comps = ['Zxx','Zxy','Zyx','Zyy']
        for i, comp in enumerate(comps):
            
            pc=axs[i].scatter(self.mt_coords['east'],self.mt_coords['north'],c=self.rms_breakDown[comp],
                        vmin=vmin,vmax=vmax,cmap=cmap,
                        marker='o', s=ms_size,linewidths=0.5,edgecolors='k',)
            plt.colorbar(pc, ax=axs[i])
            axs[i].set_title('%s'%comps[i])
            axs[i].axis('equal')
            # plt.grid(lw=0.5)
            # plt.xlabel('East (km)')
            # plt.ylabel('North (km)')
        plt.tight_layout()    
        
        if save_plot:
            plt.savefig('%s/rms_map_comp.pdf'%(self.plot_dir),dpi=300,format='pdf', bbox_inches='tight')
            # plt.close('all')

    
    def plot_Z_fit(self, 
                 plot_Z = False,
                 xlim = [0,5],
                 ylim = [0,4],
                 ylim_diag = [-2,2],
                 save_plot = True,
                 add_map_stats = False):
        
        
        if self.resp_Z is None:
            print("No Z data was inverted!")
            return 
        
        def format_fit(ax, 
                       xlim, 
                       ylim, 
                       xlabel='Log$_{10}$ Period (sec)', 
                       ylabel='Log$_{10}$ $\\rho_{app}$ ($\Omega$m)'):
        
            ax.set_xlim(xlim)  # <<< THIS ensures x-axis is locked
            ax.set_xlabel(xlabel)
            ax.set_xticks(np.arange(xlim[0], xlim[1]))
            ax.set_ylim(ylim[0], ylim[1])
            ax.set_ylabel(ylabel)
            ax.grid(lw=0.3)
            ax.legend(loc='upper right', fontsize='xx-small', borderaxespad=-0.7, framealpha=1)
        
        
        print('Plotting inversion responses ...')
        for site in range(len(self.ids_Z)):
            print('   ...MT site ', self.ids_Z[site])
            
            data = self.resp_Z[self.resp_Z['StaID'] == self.ids_Z[site]]

            f = data['Freq[Hz]']                
            print(f"Site {self.ids_Z[site]}: Frequency values =")
            print(data['Freq[Hz]'].values)
            fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8,6),
                            sharex='all')
            axs[0,0].set_title('%s'%self.ids_Z[site])
            
            if plot_Z: 
                
                self.plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Obs']), data['Re(Zxx)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[0,0],f,abs(data['Im(Zxx)_Obs']), data['Im(Zxx)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[0,0],f,abs(data['Re(Zxx)_Cal']),c='r',obs=False)
                self.plot_Z(axs[0,0],f,abs(data['Im(Zxx)_Cal']),c='b',obs=False)
                format_fit(axs[0,0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Log$_{10}$ Zxx (ohm)')
                axs[0,0].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                
                self.plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Obs']), data['Re(Zxy)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Obs']), data['Im(Zxy)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[0,1],f,abs(data['Re(Zxy)_Cal']),c='r',obs=False)
                self.plot_Z(axs[0,1],f,abs(data['Im(Zxy)_Cal']),c='b',obs=False)
                format_fit(axs[0,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zxy (ohm)')
                axs[0,1].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                    
                self.plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Obs']), data['Re(Zyx)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Obs']), data['Im(Zyx)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[1,0],f,abs(data['Re(Zyx)_Cal']),c='r',obs=False)
                self.plot_Z(axs[1,0],f,abs(data['Im(Zyx)_Cal']),c='b',obs=False)
                format_fit(axs[1,0], xlim =xlim, ylim =ylim, xlabel = '' , ylabel = 'Log$_{10}$ Zyx (ohm)')
                axs[1,0].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                
                self.plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Obs']), data['Re(Zyy)_SD'],c='r',label='Re',obs=True)
                self.plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Obs']), data['Im(Zyy)_SD'],c='b',label='Im',obs=True)
                self.plot_Z(axs[1,1],f,abs(data['Re(Zyy)_Cal']),c='r',obs=False)
                self.plot_Z(axs[1,1],f,abs(data['Im(Zyy)_Cal']),c='b',obs=False)
                format_fit(axs[1,1], xlim =xlim, ylim =ylim , ylabel = 'Log$_{10}$ Zyy (ohm)')
                axs[1,1].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                
                plt.tight_layout()
                
                
                if save_plot:
                    plt.savefig('%s/Z_%s.pdf'%(self.plot_dir,self.ids_Z[site]),dpi=300,format='pdf', bbox_inches='tight')
                    plt.close('all')
                
            else:
 
                rhoXX, phyXX, drhoXX, dphyXX, log_drhoXX = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxx)_Obs'],data['Im(Zxx)_Obs'],data['Re(Zxx)_SD'])
                rhoXY, phyXY, drhoXY, dphyXY, log_drhoXY = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxy)_Obs'],data['Im(Zxy)_Obs'],data['Re(Zxy)_SD'])
                rhoYX, phyYX, drhoYX, dphyYX, log_drhoYX = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyx)_Obs'],data['Im(Zyx)_Obs'],data['Re(Zyx)_SD'])
                rhoYY, phyYY, drhoYY, dphyYY, log_drhoYY = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyy)_Obs'],data['Im(Zyy)_Obs'],data['Re(Zyy)_SD'])
                    
                rhoXX_calc, phyXX_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxx)_Cal'],data['Im(Zxx)_Cal'],data['Re(Zxx)_SD'])
                rhoXY_calc, phyXY_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zxy)_Cal'],data['Im(Zxy)_Cal'],data['Re(Zxy)_SD'])
                rhoYX_calc, phyYX_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyx)_Cal'],data['Im(Zyx)_Cal'],data['Re(Zyx)_SD'])
                rhoYY_calc, phyYY_calc, _, _, _ = self.z2rhophy(data['Freq[Hz]'], data['Re(Zyy)_Cal'],data['Im(Zyy)_Cal'],data['Re(Zyy)_SD'])



                self.plot_rho(axs[0,0],f,rhoXY, log_drhoXY,c='r',label='xy',obs=True)
                self.plot_rho(axs[0,0],f,rhoYX, log_drhoYX,c='b',label='yx',obs=True)
                self.plot_rho(axs[0,0],f,rhoXY_calc,c='r',obs=False)
                self.plot_rho(axs[0,0],f,rhoYX_calc,c='b',obs=False)
                format_fit(axs[0,0], xlim =xlim, ylim = ylim, xlabel = '' , ylabel = 'Log$_{10}$ $\\rho_{app}$ ($\Omega$m)')
                axs[0,0].legend(
                    loc='upper left',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                
                self.plot_phy(axs[1,0],f,phyXY, dphyXY,c='r',label='xy',obs=True, fold_phase=True)
                self.plot_phy(axs[1,0],f,phyYX, dphyYX,c='b',label='yx',obs=True, fold_phase=True)
                self.plot_phy(axs[1,0],f,phyXY_calc,c='r',label='',obs=False, fold_phase=True)
                self.plot_phy(axs[1,0],f,phyYX_calc,c='b',label='',obs=False, fold_phase=True)
                format_fit(axs[1,0], xlim =xlim, ylim = [0,90] , ylabel = 'Phase ($^\circ$)')
                axs[1,0].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                    
                self.plot_rho(axs[0,1],f,rhoXX, log_drhoXX,c='r',label='xx',obs=True)
                self.plot_rho(axs[0,1],f,rhoYY, log_drhoYY,c='b',label='yy',obs=True)
                self.plot_rho(axs[0,1],f,rhoXX_calc,c='r',obs=False)
                self.plot_rho(axs[0,1],f,rhoYY_calc,c='b',obs=False)
                format_fit(axs[0,1], xlim =xlim, ylim = ylim_diag , xlabel = '', ylabel = '')
                axs[0,1].legend(
                    loc='upper left',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                
                self.plot_phy(axs[1,1],f,phyXX, dphyXX,c='r',label='xx',obs=True, fold_phase=False)
                self.plot_phy(axs[1,1],f,phyYY, dphyYY,c='b',label='yy',obs=True, fold_phase=False)
                self.plot_phy(axs[1,1],f,phyXX_calc,c='r',label='',obs=False, fold_phase=False)
                self.plot_phy(axs[1,1],f,phyYY_calc,c='b',label='',obs=False, fold_phase=False)
                format_fit(axs[1,1], xlim =xlim, ylim = [-180,180] , ylabel = '')
                axs[1,1].legend(
                    loc='upper right',
                    fontsize='small',
                    framealpha=1,
                    borderaxespad=0.2
                    )
                plt.tight_layout()
                plt.subplots_adjust(wspace=0.3)  # Increase vertical spacing between columns
                    
            if add_map_stats:
                plt.subplots_adjust(right=0.7)
                ax_map = fig.add_axes([0.75,.55,.2,.3])
                ax_stats = fig.add_axes([0.75,.1,.2,.3])
                
                self.plot_loc_map(self.ids_Z[site], ax_map)
                self.plot_inv_stats_Z(self.ids_Z[site], ax_stats)
  
                
            if save_plot:
                plt.savefig('%s/Z_%s.pdf'%(self.plot_dir,self.ids_Z[site]),dpi=300,format='pdf', bbox_inches='tight')
                plt.close('all')
                

    def plot_loc_map(self, site_id, ax):
        ax.axis('equal')
        ax.set_xticks([]) 
        ax.set_yticks([])
        ax.plot(self.mt_coords['east'], self.mt_coords['north'], 'ko',ms=2)
        coord_site = self.mt_coords[self.mt_coords['id'] == site_id]
        ax.plot(coord_site['east'], coord_site['north'], 'ro',ms=3)
        
    
    def makeloctxt(self):
        #this makes the loc.txt file of station locations required by makeCutawayForFemtic
        loc_for_export = self.mt_coords[['north', 'east', 'z', 'id']]
        loc_for_export.to_csv('loc.txt', sep=' ', index=False, header=False)
        
        
    def plot_inv_stats_Z(self, site_id, ax):
        
        ax.axis("off")
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)
        
        rms_site = self.rms_breakDown_Z[self.rms_breakDown_Z['StaID'] == site_id]
        
        ax.text(0, 9, 'RMS site = %.2f'%rms_site['Total'].values)
        ax.text(0, 7, 'RMS Zxx  = %.2f'%rms_site['Zxx'].values)
        ax.text(0, 5, 'RMS Zxy  = %.2f'%rms_site['Zxy'].values)
        ax.text(0, 3, 'RMS Zyx  = %.2f'%rms_site['Zyx'].values)
        ax.text(0, 1, 'RMS Zyy  = %.2f'%rms_site['Zyy'].values)
 

    def plot_inv_stats_VTF(self, site_id, ax):
        
        ax.axis("off")
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)
        
        rms_site = self.rms_breakDown_VTF[self.rms_breakDown_VTF['StaID'] == site_id]
        
        ax.text(0, 9, 'RMS VTF site = %.2f'%rms_site['Total'].values)
        ax.text(0, 7, 'RMS Tzx  = %.2f'%rms_site['Tzx'].values)
        ax.text(0, 5, 'RMS Tzy  = %.2f'%rms_site['Tzy'].values)
     

    def plot_induction_arrows(self,
                              frequencies,
                              real=True,
                              imag=False,
                              inv_response=False,
                              scale=1,
                              xlim=[-100, 100],
                              ylim=[-100, 100], 
                              save_plot=True):
        
        from scipy.spatial.distance import pdist
        from matplotlib.lines import Line2D
        
        if self.resp_VTF is None:
            print("No VTF data was inverted!")
            return
    
        # Compute scaling factor based on station spacing
        dist_sites = min(pdist(np.array([self.mt_coords['east'], self.mt_coords['north']]).T))
        scale = 0.5 * dist_sites * scale
    
        # Match station coordinates to VTF responses
        self.mt_coords_VTF = pd.DataFrame(columns=self.mt_coords.columns)
        for rx in range(len(self.ids_VTF)):
            matches = np.where(self.ids_VTF[rx] == self.mt_coords['id'].values)[0]
            if len(matches) == 0:
                print(f"Warning: Station ID {self.ids_VTF[rx]} not found in coordinates.")
                continue
            ind = matches[0]
            self.mt_coords_VTF = pd.concat(
                [self.mt_coords_VTF, pd.DataFrame([self.mt_coords.loc[ind]])],
                ignore_index=True
            )
    
        # Determine how many rows to use
        num_rows = int(real) + int(imag)
        row_labels = []
        if real:
            row_labels.append("Re")
        if imag:
            row_labels.append("Im")
    
        fig, axs = plt.subplots(num_rows, len(frequencies), figsize=(3.5* len(frequencies), 4.5 * num_rows), sharex=False, sharey=False)
        
        if num_rows == 1:
            axs = np.atleast_1d(axs)  # flatten if single row
        else:
            axs = np.array(axs)
    
        for fr_idx, freq in enumerate(frequencies):
            ia = self.resp_VTF[self.resp_VTF['Freq[Hz]'] == freq].reset_index()
    
            for row_idx, label in enumerate(row_labels):
                ax = axs[row_idx, fr_idx] if num_rows > 1 else axs[fr_idx]
                ax.scatter(self.mt_coords['east'], self.mt_coords['north'], c='k', marker='o', s=1)
    
                for i in range(len(ia)):
                    station_id = ia['StaID'][i]
                    match = self.mt_coords_VTF[self.mt_coords_VTF['id'] == station_id]
                    if match.empty:
                        continue
    
                    x = match['east'].values[0]
                    y = match['north'].values[0]
    
                    if label == "Re":
                        if real:
                            # Observed
                            ax.arrow(
                                x, y,
                                ia['Re(Tzy)_Obs'][i] * scale,
                                ia['Re(Tzx)_Obs'][i] * scale,
                                width=0.06,
                                head_width=0.2,
                                head_length=0.2,
                                ec='k', fc='k'
                            )
                            # Inverted
                            if inv_response:
                                ax.arrow(
                                    x, y,
                                    ia['Re(Tzy)_Cal'][i] * scale,
                                    ia['Re(Tzx)_Cal'][i] * scale,
                                    width=0.06,
                                    head_width=0.2,
                                    head_length=0.2,
                                    ec='r', fc='r', alpha=0.8
                                )
                            ax.set_title(f'Re Induction Arrows - f = {freq:.3f} Hz')
                            
                    elif label == "Im":
                        if imag:
                            # Observed
                            ax.arrow(
                                x, y,
                                ia['Im(Tzy)_Obs'][i] * scale,
                                ia['Im(Tzx)_Obs'][i] * scale,
                                width=0.06,
                                head_width=0.2,
                                head_length=0.2,
                                ec='k', fc='k'
                            )
                            # Inverted
                            if inv_response:
                                ax.arrow(
                                    x, y,
                                    ia['Im(Tzy)_Cal'][i] * scale,
                                    ia['Im(Tzx)_Cal'][i] * scale,
                                    width=0.06,
                                    head_width=0.2,
                                    head_length=0.2,
                                    ec='r', fc='r', alpha=0.8
                                )
                            ax.set_title(f'Im Induction Arrows - f = {freq:.3f} Hz')
    
                # Plot formatting
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.set_aspect('equal')
                ax.grid(lw=0.1)
                ax.set_xlabel('East (km)')
                ax.set_ylabel('North (km)')
                # Add reference arrow for scale (e.g., unit arrow length = 1)
                ref_length = 0.1 * scale
                ref_x = ax.get_xlim()[1] - 5  # place arrow 0.5 unit right of left x-limit
                ref_y = ax.get_ylim()[0] + 0.5  # place arrow 0.5 unit up from bottom y-limit
                
                # Draw the reference arrow
                ax.arrow(
                    ref_x, ref_y,
                    ref_length, 0,
                    width=0.05,
                    head_width=0.2,
                    head_length=0.3,
                    color='black'
                )
                
                # Label for the reference arrow
                ax.text(
                    ref_x + ref_length + 1,
                    ref_y,
                    f"Length = 0.1",
                    va='center',
                    fontsize=8,
                    color='black'
                )
                # Legend
                legend_elements = [
                    Line2D([0], [0], color='k', lw=2, label='Observed'),
                    Line2D([0], [0], color='r', lw=2, label='Modelled', alpha=0.8)
                ]
                ax.legend(handles=legend_elements, loc='upper left')
    
        if save_plot:
            plt.tight_layout()
            plt.subplots_adjust(wspace=0.1, hspace=0.2)  # reduce spacing between subplots

            if not os.path.exists(self.plot_dir):
                os.makedirs(self.plot_dir)
            plt.savefig(f"{self.plot_dir}/IA.pdf", dpi=300, format='pdf', bbox_inches='tight')
            plt.close('all')
    
         
        
    def plot_VTF_fit(self,   xlim = [-3,2],
                         ylim = [-0.5,0.5],
                         save_plot = True,
                         add_map_stats = False):
        
        
        if self.resp_VTF is None:
            print("No VTF data was inverted!")
            return 
        
        def format_fit(ax, 
                           xlim = [0,5], 
                           ylim = [-1.1,1.1], 
                           xlabel = 'Log$_{10}$ Period (sec)', 
                           ylabel = 'Tzx'):
                
                ax.set_xlabel(xlabel)
                ax.set_xticks(np.arange(xlim[0],xlim[1] ))
                ax.set_xlim(xlim)
                # ax.set_yticks(np.arange(ylim[0], ylim[1]))
                ax.set_ylim(ylim[0], ylim[1])
                ax.set_ylabel(ylabel)
                ax.grid(lw=0.3)
                ax.legend(loc='upper right',fontsize='xx-small',borderaxespad=-0.7,framealpha=1)
        
        
        print('Plotting VTF responses ...')
        for site in range(len(self.ids_VTF)):
            print('   ...MT site ', self.ids_VTF[site])
            
            data = self.resp_VTF[self.resp_VTF['StaID'] == self.ids_VTF[site]]

            f = data['Freq[Hz]']                
            
            fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(5,5),
                            sharex='all')
            axs[0].set_title('%s'%self.ids_VTF[site])
            
            self.plot_tz(axs[0], f, data['Re(Tzx)_Obs'], data['Re(Tzx)_SD'], c='r',label='Re',obs=True)
            self.plot_tz(axs[0], f, data['Im(Tzx)_Obs'], data['Im(Tzx)_SD'], c='b',label='Im',obs=True)
            self.plot_tz(axs[0], f, data['Re(Tzx)_Cal'], c='r',obs=False)
            self.plot_tz(axs[0], f, data['Im(Tzx)_Cal'], c='b',obs=False)
            format_fit(axs[0], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzx')
            axs[0].legend(
                loc='upper right',
                #bbox_to_anchor=(0.98, 0.98),  # x, y in axes fraction coords (1.0 = top/right edge)
                fontsize='small',
                framealpha=1,
                borderaxespad=0.2
                )
            
            self.plot_tz(axs[1], f, data['Re(Tzy)_Obs'], data['Re(Tzy)_SD'], c='r',label='Re',obs=True)
            self.plot_tz(axs[1], f, data['Im(Tzy)_Obs'], data['Im(Tzy)_SD'], c='b',label='Im',obs=True)
            self.plot_tz(axs[1], f, data['Re(Tzy)_Cal'], c='r',obs=False)
            self.plot_tz(axs[1], f, data['Im(Tzy)_Cal'], c='b',obs=False)
            format_fit(axs[1], xlim =xlim, ylim =ylim , xlabel = '', ylabel = 'Tzy')
            axs[1].legend(
                loc='upper right',
                #bbox_to_anchor=(0.98, 0.98),  # x, y in axes fraction coords (1.0 = top/right edge)
                fontsize='small',
                framealpha=1,
                borderaxespad=0.2
                )
            plt.tight_layout()

            if add_map_stats:
                plt.subplots_adjust(right=0.7)
                ax_map = fig.add_axes([0.75,.55,.2,.3])
                ax_stats = fig.add_axes([0.75,.1,.2,.3])
                
                self.plot_loc_map(self.ids_VTF[site], ax_map)
                self.plot_inv_stats_VTF(self.ids_VTF[site], ax_stats)
  
                
            if save_plot:
                plt.savefig('%s/VTF_%s.pdf'%(self.plot_dir,self.ids_VTF[site]),dpi=300,format='pdf', bbox_inches='tight')
                plt.close('all')        
        
    def plot_Z(self, ax, f, Z, dZ=None, c='r', label='xy', obs=True):
        eps = 1e-10  # Small value to avoid log10(0)
    
        if obs:
            # Clip Z, Z+dZ, Z-dZ to ensure valid log10 inputs
            Z_safe = np.clip(Z, eps, None)
            Z_plus = np.clip(Z + dZ, eps, None)
            Z_minus = np.clip(Z - dZ, eps, None)
    
            y = np.log10(Z_safe)
            yerr = np.log10(Z_plus) - np.log10(Z_minus)
    
            ax.errorbar(np.log10(1/f), y, yerr=yerr,
                        fmt=f'{c}.', label=label, zorder=32,
                        elinewidth=0.6, markersize=8,
                        capsize=2, capthick=0.6, mec='k', mew=0.5, alpha=0.5)
        else:
            Z_safe = np.clip(Z, eps, None)
            ax.plot(np.log10(1/f), np.log10(Z_safe), c=c)
    
    
    # def plot_Z(self, ax,f,Z,dZ=None,c='r',label='xy',obs=True):
    #     if obs:
    #         ax.errorbar(np.log10(1/f), np.log10(Z), yerr = np.log10(Z+dZ) - np.log10(Z-dZ), 
    #                         fmt='%s.'%c,label= label,zorder=32, 
    #                         elinewidth=0.6,markersize=8 ,
    #                         capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
    #     else:
    #         ax.plot(np.log10(1/f), np.log10(Z),c=c)


    def plot_Z_VTF_fit(self,
                        plot_Z=False,
                        xlim_Z=[0, 5],
                        ylim_Z=[0, 4],
                        ylim_diag=[-2, 2],
                        xlim_VTF=[-3, 2],
                        ylim_VTF=[-0.5, 0.5],
                        save_plot=True,
                        auto_ylim=True):

        from matplotlib.lines import Line2D
        for site in self.ids_Z:
        
            # Map Z ID → VTF ID
            vtf_site = str(int(site) + 1000)
        
            if vtf_site not in self.ids_VTF:
                print(f"Skipping site {site}: VTF site {vtf_site} not found.")
                continue
        
            data_Z = self.resp_Z[self.resp_Z['StaID'] == site]
            data_VTF = self.resp_VTF[self.resp_VTF['StaID'] == vtf_site]

            f = data_Z['Freq[Hz]']
    
            fig, axs = plt.subplots(2, 4, figsize=(18, 8), gridspec_kw={'width_ratios': [1, 1, 1, 0.5]})
            axs = axs.flatten()
            axs[0].set_title(f"{site}", fontsize=14)
    
            # Z Component: Zxy & Zyx (App. Res + Phase)
            rhoXY, phyXY, drhoXY, dphyXY, log_drhoXY = self.z2rhophy(f, data_Z['Re(Zxy)_Obs'], data_Z['Im(Zxy)_Obs'], data_Z['Re(Zxy)_SD'])
            rhoYX, phyYX, drhoYX, dphyYX, log_drhoYX = self.z2rhophy(f, data_Z['Re(Zyx)_Obs'], data_Z['Im(Zyx)_Obs'], data_Z['Re(Zyx)_SD'])
            rhoXY_calc, phyXY_calc, *_ = self.z2rhophy(f, data_Z['Re(Zxy)_Cal'], data_Z['Im(Zxy)_Cal'], data_Z['Re(Zxy)_SD'])
            rhoYX_calc, phyYX_calc, *_ = self.z2rhophy(f, data_Z['Re(Zyx)_Cal'], data_Z['Im(Zyx)_Cal'], data_Z['Re(Zyx)_SD'])
    
            rhoXX, phyXX, drhoXX, dphyXX, log_drhoXX = self.z2rhophy(f, data_Z['Re(Zxx)_Obs'], data_Z['Im(Zxx)_Obs'], data_Z['Re(Zxx)_SD'])
            rhoYY, phyYY, drhoYY, dphyYY, log_drhoYY = self.z2rhophy(f, data_Z['Re(Zyy)_Obs'], data_Z['Im(Zyy)_Obs'], data_Z['Re(Zyy)_SD'])
            rhoXX_calc, phyXX_calc, *_ = self.z2rhophy(f, data_Z['Re(Zxx)_Cal'], data_Z['Im(Zxx)_Cal'], data_Z['Re(Zxx)_SD'])
            rhoYY_calc, phyYY_calc, *_ = self.z2rhophy(f, data_Z['Re(Zyy)_Cal'], data_Z['Im(Zyy)_Cal'], data_Z['Re(Zyy)_SD'])
    
            offdiag_rhos = np.concatenate([rhoXY, rhoYX, rhoXY_calc, rhoYX_calc])
            diag_rhos = np.concatenate([rhoXX, rhoYY, rhoXX_calc, rhoYY_calc])
    
            def compute_ylim(log_vals):
                log_vals = np.log10(log_vals[log_vals > 0])
                min_decade = int(np.floor(log_vals.min()))
                max_decade = int(np.ceil(log_vals.max()))
                range_decades = max(2, min(8, max_decade - min_decade))
                return [min_decade, min_decade + range_decades]
    
            if auto_ylim:
                ylim_rho_offdiag = compute_ylim(offdiag_rhos)
                ylim_rho_diag = compute_ylim(diag_rhos)
            else:
                ylim_rho_offdiag = ylim_Z
                ylim_rho_diag = ylim_diag
    
            self.plot_rho(axs[0], f, rhoXY, log_drhoXY, 'r', 'xy', True)
            self.plot_rho(axs[0], f, rhoYX, log_drhoYX, 'b', 'yx', True)
            self.plot_rho(axs[0], f, rhoXY_calc, obs=False, c='r')
            self.plot_rho(axs[0], f, rhoYX_calc, obs=False, c='b')
            axs[0].set_ylabel('Log10 Rho ($Ω$m)', fontsize=12)
            axs[0].set_ylim(ylim_rho_offdiag)
            axs[0].set_xlim(xlim_Z)
            axs[0].legend(handles=[
                Line2D([0], [0], color='r', label='ρxy Obs/Cal'),
                Line2D([0], [0], color='b', label='ρyx Obs/Cal')
            ], loc='upper left', fontsize=11)
    
            self.plot_rho(axs[1], f, rhoXX, log_drhoXX, 'r', 'xx', True)
            self.plot_rho(axs[1], f, rhoYY, log_drhoYY, 'b', 'yy', True)
            self.plot_rho(axs[1], f, rhoXX_calc, obs=False, c='r')
            self.plot_rho(axs[1], f, rhoYY_calc, obs=False, c='b')
            axs[1].set_ylabel('Log10 Rho ($Ω$m)', fontsize=12)
            axs[1].set_ylim(ylim_rho_diag)
            axs[1].set_xlim(xlim_Z)
            axs[1].legend(handles=[
                Line2D([0], [0], color='r', label='ρxx Obs/Cal'),
                Line2D([0], [0], color='b', label='ρyy Obs/Cal')
            ], loc='upper left', fontsize=11)
    
            f_vtf = data_VTF['Freq[Hz]']
            self.plot_tz(axs[2], f_vtf, data_VTF['Re(Tzx)_Obs'], data_VTF['Re(Tzx)_SD'], 'r', 'Re', True)
            self.plot_tz(axs[2], f_vtf, data_VTF['Im(Tzx)_Obs'], data_VTF['Im(Tzx)_SD'], 'b', 'Im', True)
            self.plot_tz(axs[2], f_vtf, data_VTF['Re(Tzx)_Cal'], c='r', obs=False)
            self.plot_tz(axs[2], f_vtf, data_VTF['Im(Tzx)_Cal'], c='b', obs=False)
            axs[2].set_ylabel('Tzx', fontsize=12)
            axs[2].set_ylim(ylim_VTF)
            axs[2].set_xlim(xlim_VTF)
            axs[2].legend(handles=[
                Line2D([0], [0], color='r', label='Tzx Re Obs/Cal'),
                Line2D([0], [0], color='b', label='Tzx Im Obs/Cal')
            ], loc='upper left', fontsize=11)
    
            self.plot_loc_map(site, axs[3])
    
            self.plot_phy(axs[4], f, phyXY, dphyXY, 'r', 'xy', True, fold_phase=True)
            self.plot_phy(axs[4], f, phyYX, dphyYX, 'b', 'yx', True, fold_phase=True)
            self.plot_phy(axs[4], f, phyXY_calc, obs=False, c='r', fold_phase=True)
            self.plot_phy(axs[4], f, phyYX_calc, obs=False, c='b', fold_phase=True)
            axs[4].set_ylabel('Phase ($^\circ$)', fontsize=12)
            axs[4].set_ylim([0, 90])
            axs[4].set_xlim(xlim_Z)
            axs[4].legend(handles=[
                Line2D([0], [0], color='r', label='φxy Obs/Cal'),
                Line2D([0], [0], color='b', label='φyx Obs/Cal')
            ], loc='upper right', fontsize=11)
    
            self.plot_phy(axs[5], f, phyXX, dphyXX, 'r', 'xx', True, fold_phase=False)
            self.plot_phy(axs[5], f, phyYY, dphyYY, 'b', 'yy', True, fold_phase=False)
            self.plot_phy(axs[5], f, phyXX_calc, obs=False, c='r', fold_phase=False)
            self.plot_phy(axs[5], f, phyYY_calc, obs=False, c='b', fold_phase=False)
            axs[5].set_ylabel('Phase ($^\circ$)', fontsize=12)
            axs[5].set_ylim([-180, 180])
            axs[5].set_xlim(xlim_Z)
            axs[5].legend(handles=[
                Line2D([0], [0], color='r', label='φxx Obs/Cal'),
                Line2D([0], [0], color='b', label='φyy Obs/Cal')
            ], loc='upper right', fontsize=11)
    
            self.plot_tz(axs[6], f_vtf, data_VTF['Re(Tzy)_Obs'], data_VTF['Re(Tzy)_SD'], 'r', 'Re', True)
            self.plot_tz(axs[6], f_vtf, data_VTF['Im(Tzy)_Obs'], data_VTF['Im(Tzy)_SD'], 'b', 'Im', True)
            self.plot_tz(axs[6], f_vtf, data_VTF['Re(Tzy)_Cal'], c='r', obs=False)
            self.plot_tz(axs[6], f_vtf, data_VTF['Im(Tzy)_Cal'], c='b', obs=False)
            axs[6].set_ylabel('Tzy', fontsize=12)
            axs[6].set_ylim(ylim_VTF)
            axs[6].set_xlim(xlim_VTF)
            axs[6].legend(handles=[
                Line2D([0], [0], color='r', label='Tzy Re Obs/Cal'),
                Line2D([0], [0], color='b', label='Tzy Im Obs/Cal')
            ], loc='upper left', fontsize=11)
    
            axs[7].axis('off')
            z_row = self.rms_breakDown_Z[self.rms_breakDown_Z['StaID'] == site]
            vtf_row = self.rms_breakDown_VTF[self.rms_breakDown_VTF['StaID'] == vtf_site]
            site_rms_total = self.rms_total_by_site.get(site, np.nan)
            table_data = [
                ["Zxx", f"{z_row['Zxx'].values[0]:.2f}"],
                ["Zxy", f"{z_row['Zxy'].values[0]:.2f}"],
                ["Zyx", f"{z_row['Zyx'].values[0]:.2f}"],
                ["Zyy", f"{z_row['Zyy'].values[0]:.2f}"],
                ["Tzx", f"{vtf_row['Tzx'].values[0]:.2f}"],
                ["Tzy", f"{vtf_row['Tzy'].values[0]:.2f}"],
                ["Site total", f"{site_rms_total:.2f}"]
            ]
            table = axs[7].table(cellText=table_data,
                     colLabels=["Component", "RMS"],
                     loc='center', cellLoc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(12)
            
            for key, cell in table.get_celld().items():
                cell.set_height(0.1)

            
            for i, ax in enumerate(axs):
                ax.grid(True)
                if i != 3 and i!=7:  # Skip map axis
                    ax.set_xlabel('Log10 Period (s)', fontsize=12)
                elif i==3:
                    ax.set_xlabel('Site map', fontsize=12)
                ax.tick_params(labelsize=11)
    
            fig.tight_layout()
            if save_plot:
                plt.savefig(f"{self.plot_dir}/ZVTF_{site}.pdf", dpi=300, format='pdf', bbox_inches='tight')
                plt.close(fig)

            

    # def plot_Z_VTF_fit(self,
    #                     plot_Z=False,
    #                     xlim_Z=[0, 5],
    #                     ylim_Z=[0, 4],
    #                     ylim_diag=[-2, 2],
    #                     xlim_VTF=[-3, 2],
    #                     ylim_VTF=[-0.5, 0.5],
    #                     save_plot=True,
    #                     auto_ylim=True):
    
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from matplotlib.lines import Line2D
    
    #     for site in self.ids_Z:
    #         if site not in self.ids_VTF:
    #             print(f"Skipping {site}: not found in VTF data.")
    #             continue
    
    #         data_Z = self.resp_Z[self.resp_Z['StaID'] == site]
    #         data_VTF = self.resp_VTF[self.resp_VTF['StaID'] == site]
    #         f = data_Z['Freq[Hz]']
    
    #         fig, axs = plt.subplots(2, 4, figsize=(18, 8))
    #         axs = axs.flatten()
    #         axs[0].set_title(f"{site}")
    
    #         # Z Component: Zxy & Zyx (App. Res + Phase)
    #         rhoXY, phyXY, drhoXY, dphyXY, log_drhoXY = self.z2rhophy(f, data_Z['Re(Zxy)_Obs'], data_Z['Im(Zxy)_Obs'], data_Z['Re(Zxy)_SD'])
    #         rhoYX, phyYX, drhoYX, dphyYX, log_drhoYX = self.z2rhophy(f, data_Z['Re(Zyx)_Obs'], data_Z['Im(Zyx)_Obs'], data_Z['Re(Zyx)_SD'])
    #         rhoXY_calc, phyXY_calc, *_ = self.z2rhophy(f, data_Z['Re(Zxy)_Cal'], data_Z['Im(Zxy)_Cal'], data_Z['Re(Zxy)_SD'])
    #         rhoYX_calc, phyYX_calc, *_ = self.z2rhophy(f, data_Z['Re(Zyx)_Cal'], data_Z['Im(Zyx)_Cal'], data_Z['Re(Zyx)_SD'])
    
    #         rhoXX, phyXX, drhoXX, dphyXX, log_drhoXX = self.z2rhophy(f, data_Z['Re(Zxx)_Obs'], data_Z['Im(Zxx)_Obs'], data_Z['Re(Zxx)_SD'])
    #         rhoYY, phyYY, drhoYY, dphyYY, log_drhoYY = self.z2rhophy(f, data_Z['Re(Zyy)_Obs'], data_Z['Im(Zyy)_Obs'], data_Z['Re(Zyy)_SD'])
    #         rhoXX_calc, phyXX_calc, *_ = self.z2rhophy(f, data_Z['Re(Zxx)_Cal'], data_Z['Im(Zxx)_Cal'], data_Z['Re(Zxx)_SD'])
    #         rhoYY_calc, phyYY_calc, *_ = self.z2rhophy(f, data_Z['Re(Zyy)_Cal'], data_Z['Im(Zyy)_Cal'], data_Z['Re(Zyy)_SD'])
    
    #         offdiag_rhos = np.concatenate([rhoXY, rhoYX, rhoXY_calc, rhoYX_calc])
    #         diag_rhos = np.concatenate([rhoXX, rhoYY, rhoXX_calc, rhoYY_calc])
    
    #         def compute_ylim(log_vals):
    #             log_vals = np.log10(log_vals[log_vals > 0])
    #             min_decade = int(np.floor(log_vals.min()))
    #             max_decade = int(np.ceil(log_vals.max()))
    #             range_decades = max(2, min(4, max_decade - min_decade))
    #             return [min_decade, min_decade + range_decades]
    
    #         if auto_ylim:
    #             ylim_rho_offdiag = compute_ylim(offdiag_rhos)
    #             ylim_rho_diag = compute_ylim(diag_rhos)
    #         else:
    #             ylim_rho_offdiag = ylim_Z
    #             ylim_rho_diag = ylim_diag
    
    #         self.plot_rho(axs[0], f, rhoXY, log_drhoXY, 'r', 'xy', True)
    #         self.plot_rho(axs[0], f, rhoYX, log_drhoYX, 'b', 'yx', True)
    #         self.plot_rho(axs[0], f, rhoXY_calc, obs=False, c='r')
    #         self.plot_rho(axs[0], f, rhoYX_calc, obs=False, c='b')
    #         axs[0].set_ylabel('Log10 Rho ($Ω$m)')
    #         axs[0].set_ylim(ylim_rho_offdiag)
    #         axs[0].set_xlim(xlim_Z)
    #         axs[0].legend(handles=[
    #             Line2D([0], [0], color='r', label='ρxy Obs/Cal'),
    #             Line2D([0], [0], color='b', label='ρyx Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         self.plot_rho(axs[1], f, rhoXX, log_drhoXX, 'r', 'xx', True)
    #         self.plot_rho(axs[1], f, rhoYY, log_drhoYY, 'b', 'yy', True)
    #         self.plot_rho(axs[1], f, rhoXX_calc, obs=False, c='r')
    #         self.plot_rho(axs[1], f, rhoYY_calc, obs=False, c='b')
    #         axs[1].set_ylabel('Log10 Rho ($Ω$m)')
    #         axs[1].set_ylim(ylim_rho_diag)
    #         axs[1].set_xlim(xlim_Z)
    #         axs[1].legend(handles=[
    #             Line2D([0], [0], color='r', label='ρxx Obs/Cal'),
    #             Line2D([0], [0], color='b', label='ρyy Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         f_vtf = data_VTF['Freq[Hz]']
    #         self.plot_tz(axs[2], f_vtf, data_VTF['Re(Tzx)_Obs'], data_VTF['Re(Tzx)_SD'], 'r', 'Re', True)
    #         self.plot_tz(axs[2], f_vtf, data_VTF['Im(Tzx)_Obs'], data_VTF['Im(Tzx)_SD'], 'b', 'Im', True)
    #         self.plot_tz(axs[2], f_vtf, data_VTF['Re(Tzx)_Cal'], c='r', obs=False)
    #         self.plot_tz(axs[2], f_vtf, data_VTF['Im(Tzx)_Cal'], c='b', obs=False)
    #         axs[2].set_ylabel('Tzx')
    #         axs[2].set_ylim(ylim_VTF)
    #         axs[2].set_xlim(xlim_VTF)
    #         axs[2].legend(handles=[
    #             Line2D([0], [0], color='r', label='Tzx Re Obs/Cal'),
    #             Line2D([0], [0], color='b', label='Tzx Im Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         self.plot_loc_map(site, axs[3])
    
    #         self.plot_phy(axs[4], f, phyXY, dphyXY, 'r', 'xy', True)
    #         self.plot_phy(axs[4], f, phyYX, dphyYX, 'b', 'yx', True)
    #         self.plot_phy(axs[4], f, phyXY_calc, obs=False, c='r')
    #         self.plot_phy(axs[4], f, phyYX_calc, obs=False, c='b', )
    #         axs[4].set_ylabel('Phase ($^\circ$)')
    #         axs[4].set_ylim([-180, 180])
    #         axs[4].set_xlim(xlim_Z)
    #         axs[4].legend(handles=[
    #             Line2D([0], [0], color='r', label='φxy Obs/Cal'),
    #             Line2D([0], [0], color='b', label='φyx Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         self.plot_phy(axs[5], f, phyXX, dphyXX, 'r', 'xx', True, fold_phase=False)
    #         self.plot_phy(axs[5], f, phyYY, dphyYY, 'b', 'yy', True, fold_phase=False)
    #         self.plot_phy(axs[5], f, phyXX_calc, obs=False, c='r', fold_phase=False)
    #         self.plot_phy(axs[5], f, phyYY_calc, obs=False, c='b', fold_phase=False)
    #         axs[5].set_ylabel('Phase ($^\circ$)')
    #         axs[5].set_ylim([-180, 180])
    #         axs[5].set_xlim(xlim_Z)
    #         axs[5].legend(handles=[
    #             Line2D([0], [0], color='r', label='φxx Obs/Cal'),
    #             Line2D([0], [0], color='b', label='φyy Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         self.plot_tz(axs[6], f_vtf, data_VTF['Re(Tzy)_Obs'], data_VTF['Re(Tzy)_SD'], 'r', 'Re', True)
    #         self.plot_tz(axs[6], f_vtf, data_VTF['Im(Tzy)_Obs'], data_VTF['Im(Tzy)_SD'], 'b', 'Im', True)
    #         self.plot_tz(axs[6], f_vtf, data_VTF['Re(Tzy)_Cal'], c='r', obs=False)
    #         self.plot_tz(axs[6], f_vtf, data_VTF['Im(Tzy)_Cal'], c='b', obs=False)
    #         axs[6].set_ylabel('Tzy')
    #         axs[6].set_ylim(ylim_VTF)
    #         axs[6].set_xlim(xlim_VTF)
    #         axs[6].legend(handles=[
    #             Line2D([0], [0], color='r', label='Tzy Re Obs/Cal'),
    #             Line2D([0], [0], color='b', label='Tzy Im Obs/Cal')
    #         ], loc='upper left', fontsize='small')
    
    #         axs[7].axis('off')
    #         z_row = self.rms_breakDown_Z[self.rms_breakDown_Z['StaID'] == site]
    #         vtf_row = self.rms_breakDown_VTF[self.rms_breakDown_VTF['StaID'] == site]
    #         rms_total = self.rms_total if hasattr(self, 'rms_total') else np.nan
    #         table_data = [
    #             ["Zxx", z_row['Zxx'].values[0]],
    #             ["Zxy", z_row['Zxy'].values[0]],
    #             ["Zyx", z_row['Zyx'].values[0]],
    #             ["Zyy", z_row['Zyy'].values[0]],
    #             ["Tzx", vtf_row['Tzx'].values[0]],
    #             ["Tzy", vtf_row['Tzy'].values[0]],
    #             ["Total", rms_total]
    #         ]
    #         axs[7].table(cellText=table_data,
    #                     colLabels=["Component", "RMS"],
    #                     loc='center', cellLoc='center')
    
    #         for ax in axs:
    #             ax.grid(True)
    #             ax.set_xlabel('Log10 Period (s)')
    
    #         fig.tight_layout()
    #         if save_plot:
    #             plt.savefig(f"{self.plot_dir}/ZVTF_{site}.png", dpi=300, bbox_inches='tight')
    #             plt.close(fig)




            
    
    def plot_rho(self, ax,f,rho,log_drho=None,c='r',label='xy',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), np.log10(rho), yerr = log_drho, 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
            # Shaded envelope
            if log_drho is not None:
                ax.fill_between(np.log10(1/f), np.log10(rho) - log_drho, np.log10(rho) + log_drho, color=c, alpha=0.1)
            
        else:
            ax.plot(np.log10(1/f), np.log10(rho),c=c)
    

    def plot_phy(self, ax, f, phy, dphy=None, c='r', label='xy', obs=True, fold_phase=False):
        
        def fold_phase_to_0_90(phase_array):
            """
            Map phase values from [-180, 180] into [0, 90] degrees.
            """
            phase_mod = np.mod(phase_array, 360)        # convert to [0, 360]
            phase_folded = np.abs((phase_mod + 90) % 180 - 90)  # fold into [0, 90]
            return phase_folded
        
        if fold_phase:
            phy = fold_phase_to_0_90(phy)
            if dphy is not None:
                # Approximate transformation of the error bar — not mathematically rigorous
                dphy = np.abs(dphy)  # Keep error bars positive in folded range
    
        if obs:
            ax.errorbar(np.log10(1/f), phy, yerr=dphy,
                        fmt='%s.' % c, label=label, zorder=32,
                        elinewidth=0.6, markersize=8,
                        capsize=2, capthick=0.6, mec='k', mew=0.5, alpha=0.5)
            if dphy is not None:
                ax.fill_between(np.log10(1/f), phy - dphy, phy + dphy, color=c, alpha=0.1)
        else:
            ax.plot(np.log10(1/f), phy, c=c, label=label)

    
    
    # def plot_phy(self, ax,f,phy,dphy=None,c='r',label='xy',obs=True):
    #     if obs:
    #         ax.errorbar(np.log10(1/f), phy, yerr = dphy, 
    #                         fmt='%s.'%c,label= label,zorder=32, 
    #                         elinewidth=0.6,markersize=8 ,
    #                         capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
    #         # Shaded envelope
    #         if dphy is not None:
    #             ax.fill_between(np.log10(1/f), phy - dphy, phy + dphy, color=c, alpha=0.1)

    #     else:
    #         ax.plot(np.log10(1/f), phy,c=c)
    

    # def plot_tz(self, ax, f, tz, dtz=None, c='r', label='tzx', obs=True):
    #     x = np.log10(1 / f)
        
    #     if obs:
    #         # Plot central line
    #         ax.plot(x, tz, f'{c}o', label=label, markersize=6, alpha=0.7, markerfacecolor='none')
            
    #         # Shaded envelope
    #         if dtz is not None:
    #             ax.fill_between(x, tz - dtz, tz + dtz, color=c, alpha=0.2, label=f"{label} ±SD")
        
    #     else:
    #         ax.plot(x, tz, color=c, linestyle='-', label=f"{label} model")
    
    def plot_tz(self, ax,f,tz,dtz=None,c='r',label='tzx',obs=True):
        if obs:
            ax.errorbar(np.log10(1/f), tz, yerr = dtz, 
                            fmt='%s.'%c,label= label,zorder=32, 
                            elinewidth=0.6,markersize=8 ,
                            capsize=3,capthick=0.6,mec='k',mew=0.5, alpha=0.5)
            # Shaded envelope
            if dtz is not None:
                ax.fill_between(np.log10(1/f), tz - dtz, tz + dtz, color=c, alpha=0.1)
        else:
            ax.plot(np.log10(1/f), tz,c=c)  


