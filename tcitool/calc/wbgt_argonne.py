import numpy as np
import pandas as pd
import xarray as xr
import dask.diagnostics as daskdiag
import scipy.optimize


import tcitool
import tcitool.func as tf

class WBGT_ArgonneCalculator(tcitool.OptimizationCalculator):
    def __init__(self,tool,**kwargs):
        super().__init__(tool)
        self.export_params = {'wbgt':'wbgt_argonne',
                              'tg':'tg_argonne',
                              'tnw':'tnw_argonne'}
        self.tool.require_data('t2m','skt','rh','e_kPa','msl_kPa',
            'ws10','fsr','Isw_in','Ibeam','solza','soldist','fal')

        self.hyperparams = kwargs
        self.fn = {}
        self.const = {}

    def preface(self):
        self.default_hyperparams()
        self.constants()
        self.initial_calculations()
        self.functions()

    def default_hyperparams(self):
        if 'Tg_lim' not in self.hyperparams:
            self.hyperparams['Tg_lim'] = (-60,120)
        if 'Tnw_lim' not in self.hyperparams:
            self.hyperparams['Tnw_lim'] = (-60,90)
        if 'xtol' not in self.hyperparams:
            self.hyperparams['xtol'] = 0.01
        if 'calcTpsy' not in self.hyperparams:
            self.hyperparams['calcTpsy'] = False
        if 'useTg150' not in self.hyperparams:
            self.hyperparams['useTg150'] = True
        if 'daytime' not in self.hyperparams:
            self.hyperparams['daytime'] = False
        if 'cza_min' not in self.hyperparams:
            self.hyperparams['cza_min'] = np.cos(np.deg2rad(87.5))
        if 'height' not in self.hyperparams:
            self.hyperparams['height'] = 2
        if 'min_ws' not in self.hyperparams:
            self.hyperparams['min_ws'] = self.tool.options.get(
                'windspeed_lowlimit',0.1)

    def constants(self):
        self.const['GRAVITY'] = tf.td.GRAVITATIONAL_AC
        self.const['STEFANB'] = tf.m.STEFAN_BOLTZMANN
        self.const['SOLAR_CONST'] = 1367.
        self.const['Cp'] = 1003.5
        self.const['M_AIR'] = tf.td.MOLAR_MASS_AIR
        self.const['M_H2O'] = tf.td.MOLAR_MASS_H2O
        self.const['RATIO'] = ( self.const['Cp']*self.const['M_AIR']/self.const['M_H2O'] )
        self.const['R_GAS'] = tf.td.GAS_CONSTANT
        self.const['R_AIR'] = ( self.const['R_GAS']/self.const['M_AIR'] )
        # define wick constants
        self.const['EMIS_WICK'] = 0.95
        self.const['ALB_WICK'] = 0.4
        self.const['D_WICK'] = 0.007
        self.const['L_WICK'] = 0.0254
        # define globe constants
        self.const['EMIS_GLOBE'] = 0.95
        self.const['ALB_GLOBE'] = 0.05
        self.const['D_GLOBE'] = 0.0508
        # define surface constants
        self.const['EMIS_SFC'] = 0.999
        # define computational and physical limits
        self.const['CZA_MIN'] = self.hyperparams['cza_min']
        self.const['MIN_SPEED'] = self.hyperparams['min_ws']

    def initial_calculations(self):
        self.data = self.tool.data[['t2m','skt','rh','e_kPa','msl_kPa','fal']]

        solza_fill = xr.where(
            self.tool.data['solza']>1.57079615,
            np.deg2rad(90),
            self.tool.data['solza']
        )
        self.data['solcza'] = np.cos(solza_fill)
        self.data['solcza'].attrs = {'units':'',
            'long_name':'cosine of solar zenith angle'}

        self.data['Isw_in'] = xr.where(
            self.data['solcza']>self.const['CZA_MIN'],
            self.tool.data['Isw_in'],0)
        self.data['Isw_in'].attrs = {'units':'W/m2',
            'long_name':'solar irradiance'}
        self.data['Isw_frac'] = np.clip(
            xr.where(self.data['Isw_in']<1,0,
                self.tool.data['Ibeam']/self.data['Isw_in']),0,0.9)
        self.data['ws'] = np.clip(tf.m.wind_at_height_using_fsr(
            self.tool.data['ws10'],
            self.tool.data['fsr'],
            self.hyperparams['height']
        ),self.const['MIN_SPEED'],None)
        self.data['ws'].attrs = {
            'units': self.tool.data['ws10'].attrs['units'],
            'long_name': ('Estimated windspeed at %3.1fm' %
                self.hyperparams['height'])
        }
        self.data = self.data[
                ['t2m','skt','rh','e_kPa','msl_kPa','ws','Isw_in',
                'Isw_frac','solcza','fal']
            ].rename_vars({'msl_kPa':'P_kPa'})
        self.data = self.data.persist()

    def functions(self):
        self.fn = {
            'tref': lambda temp1,temp2: 0.5*(temp1+temp2),
            'esat': tf.m.saturated_vapor_pressure,
            'evap': tf.td.heat_of_evaporation,
            'emis': lambda temp_K, rh: tf.td.atmospheric_emissivity(
                rh*tf.m.saturated_vapor_pressure(temp_K)),
            'dens': tf.td.density,
            'diff': tf.td.thermal_diffusivity,
            'visc': tf.td.dynamic_viscosity,
            'cond': tf.td.thermal_conductivity,
            'Sc': tf.td.schmidt_number,
            'Pr': tf.td.prandtl_number
        }

        self.fn['BG_Re'] = lambda temp_K, P_kPa, ws: (
            ws * self.fn['dens'](temp_K,P_kPa) * self.const['D_GLOBE'] /
                self.fn['visc'](temp_K))
        self.fn['BG_Nu'] = lambda temp_K, P_kPa, ws: (
            2 + 0.6 * np.sqrt(self.fn['BG_Re'](temp_K, P_kPa, ws)) *
                np.power(self.fn['Pr'](temp_K),0.3333))
        self.fn['BG_h'] = lambda temp_K, P_kPa, ws: (
            self.fn['BG_Nu'](temp_K, P_kPa, ws) * self.fn['cond'](temp_K) *
                self.const['D_GLOBE'])
        self.fn['h_sphere_in_air'] = lambda Tg, t2m, P_kPa, ws: (
            self.fn['BG_h'](self.fn['tref'](Tg, t2m), P_kPa, ws))
        self.fn['Tg'] = ( lambda Tg, t2m, skt, rh, e_kPa, P_kPa, ws, Isw_in,
            Isw_frac, solcza, fal:
            0.5 * self.fn['emis'](self.fn['tref'](Tg, t2m),rh) * np.power(t2m,4.)
            + 0.5 * self.const['EMIS_SFC'] * np.power(skt,4.)
            - (self.fn['h_sphere_in_air'](Tg,t2m,P_kPa,ws) / (
                self.const['STEFANB'] * self.const['EMIS_GLOBE']) * (Tg - t2m) )
            + (Isw_in / (2. * self.const['STEFANB'] * self.const['EMIS_GLOBE']))
                * (1. - self.const['ALB_GLOBE'])
                * (Isw_frac*(1./(2.*solcza)-1.)+1.+fal) \
            - np.power(Tg,4) )

        self.fn['Fatm'] = ( lambda Tnw, t2m, skt, rh, Isw_in, Isw_frac, solcza,
            fal:
            self.const['STEFANB'] * self.const['EMIS_WICK'] * (
                0.5*( self.fn['emis'](t2m,rh)*np.power(t2m,4.)
                    + self.const['EMIS_SFC']*np.power(skt,4.))
                - np.power(Tnw,4.))
            + (1. - self.const['ALB_WICK']) * Isw_in * (
                (1. - Isw_frac) * (
                    1. + 0.25 * self.const['D_WICK'] / self.const['L_WICK'])
                + Isw_frac * (
                    (np.tan(np.arccos(solcza))/np.pi)
                    + 0.25 * self.const['D_WICK'] / self.const['L_WICK'])
                + fal) )
        self.fn['NW_Re'] = lambda temp_K, P_kPa, ws: (
            ws * self.fn['dens'](temp_K,P_kPa) * self.const['D_WICK'] /
                self.fn['visc'](temp_K) )
        self.fn['NW_Nu'] = lambda temp_K, P_kPa, ws: (
            0.281 * np.power(self.fn['NW_Re'](temp_K, P_kPa, ws),1-0.4) *
                np.power(self.fn['Pr'](temp_K),1-0.56) )
        self.fn['NW_h'] = lambda temp_K, P_kPa, ws: (
            self.fn['NW_Nu'](temp_K, P_kPa, ws) * self.fn['cond'](temp_K) /
                self.const['D_WICK'] )
        self.fn['h_cylinder_in_air'] = lambda Tnw, t2m, P_kPa, ws: (
            self.fn['NW_h'](self.fn['tref'](Tnw, t2m), P_kPa, ws) )

        self.fn['wetfrac'] = lambda Tnw, e_kPa, P_kPa: (
            (self.fn['esat'](Tnw) - e_kPa)/(P_kPa - self.fn['esat'](Tnw)) )
        self.fn['Tnw'] = ( lambda Tnw, t2m, skt, rh, e_kPa, P_kPa, ws,
            Isw_in, Isw_frac, solcza, fal:
            Tnw - (t2m - self.fn['evap'](self.fn['tref'](Tnw,t2m))
                / self.const['RATIO']
                * self.fn['wetfrac'](Tnw, rh*self.fn['esat'](t2m), P_kPa)
                * np.power(self.fn['Pr'](self.fn['tref'](Tnw, t2m))
                    / self.fn['Sc'](self.fn['tref'](Tnw, t2m),P_kPa),0.56)
                + self.fn['Fatm'](Tnw,t2m,skt,rh,Isw_in,Isw_frac,solcza,fal)
                    / self.fn['h_cylinder_in_air'](Tnw, t2m, P_kPa, ws)) )

    def optimize_params(self):
        return self.data[['t2m','skt','rh','e_kPa','P_kPa','ws','Isw_in',
            'Isw_frac','solcza','fal']]

    def optimize_globe_temperature(self,params):
        param_dict = dict(params)
        param_list = ['t2m','skt','rh','e_kPa','P_kPa','ws','Isw_in',
            'Isw_frac','solcza','fal']
        params_tuple = tuple([param_dict[p] for p in param_list])

        tg = np.nan
        if self.hyperparams['daytime'] and (
                params['solcza']<=self.const['CZA_MIN'] or
                params['Isw_in']<=1 or params['Isw_frac']<=0.01):
            return np.nan
        try:
            tg = scipy.optimize.brentq(
                self.fn['Tg'],
                tf.u.tempC2K(self.hyperparams['Tg_lim'][0]),
                tf.u.tempC2K(self.hyperparams['Tg_lim'][1]),
                xtol=self.hyperparams['xtol'],
                args=params_tuple,
                disp=False,
                full_output=False)
        except ValueError as e:
            self.postface()
            if err.args[0] == 'f(a) and f(b) must have different signs':
                olist = []
                for t in range(self.hyperparams['Tg_lim'][0],
                        self.hyperparams['Tg_lim'][1]+1,10):
                    olist.append('    %d: %f'%(
                        t,self.fn['Tg'](tf.u.tempC2K(t),*params_tuple)))
                msg = ('\nNo Tg could be found, that satisfies the conditions. '
                    '\nTg ranges between %+06.1f °C and %+06.1f °C'
                    '(tollerance %05.3f °C).'
                    '\nParameters: %s \nOptimizor limits: \n')
                msg = msg % (
                    tf.u.tempC2K(self.hyperparams['Tg_lim'][0]),
                    tf.u.tempC2K(self.hyperparams['Tg_lim'][1]),
                    self.hyperparams['xtol'],
                    str(dict(params)),
                ) + '\n'.join(olist)
                e.args = tuple(err.args[0]+msg,)
            raise
        return tg

    def optimize_natural_wetbulb_temperature(self,params):
        param_dict = dict(params)
        param_list = ['t2m','skt','rh','e_kPa','P_kPa','ws','Isw_in',
            'Isw_frac','solcza','fal']
        params_tuple = tuple([param_dict[p] for p in param_list])

        tnw = np.nan
        if self.hyperparams['daytime'] and (
                params['solcza']<=self.const['CZA_MIN'] or
                params['Isw_in']<=1 or params['Isw_frac']<=0.01):
            return np.nan
        try:
            tnw = scipy.optimize.brentq(
                self.fn['Tnw'],
                tf.u.tempC2K(self.hyperparams['Tnw_lim'][0]),
                tf.u.tempC2K(self.hyperparams['Tnw_lim'][1]),
                xtol=self.hyperparams['xtol'],
                args=params_tuple,
                disp=False,
                full_output=False)
        except ValueError as e:
            self.postface()
            if err.args[0] == 'f(a) and f(b) must have different signs':
                olist = []
                for t in range(self.hyperparams['Tg_lim'][0],
                        self.hyperparams['Tg_lim'][1]+1,10):
                    olist.append('    %d: %f'%(
                        t,self.fn['Tg'](tf.u.tempC2K(t),*params_tuple)))
                msg = ('\nNo Tnw could be found, that satisfies the conditions. '
                    '\nTnw ranges between %+06.1f °C and %+06.1f °C'
                    '(tollerance %05.3f °C).'
                    '\nParameters: %s \nOptimizor limits: \n')
                msg = msg % (
                    tf.u.tempC2K(self.hyperparams['Tnw_lim'][0]),
                    tf.u.tempC2K(self.hyperparams['Tnw_lim'][1]),
                    self.hyperparams['xtol'],
                    str(dict(params)),
                ) + '\n'.join(olist)
                e.args = tuple(err.args[0]+msg,)
            raise
        return tnw

    def optimize(self):
        dask = False
        shape = self.tool.data.shape
        if self.tool.data.shape != self.tool.data.chunks:
            dask = True
            df = self.optimize_params().to_dask_dataframe(
                dim_order=['time','latitude','longitude'])
            df['tg_5cm'] = df.apply(self.optimize_globe_temperature,
                axis=1,meta=('tg_5cm','f8'))
            df['tnw'] = df.apply(self.optimize_natural_wetbulb_temperature,
                axis=1,meta=('tnw','f8'))
            pbar = daskdiag.ProgressBar()
            pbar.register()
            df['tg_5´cm'] = df['tg_5cm'].compute()
            df['tnw'] = df['tnw'].compute()
            pbar.unregister()
            tg_arr = df['tg_5cm'].to_dask_array(True)
            tnw_arr = df[['tnw']].to_dask_array(True)
            tg_arr = tg_arr.compute_chunk_sizes()
            tnw_arr = tnw_arr.compute_chunk_sizes()
        else:
            df = self.optimize_params().to_dataframe()
            df = df.reorder_levels(['time','latitude','longitude'])
            df['tg_5cm'] = df.apply(self.optimize_globe_temperature,
                axis=1)
            df['tnw'] = df.apply(self.optimize_natural_wetbulb_temperature,
                axis=1)
            tg_arr = df['tg_5cm'].to_numpy()
            tnw_arr = df['tnw'].to_numpy()
        tg_arr = tg_arr.reshape(*shape)
        tnw_arr = tnw_arr.reshape(*shape)
        self.data['tg_5cm'] = ['time','latitude','longitude'], tg_arr
        self.data['tnw'] = ['time','latitude','longitude'], tnw_arr

    def postface(self):
        self.closing_calculations()
        self.export()

    def closing_calculations(self):
        tg_5cmC = tf.u.tempK2C(self.data['tg_5cm'])
        tg_5cmC.attrs.update({'units':'deg C'})
        t2mC = tf.u.tempK2C(self.tool.data['t2m'])
        t2mC.attrs.update({'units':'deg C'})
        tgC = (
            tg_5cmC + (
                1 + 1.13 * np.power(self.data['ws'], 0.6)
                * np.power(self.const['D_GLOBE']*1000, -0.4) * (tg_5cmC-t2mC)
            ) / (1 + 2.41 * np.power(self.data['ws'], 0.6) ) )
        self.data['tg'] = tf.u.tempC2K(tgC)
        self.data['wbgt'] = (0.1 * self.data['t2m'] + 0.2 * self.data['tg']
            + 0.7 * self.data['tnw'])

        self.data['wbgt'].attrs = {'units':'K',
            'long_name': 'Wet Bulb Globe Temperature (using the Argonne model)',
            'source': 'Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp,'
                ' en R. Sharp, 2008: Modeling the Wet Bulb Globe Temperature '
                'Using Standard Meteorological Measurements. Journal of '
                'Occupational and Environmental Hygiene, 5 (10), 645–655, '
                'doi: 10.1080/15459620802310770'
        }
        self.data['tg'].attrs = {'units':'K',
            'long_name': 'Globe temperature (using the Argonne model), corrected'
                'for a globe with a diameter of 5 cm',
            'source': 'Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp,'
                ' en R. Sharp, 2008: Modeling the Wet Bulb Globe Temperature '
                'Using Standard Meteorological Measurements. Journal of '
                'Occupational and Environmental Hygiene, 5 (10), 645–655, '
                'doi: 10.1080/15459620802310770',
            'source_diameter_correction': 'ISO 7243:2017, 2017: Ergonomics of '
                'the thermal environment – Assessment of heat stress using the '
                'WBGT (wet bulb globe temperature) index. Standard, '
                'International Organization for Standardization, Geneva, CH.'
        }
        self.data['tnw'].attrs = {'units':'K',
            'long_name': 'Natural Wet Bulb temperature (using the Argonne '
                'model)',
            'source': 'Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp,'
                ' en R. Sharp, 2008: Modeling the Wet Bulb Globe Temperature '
                'Using Standard Meteorological Measurements. Journal of '
                'Occupational and Environmental Hygiene, 5 (10), 645–655, '
                'doi: 10.1080/15459620802310770'
        }
