#!/usr/bin/env -S python -u
import os
import time
from contextlib import contextmanager
import argparse

import numpy as np
import scipy.optimize

class ArgonneModel(object):
    def __init__(self,**kwargs):
        self.const = self.constants()
        self.fn, self.bgfn, self.nwfn = self.functions()
        self.hyperparams = {
            'Tg_lim': (-60,120),
            'Tnw_lim': (-60,90),
            'xtol': 0.1,
            'daytime': False,
            'no_oob': True
        }
        self.hyperparams.update(kwargs)
    def constants(self):
        const = {}
        const['GRAVITY'] = 9.81
        const['STEFANB'] = 5.67e-8
        const['SOLAR_CONST'] = 1367.
        const['Cp'] = 1003.5
        const['M_AIR'] = 28.965
        const['M_H2O'] = 18.015
        const['RATIO'] = ( const['Cp']*const['M_AIR']/const['M_H2O'] )
        const['R_GAS'] = 8314.34
        const['R_AIR'] = ( const['R_GAS']/const['M_AIR'] )
        # define wick constants
        const['EMIS_WICK'] = 0.95
        const['ALB_WICK'] = 0.4
        const['D_WICK'] = 0.007
        const['L_WICK'] = 0.0254
        # define globe constants
        const['EMIS_GLOBE'] = 0.95
        const['ALB_GLOBE'] = 0.05
        const['D_GLOBE'] = 0.0508
        # define surface constants
        const['EMIS_SFC'] = 0.999
        # define computational and physical limits
        const['CZA_MIN'] = np.cos(np.deg2rad(87.5))
        const['MIN_SPEED'] = 0.1
        return const
    def functions(self):
        const = self.const
        fn = {}
        fn['tref'] = lambda temp1,temp2: 0.5*(temp1+temp2)
        fn['esat'] = lambda temp_K: 0.611 * np.exp(17.2694*(temp_K-273.16)/(temp_K-35.86))
        fn['evap'] = lambda temp_K: (313.15-temp_K)/30*-71100+2.4073e6
        fn['emis_e'] = lambda e_kPa: 0.575 * np.power(e_kPa,0.143)
        fn['emis'] = lambda temp_K, rh: fn['emis_e'](rh*fn['esat'](temp_K))
        fn['dens'] = lambda temp_K, P_kPa: (P_kPa*1e3)/((const['R_GAS']/const['M_AIR'])*temp_K)
        fn['cond'] = lambda temp_K: 0.02624 * np.power(temp_K/300,0.8646)
        fn['capp'] = lambda temp_K: 1002.5+275e-6 * np.power(temp_K-200,2)
        fn['vhca'] = lambda temp_K, P_kPa: fn['dens'](temp_K,P_kPa) * fn['capp'](temp_K)
        fn['diff'] = lambda temp_K, P_kPa: fn['cond'](temp_K) / fn['vhca'](temp_K,P_kPa)
        fn['visc'] = lambda temp_K: 1.458e-6*(temp_K**1.5) / (temp_K +110.4)
        fn['Sc'] = lambda temp_K, P_kPa: fn['visc'](temp_K) / (fn['dens'](temp_K,P_kPa) * fn['diff'](temp_K,P_kPa))
        fn['Pr'] = lambda temp_K: fn['capp'](temp_K) * fn['visc'](temp_K) / fn['cond'](temp_K)

        bgfn = {}
        bgfn['Re'] = lambda temp_K, P_kPa, ws: \
            ws * fn['dens'](temp_K,P_kPa) * const['D_GLOBE'] / self.fn['visc'](temp_K)
        bgfn['Nu'] = lambda temp_K, P_kPa, ws: \
            2 + 0.6 * np.sqrt(bgfn['Re'](temp_K, P_kPa, ws)) * np.power(fn['Pr'](temp_K),0.3333)
        bgfn['h'] = lambda temp_K, P_kPa, ws: \
            bgfn['Nu'](temp_K, P_kPa, ws) * fn['cond'](temp_K) * const['D_GLOBE']
        bgfn['hs'] = lambda Tg, t2m, P_kPa, ws: bgfn['h'](fn['tref'](Tg, t2m), P_kPa, ws)
        bgfn['Tg'] = ( lambda Tg, t2m, skt, rh, e_kPa, P_kPa, ws, Isw_in, Isw_frac, solcza, fal:
            0.5 * fn['emis'](fn['tref'](Tg, t2m),rh) * np.power(t2m,4.)
            + 0.5 * const['EMIS_SFC'] * np.power(skt,4.)
            - (bgfn['hs'](Tg,t2m,P_kPa,ws) / (const['STEFANB'] * const['EMIS_GLOBE']) * (Tg - t2m) )
            + (Isw_in / (2. * const['STEFANB'] * const['EMIS_GLOBE']))
                * (1. - const['ALB_GLOBE']) * (Isw_frac*(1./(2.*solcza)-1.)+1.+fal) \
            - np.power(Tg,4) )

        nwfn = {}
        nwfn['Fatm'] = ( lambda Tnw, t2m, skt, rh, Isw_in, Isw_frac, solcza, fal:
            const['STEFANB'] * const['EMIS_WICK'] * (
                0.5*( fn['emis'](t2m,rh)*np.power(t2m,4.) + const['EMIS_SFC']*np.power(skt,4.))
                - np.power(Tnw,4.))
            + (1. - const['ALB_WICK']) * Isw_in * (
                (1. - Isw_frac) * (1. + 0.25 * const['D_WICK'] / const['L_WICK'])
                + Isw_frac * ( (np.tan(np.arccos(solcza))/np.pi) + 0.25 * const['D_WICK'] / const['L_WICK'])
                + fal)
            )
        nwfn['Re'] = lambda temp_K, P_kPa, ws: \
            ws * fn['dens'](temp_K,P_kPa) * const['D_WICK'] / fn['visc'](temp_K)
        nwfn['Nu'] = lambda temp_K, P_kPa, ws: \
            0.281 * np.power(nwfn['Re'](temp_K, P_kPa, ws),1-0.4) * np.power(fn['Pr'](temp_K),1-0.56)
        nwfn['h'] = lambda temp_K, P_kPa, ws: \
            nwfn['Nu'](temp_K, P_kPa, ws) * fn['cond'](temp_K) / const['D_WICK']
        nwfn['hc'] = lambda Tnw, t2m, P_kPa, ws: nwfn['h'](fn['tref'](Tnw, t2m), P_kPa, ws)
        nwfn['wf'] = lambda Tnw, t2m, rh, P_kPa: (fn['esat'](Tnw) - (rh*fn['esat'](t2m)))/(P_kPa - fn['esat'](Tnw))
        nwfn['Tnw'] = ( lambda Tnw, t2m, skt, rh, e_kPa, P_kPa, ws, Isw_in, Isw_frac, solcza, fal:
            Tnw - ( t2m - (
                fn['evap'](fn['tref'](Tnw,t2m))
                / const['RATIO']
                * nwfn['wf'](Tnw, t2m, rh, P_kPa)
                * np.power( fn['Pr'](fn['tref'](Tnw, t2m)) / fn['Sc'](fn['tref'](Tnw, t2m),P_kPa),0.56)
                + ( nwfn['Fatm'](Tnw,t2m,skt,rh,Isw_in,Isw_frac,solcza,fal) / nwfn['hc'](Tnw, t2m, P_kPa, ws) )
            ) ) )
        return fn, bgfn, nwfn
    def optimize_globe_temperature(self,params):
        param_list = ['t2m','skt','rh','e_kPa','P_kPa','ws','Isw_in','Isw_frac','solcza','fal']
        params_tuple = tuple(params)[-len(param_list):]

        tg = np.nan
        if self.hyperparams['daytime'] and (
                params_tuple[6]<=1 or #Isw_in
                params_tuple[7]<=0.01 or #Isw_in
                params_tuple[8]<=self.const['CZA_MIN']): #solcza
            return np.nan
        if len(params_tuple)<len(param_list):
            raise ValueError('Not enough parameters')
        if self.hyperparams['no_oob']:
            a = self.bgfn['Tg'](self.hyperparams['Tg_lim'][0]+273.15,*params_tuple)
            b = self.bgfn['Tg'](self.hyperparams['Tg_lim'][1]+273.15,*params_tuple)
            if np.sign(a)==np.sign(b):
                return np.nan
        try:
            tg = scipy.optimize.brentq(
                f=self.bgfn['Tg'],
                a=self.hyperparams['Tg_lim'][0]+273.15,
                b=self.hyperparams['Tg_lim'][1]+273.15,
                xtol=self.hyperparams['xtol'],
                args=params_tuple,
                full_output=False,
                disp=False)
        except ValueError as err:
            if err.args[0] == 'f(a) and f(b) must have different signs':
                raise self.errormessage(err,'Tg',zip(param_list,params_tuple)) from err
            raise
        return tg
    def optimize_natural_wetbulb_temperature(self,params):
        param_list = ['t2m','skt','rh','e_kPa','P_kPa','ws','Isw_in','Isw_frac','solcza','fal']
        params_tuple = tuple(params)[-len(param_list):]

        tnw = np.nan
        if self.hyperparams['daytime'] and (
                params_tuple[6]<=1 or #Isw_in
                params_tuple[7]<=0.01 or #Isw_in
                params_tuple[8]<=self.const['CZA_MIN']): #solcza
            return np.nan
        if len(params_tuple)<len(param_list):
            raise ValueError('Not enough parameters')
        if self.hyperparams['no_oob']:
            a = self.nwfn['Tnw'](self.hyperparams['Tnw_lim'][0]+273.15,*params_tuple)
            b = self.nwfn['Tnw'](self.hyperparams['Tnw_lim'][1]+273.15,*params_tuple)
            if np.sign(a)==np.sign(b):
                return np.nan
        try:
            tnw = scipy.optimize.brentq(
                f=self.nwfn['Tnw'],
                a=self.hyperparams['Tnw_lim'][0]+273.15,
                b=self.hyperparams['Tnw_lim'][1]+273.15,
                xtol=self.hyperparams['xtol'],
                args=params_tuple,
                full_output=False,
                disp=False)
        except ValueError as err:
            if err.args[0] == 'f(a) and f(b) must have different signs':
                raise self.errormessage(err,'Tnw',zip(param_list,params_tuple)) from err
            raise
        return tnw
    def errormessage(self,err,tg_or_tnw,param):
        param = dict(param)
        msg = (
            '\nNo {tg_or_tnw} could be found, that satisfies the conditions.'
            '\n{tg_or_tnw} ranges between {lim[0]:06.1f}°C and {lim[1]:06.1f} °C (tollerance {xtol:06.1f}°C).'
            '\nParameters: {p} \nOptimizor limits: \n')
        msg = msg.format(**{
        	'tg_or_tnw': tg_or_tnw,
        	'lim': self.hyperparams[tg_or_tnw+'_lim'],
        	'xtol': self.hyperparams['xtol'],
        	'p': str(param) })
        olist = []
        if tg_or_tnw=='Tg':
            for t in range(self.hyperparams['Tg_lim'][0],self.hyperparams['Tg_lim'][1],10):
                olist.append('    %d: %f'%( t,self.bgfn['Tg'](t+273.15,*tuple(param.values())) ))
        elif tg_or_tnw=='Tnw':
            for t in range(self.hyperparams['Tnw_lim'][0],self.hyperparams['Tnw_lim'][1],10):
                olist.append('    %d: %f'%( t,self.nwfn['Tnw'](t+273.15,*tuple(param.values())) ))
        return ValueError(msg+'\n'.join(olist))

@contextmanager
def timeit(premsg='',postmsg='',verbose=False):
    vp(verbose,premsg,flush=True)
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    timems = format(int(elapsedTime * 1000),',d').replace(',',' ')
    vp(verbose,postmsg+' finished in {} ms'.format(timems),flush=True)

def vp(verbose,msg,**kwargs):
    if verbose:
    	print(msg,**kwargs)

def main(infile,outfile,verbose=True,pid=''):
    vp(verbose,'[%s ] Starting. Using %s as data'%(pid,infile),flush=True)
    model = ArgonneModel()
    data = np.load(infile,mmap_mode='r')
    datasize = format(data.size//10,',d').replace(',',' ')
    with timeit(
            premsg='[%sw] Calculating %s natural wetbulb temperatures...'%(pid,datasize),
            postmsg='[%sw] Done,'%(pid),
            verbose=verbose):
        tnwa = np.apply_along_axis(func1d=model.optimize_natural_wetbulb_temperature,axis=0,arr=data)
    with timeit(
            premsg='[%sg] Calculating %s globe temperatures... '%(pid,datasize),
            postmsg='[%sg] Done,'%(pid),
            verbose=verbose):
    	tga = np.apply_along_axis(func1d=model.optimize_globe_temperature,axis=0,arr=data)
    outdata = np.stack([tga,tnwa],axis=0)
    np.save(outfile,outdata)
    del data, datasize, tga, tnwa, outdata, model
    vp(verbose,'[%s ] Stored output to %s'%(pid,outfile),flush=True)

def is_valid_input_file(parser,arg):
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist!"%arg)
    elif not arg.endswith('.npy'):
        parser.error("The file %s is not a *.npy file!"%arg)
    else:
        return arg

def is_valid_output_file(parser,arg):
    dirname = os.path.dirname(arg)
    if not os.path.isdir(dirname):
        parser.error("The directory %s does not exist!"%dirname)
    elif not os.access(dirname, os.W_OK):
        parser.error("The directory %s is not writeable!"%dirname)
    else:
        return arg

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the Argonne model on a Numpy Array file (*.npy).')
    parser.add_argument('-i', dest='inputfilename', required=True,
        type=lambda x: is_valid_input_file(parser, x),
        help='The input array as *.npy\n'
        'Different input variables need to be in the first dimension (axis=0) and the order '
        '2m temp, surface temp, rel hum, vapor pressure, pressure, 2m wind speed, global radiation, '
        'fraction direct/global radiation, cosine of solar zenith angle, albedo')
    parser.add_argument('-o', dest='outputfilename', required=True,
        type=lambda x: is_valid_output_file(parser, x),
        help='The output array as *.npy\n'
        'Different output variables will be in the first dimension (axis=0) and the order '
        'Globe temperature (5cm globe), Natural wet bulb temperature')
    parser.add_argument('-p', dest='pid', help='String to separate the output from this script, from others', action='store')
    parser.add_argument('-v', dest='verbose', help='Print more data', action='store_true')
    args = parser.parse_args()
    main(args.inputfilename,args.outputfilename,True,args.pid)
