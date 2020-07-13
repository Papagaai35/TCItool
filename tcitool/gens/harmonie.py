import numpy as np
import xarray as xr

import tcitool

class HarmonieGenerators(object):
    @classmethod
    def register_generators(cls,gr):
        gr.register(cls.globrad, 'Isw_in', 'grad', 'radiation_integration_time')
        gr.register(cls.netrad, 'Isw_net', 'nswrs', 'radiation_integration_time')
        gr.register(cls.fracrad, 'Isw_frac', ['grad','tcc'], 'radiation_integration_time')
        gr.register(cls.fracrad, 'Ibeam', ['grad','tcc'], 'radiation_integration_time')
        gr.register(cls.albedo, 'fal', ['grad','nswrs'], 'radiation_integration_time')
        gr.register(cls.wind, 'ws2', 'ws10')
        gr.register(cls.wind, 'ws2', ['u10','v10'])

    ws10 = tcitool.CommonMeteoGenerators.ws10
    cumulatives2regular = tcitool.IntegratedVarsGenerators.cumulatives2regular

    @classmethod
    def globrad(cls,tool):
        """Calculates the incomming shortwave radiation in W/m²."""
        grad = tool.data['grad']
        if ('radiation_cumulative' in tool.options and
            tool.options['radiation_cumulative']):
            grad = cls.cumulatives2regular(grad, 0)
        tool.data['Isw_in'] = grad/tool.options['radiation_integration_time']
        tool.data['Isw_in'].attrs = {'units':'W m**-2',
            'long_name':'Surface solar irradiation downwards'}

    @classmethod
    def netrad(cls,tool):
        """Calculates the incomming shortwave radiation in W/m²."""
        nswrs = tool.data['nswrs']
        if ('radiation_cumulative' in tool.options and
            tool.options['radiation_cumulative']):
            nswrs = cls.cumulatives2regular(nswrs, 0)
        tool.data['Isw_net'] = nswrs/tool.options['radiation_integration_time']
        tool.data['Isw_net'].attrs = {'units':'W m**-2',
            'long_name':'Net short-wave radiation flux'}

    @classmethod
    def fracrad(cls,tool):
        if 'Isw_in' not in tool.data and 'grad' in tool.data:
            cls.globrad(tool)
        kd = xr.where(tool.data['tcc']<=0.22,1-0.09*tool.data['tcc'],
            xr.where(tool.data['tcc']<=0.8,0.9511-0.1604*tool.data['tcc']+4.39*(tool.data['tcc']**2)-16.64*(tool.data['tcc']**3),0.165))
        kd.attrs = {'units':'-','long_name':'Fraction diffuse/total radiation'}

        tool.data['Isw_frac'] = xr.where(tool.data['Isw_in']<1,0,1-kd).clip(0,0.9)
        tool.data['Isw_frac'].attrs = {'units':'-','long_name':'Fraction direct/total radiation'}
        tool.data['Ibeam'] = tool.data['Isw_in']*tool.data['Isw_frac']
        tool.data['Ibeam'].attrs = {'units':'W m-2','long_name':'Direct normal irradiance',
            'code': 84, 'table': 253, 'institution': 'KNMI', 'source': 'calculated'}

    @classmethod
    def albedo(cls,tool):
        if 'Isw_in' not in tool.data and 'grad' in tool.data:
            cls.globrad(tool)
        if 'Isw_net' not in tool.data and 'nswrs' in tool.data:
            cls.netrad(tool)
        al = (tool.data['Isw_in']-tool.data['Isw_net'])/tool.data['Isw_in']
        al = al.median(dim='time').expand_dims({'time':tool.data['time'].size})
        tool.data['fal'] = al.dims, al.values
        tool.data['fal'].attrs = {'long_name': 'Albedo', 'units': '-',
            'code': 84, 'table': 253, 'institution': 'KNMI', 'source': 'calculated'}

    @classmethod
    def wind(cls,tool):
        if 'ws10' not in tool.data and 'u10' in tool.data and 'v10' in tool.data:
            cls.ws10(tool)
        tool.data['ws2'] = tool.data['ws10'] * ((2/10)**0.28)
        tool.data['ws2'].attrs = {'long_name': 'Wind speed at 2m', 'units': tool.data['u10'].attrs['units'],
            'source': 'calculated'}
