import tcitool

class CommonMeteoGenerators(object):
    @classmethod
    def register_generators(cls,gr):
        gr.register(cls.t2mC, 't2mC', 't2m')
        gr.register(cls.d2mC, 'd2mC', 'd2m')
        gr.register(cls.sktC, 'sktC', 'skt')
        gr.register(cls.ws10, 'ws10', ['u10','v10'])
        gr.register(cls.wdir10, 'wdir10', ['u10','v10'])

        gr.register(cls.dewpoint, 'd2m', 'e_kPa')
        gr.register(cls.from_rh, 'd2m', ['t2m','rh'])
        gr.register(cls.e, 'e_kPa', 'd2m')
        gr.register(cls.from_rh, 'e_kPa', ['t2m','rh'])
        gr.register(cls.e_sat, 'e_sat_kPa', 't2m')
        gr.register(cls.rh, 'rh', ['t2m','e_kPa'])
        gr.register(cls.rh, 'rh', ['t2m','d2m'])

        gr.register(cls.pressure_kPa, ['sp_kPa','msl_kPa'], ['sp','msl'])
        gr.register(cls.pressure_Pa,  ['sp','msl'], ['sp_kPa','msl_kPa'])

    @classmethod
    def t2mC(cls,tool):
        """Converts the 2m temperature from K to deg C"""
        tool.data['t2mC'] = tcitool.UnitFuncs.tempK2C(tool.data['t2m'])
        tool.data['t2mC'].attrs = tool.data['t2m'].attrs
        tool.data['t2mC'].attrs['units'] = 'deg C'

    @classmethod
    def d2mC(cls,tool):
        """Converts the 2m dew point from K to deg C"""
        tool.data['d2mC'] = tcitool.UnitFuncs.tempK2C(tool.data['d2m'])
        tool.data['d2mC'].attrs = tool.data['d2m'].attrs
        tool.data['d2mC'].attrs['units'] = 'deg C'

    @classmethod
    def sktC(cls,tool):
        """Converts the skin temperature from K to deg C"""
        tool.data['sktC'] = tcitool.UnitFuncs.tempK2C(tool.data['skt'])
        tool.data['sktC'].attrs = tool.data['skt'].attrs
        tool.data['sktC'].attrs['units'] = 'deg C'

    @classmethod
    def ws10(cls,tool):
        tool.data['ws10'] = tcitool.MeteoFuncs.wind_speed(
            tool.data['u10'],tool.data['v10'])
        tool.data['ws10'].attrs = tool.data['u10'].attrs
        tool.data['ws10'].attrs['long_name'] = 'Wind speed at 10 metre'
    @classmethod
    def wdir10(cls,tool):
        tool.data['wdir10'] = tcitool.MeteoFuncs.wind_direction(
            tool.data['u10'],tool.data['v10'])
        tool.data['wdir10'].attrs.update({
            'units': 'rad',
            'long_name': 'Wind direction at 10 metre'
        })

    @classmethod
    def dewpoint(cls,tool):
        tool.data['d2m'] = tcitool.MeteoFuncs.dewpoint(tool.data['e_kPa'])
        tool.data['d2m'].attrs.update({
            'units': 'K',
            'long_name': 'Dew point'
        })
    @classmethod
    def e_sat(cls,tool):
        tool.data['e_sat_kPa'] = tcitool.MeteoFuncs.saturated_vapor_pressure(
            tool.data['t2m'])
        tool.data['e_sat_kPa'].attrs.update({
            'units': 'kPa',
            'long_name': 'Saturated vapor pressure'
        })
    @classmethod
    def e(cls,tool):
        tool.data['e_kPa'] = tcitool.MeteoFuncs.saturated_vapor_pressure(
            tool.data['d2m'])
        tool.data['e_kPa'].attrs.update({
            'units': 'kPa',
            'long_name': 'Vapor pressure'
        })
    @classmethod
    def rh(cls,tool):
        """Calculates relative humidity from the temperature and dew point"""
        if 'e_sat_kPa' not in tool.data and 't2m' in tool.data:
            cls.e_sat(tool)
        if 'e_kPa' not in tool.data and 'd2m' in tool.data:
            cls.e(tool)
        tool.data['rh'] = tool.data['e_kPa']/tool.data['e_sat_kPa']
        tool.data['rh'].attrs.update({
            'units': '(0 - 1)',
            'long_name': 'Relative Humidity'
        })
    @classmethod
    def from_rh(cls,tool):
        """Calculates dewpoint and vapor_pressure from temp and rh"""
        if 'e_sat_kPa' not in tool.data and 't2m' in tool.data:
            cls.e_sat(tool)
        tool.data['e_kPa'] = tool.data['e_sat_kPa'] * tool.data['rh']
        tool.data['e_kPa'].attrs.update({
            'units': 'kPa',
            'long_name': 'Vapor pressure'
        })
        cls.dewpoint(tool)

    @classmethod
    def pressure_kPa(cls,tool):
        msl = tool.data['msl']/1000
        sp = tool.data['sp']/1000
        tool.data['msl_kPa'] = msl.assign_attrs(units='kPa',
            long_name=tool.data['msl'].attrs['long_name'])
        tool.data['sp_kPa'] = sp.assign_attrs(units='kPa',
            long_name=tool.data['sp'].attrs['long_name'])
        tool.data.ds = tool.data.ds.drop_vars(['msl','sp'])

    @classmethod
    def pressure_Pa(cls,tool):
        msl = tool.data['msl_kPa']*1000
        sp = tool.data['sp_kPa']*1000
        tool.data['msl'] = msl.assign_attrs(units='Pa',
            long_name=tool.data['msl_kPa'].attrs['long_name'])
        tool.data['sp'] = sp.assign_attrs(units='Pa',
            long_name=tool.data['sp_kPa'].attrs['long_name'])
        tool.data.ds = tool.data.ds.drop_vars(['msl_kPa','sp_kPa'])

class IntegratedVarsGenerators(object):
    @classmethod

    def register_generators(cls,gr):
        gr.register(cls.Isw_in,'Isw_in','ssrd','radiation_integration_time')
        gr.register(cls.Ilw_in,'Ilw_in','strd','radiation_integration_time')
        gr.register(cls.Ibeam,'Ibeam','fdir','radiation_integration_time')

    def cumulatives2regular(cls,data,axis):
        """Helper-function which converts cumulative metrics to averages

        Args:
            data: An np.array containing cumulative values over a certain axis
            axis: an interger describing the axis over which the data is
                cumulative, usualy the time axis.
        Returns:
            A np.array containing the de-cumulatived data
                (averaged over a period)
        """
        prev = np.roll(data,1,axis)
        slices = [(slice(0,1,None)
                   if i==axis
                   else slice(None)) for i in range(data.ndims)]
        prev[slices] = np.nan
        return data-prev

    @classmethod
    def Isw_in(cls,tool):
        """Calculates the incomming shortwave radiation in W/m²."""
        ssrd = tool.data['ssrd']
        if ('radiation_cumulative' in tool.options and
            tool.options['radiation_cumulative']):
            ssrd = cls.cumulatives2regular(ssrd, 0)
        tool.data['Isw_in'] = ssrd/tool.options['radiation_integration_time']
        tool.data['Isw_in'].attrs = {'units':'W m**-2',
            'long_name':'Surface solar irradiation downwards'}

    @classmethod
    def Ilw_in(cls,tool):
        """Calculates the incomming longwave radiation in W/m²."""
        strd = tool.data['strd']
        if ('radiation_cumulative' in tool.options and
            tool.options['radiation_cumulative']):
            strd = cls.cumulatives2regular(strd,0)
        tool.data['Ilw_in'] =  strd/tool.options['radiation_integration_time']
        tool.data['Ilw_in'].attrs = {'units':'W m**-2',
            'long_name':'Surface thermal irradiation downwards'}

    @classmethod
    def Ibeam(cls,tool):
        """Calculates the incomming direct shortwave radiation in W/m²."""
        fdir = tool.data['fdir']
        if ('radiation_cumulative' in tool.options and
            tool.options['radiation_cumulative']):
            fdir = cls.cumulatives2regular(fdir,0)
        tool.data['Ibeam'] =  fdir/tool.options['radiation_integration_time']
        tool.data['Ibeam'].attrs = {'units':'W m**-2',
            'long_name':'Total sky direct solar radiation at surface'}
