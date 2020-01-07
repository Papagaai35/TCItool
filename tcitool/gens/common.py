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
        gr.register(cls.e, 'e_kPa', 'd2m')

    @classmethod
    def t2mC(cls,tool):
        """Converts the 2m temperature from K to deg C"""
        tool.data['t2mC'] = tcitool.UnitFuncs.tempK2C(tool.data['t2m'])
        tool.data['t2mC'].attrs['units'] = 'deg C'

    @classmethod
    def d2mC(cls,tool):
        """Converts the 2m dew point from K to deg C"""
        tool.data['d2mC'] = tcitool.UnitFuncs.tempK2C(tool.data['d2m'])
        tool.data['d2mC'].attrs['units'] = 'deg C'

    @classmethod
    def sktC(cls,tool):
        """Converts the skin temperature from K to deg C"""
        tool.data['sktC'] = tcitool.UnitFuncs.tempK2C(tool.data['skt'])
        tool.data['sktC'].attrs['units'] = 'deg C'

    @classmethod
    def ws10(cls,tool):
        tool.data['ws10'] = tcitool.MeteoFuncs.wind_speed(
            tool.data['u10'],tool.data['v10'])
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
    def e(cls,tool):
        tool.data['e_kPa'] = tcitool.MeteoFuncs.saturated_vapor_pressure(
            tool.data['d2m'])
        tool.data['e_kPa'].attrs.update({
            'units': 'kPa',
            'long_name': 'Vapor pressure'
        })
