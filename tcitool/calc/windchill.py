import tcitool

class WindChill_JAGTICalculator(tcitool.Calculator):
    def __init__(self,tool):
        super().__init__(tool)
        self.export_params = {'wcet':'wcet_jagti'}
        self.require_data('t2m','ws10')

    def main(self):
        t2mC = self.tool.data.get('t2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['t2m']))
        wind_at_15dm = (3.6*self.tool.data['ws10'])**0.16
        wcet = 13.12 + 0.6215 * t2mC \
                     - 11.37 * wind_at_15dm \
                     + 0.3965 * t2mC * wind_at_15dm
        wcet.attrs.update({
            'units': 'deg C',
            'long_name': 'Wind Chill Equivalent Temperature (using JAG/TI '
                'calculation method)',
            'source': 'Eq. 2.3 in https://cdn.knmi.nl/system/downloads/files/000/000/016/original/gevoelstemperatuur.pdf?1433939065'
        })
        self.data['wcet'] = wcet
