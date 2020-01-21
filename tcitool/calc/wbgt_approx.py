import numpy as np
import xarray as xr
import tcitool

class WBGTapprox_ACSMCalculator(tcitool.Calculator):
    def __init__(self,tool):
        super().__init__(tool)
        self.export_params = {'wbgt':'wbgt_acsm'}
        self.require_data('t2m','d2m')

    def main(self):
        t2mC = self.tool.data.get('t2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['t2m']))
        d2mC = self.tool.data.get('d2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['d2m']))
        vapor_pressure = 6.112 * np.exp((17.67*d2mC)/(d2mC+243.5))
        wbgt = 0.567 * t2mC + 0.393 * vapor_pressure + 3.94
        wbgt.attrs = {
            'units': 'deg C',
            'long_name': 'Wet Bulb Globe Temperature (using ACSM '
                'calculation method)',
            'source': 'ACSM, 1984: Prevention of thermal injuries during '
                'distance running: Position stand.'
                'Medical Journal of Australia,  141 (12 – 13), 876 – 879, '
                'doi: 10.5694/j.1326-5377.1984.tb132981.x'
        }
        self.data['wbgt'] = wbgt

class WBGTapprox_BernardCalculator(tcitool.Calculator):
    def __init__(self,tool):
        super().__init__(tool)
        self.export_params = {'wbgt':'wbgt_bernard'}
        self.tool.require_data('t2m','e_kPa','solza')
    def main(self):
        t2mC = self.tool.data.get('t2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['t2m']))
        direct_sun = xr.where(self.tool.data['solza']<1.57079615,1,0)
        direct_sun.attrs = {
            'units': 'bool',
            'long_name': 'Sensor in direct sunlight (True if so)',
        }
        wbgt = (1.1 + 0.66*t2mC +
            2.9*self.tool.data['e_kPa'] +
            direct_sun * -1.8)
        wbgt.attrs = {
            'units': 'deg C',
            'long_name': 'Wet Bulb Globe Temperature (using Bernard & Barrow '
                'calculation method)',
            'source': 'Bernard, T. E., en C. A. Barrow, 2013: Empirical '
                'Approach to Outdoor WBGT from Meteorological Data and '
                'Performance of Two Different Instrument Designs. '
                'Industrial Health, 51 (1), 79–85, '
                'doi: 10.2486/indhealth.2012-0160'
        }
        self.data['wbgt'] = wbgt

class WBGTapprox_DimiceliCalculator(tcitool.Calculator):
    def __init__(self,tool):
        super().__init__(tool)
        self.export_params = {'wbgt':'wbgt_dimiceli'}
        self.tool.require_data('t2m','rh')
    def main(self):
        t2mC = self.tool.data.get('t2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['t2m']))
        rh_procent = self.tool.data.get('rh_procent',
            tcitool.UnitFuncs.rhfraction2procent(self.tool.data['rh']))
        wbgt = (-5.806
                + 0.672*t2mC
                - 0.006*np.power(t2mC,2)
                + 0.061*rh_procent
                + 0.004*rh_procent*t2mC
                + 9.9e-5*rh_procent*np.power(t2mC,2)
                - 3.3e-5*np.power(rh_procent,2)
                - 5e-6*np.power(rh_procent,2)*t2mC
                - 1e-7*np.power(rh_procent,2)*np.power(t2mC,2))
        wbgt.attrs = {
            'units': 'deg C',
            'long_name': 'Wet Bulb Globe Temperature (using Dimiceli et al. '
                'calculation method)',
            'source': 'Dimiceli, V. E., S. F. Piltz, en S. A. Amburn, 2013: '
                'Black Globe Temperature Estimate for the WBGT Index, '
                'IAENG Transactions on Engineering Technologies: Special '
                'Edition of the World Congress on Engineering and Computer '
                'Science 2011, 323–334. Springer, Dordrecht, '
                'doi: 10.1007/978-94-007-4786-9_26.'
        }
        self.data['wbgt'] = wbgt


class WBGTapprox_GommersCalculator(tcitool.Calculator):
    def __init__(self,tool):
        super().__init__(tool)
        self.export_params = {'wbgt':'wbgt_gommers'}
        self.tool.require_data('t2m','d2m','skt','ws10','fsr','Isw_in')
    def main(self):
        t2mC = self.tool.data.get('t2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['t2m']))
        d2mC = self.tool.data.get('d2mC',
            tcitool.UnitFuncs.tempK2C(self.tool.data['d2m']))
        Isw_in = np.clip(self.tool.data['Isw_in'],0,None)
        ACSM_vapor_pressure = 6.112 * np.exp((17.67*d2mC)/(d2mC+243.5))
        ACSM_wbgt = 0.567 * t2mC + 0.393 * ACSM_vapor_pressure + 3.94
        wind2m = (self.tool.data['ws10']
            * (np.log(2 / self.tool.data['fsr'])
                / np.log(10 / self.tool.data['fsr'])))
        wind2m = np.clip(wind2m,0.1,None)
        t2mClog = np.log(t2mC)

        wbgt = xr.where(ACSM_wbgt<0,np.nan,
            xr.where(t2mClog<0.01,np.nan,
            -17.250591
            + 0.4253438 * np.sqrt(self.tool.data['skt']) * np.sqrt(ACSM_wbgt)
            + 7.09012e-05 * np.power(Isw_in, 1.5) * np.reciprocal(wind2m)
            - 4.750152e-06 * np.power(Isw_in, 1.5) * np.reciprocal(wind2m ** 2)
            + 0.04459706 * t2mClog * np.sqrt(Isw_in)
            - 7.534906e-04 * Isw_in * t2mClog))
        wbgt.attrs = {
            'units': 'deg C',
            'long_name': 'WBGT using the Gommers Stepwise Approximation '
                'of Liljegren\'s Argonne model',
            'source': 'Gommers, D. J., 2019: Modelleren van de Wet Bulb Globe '
                'Temperature: Ter voorkoming van hitteletsel. MSc Internship '
                'Report. Wageningen University & Research and the Joint '
                'Meteorological Group, Royal Netherlands Air Force.'
        }
        self.data['wbgt'] = wbgt
