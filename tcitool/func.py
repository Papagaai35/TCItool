import numpy as np

class UnitFuncs(object):
    @classmethod
    def dms2deg(cls,d,m,s):
        """Converts a degree,minute,seconds to decimal degrees"""
        return d+m/60+s/3600
    @classmethod
    def dms2rad(cls,d,m,s):
        """Converts a degree,minute,seconds to radians"""
        return np.deg2rad(cls.dms2deg(d,m,s))
    @classmethod
    def tempC2K(cls,temp_C):
        """Converts a temperature in Celcius to Kelvin"""
        return temp_C + 273.15
    @classmethod
    def tempK2C(cls,temp_K):
        """Converts a temperature in Kelvin to Celcius"""
        return temp_K - 273.15
    @classmethod
    def presPa2hPa(cls,p_Pa):
        return p_Pa/100
    @classmethod
    def presPa2kPa(cls,p_Pa):
        return p_Pa/1000
    @classmethod
    def preshPa2Pa(cls,p_hPa):
        return p_hPa*100
    @classmethod
    def preshPa2kPa(cls,p_hPa):
        return p_hPa/10
    @classmethod
    def preskPa2Pa(cls,p_kPa):
        return p_kPa*1000
    @classmethod
    def preskPa2hPa(cls,p_kPa):
        return p_kPa*100
    @classmethod
    def rhfraction2procent(cls,rh):
        return rh*100
    @classmethod
    def rhprocent2fraction(cls,rh_procent):
        return rh_procent/100

class MeteoFuncs(object):
    @classmethod
    def dewpoint(cls,e_kPa):
        """Calculates the dew point, based on the specific humidity

        Args:
            e_kPa: specific humidity [kPa]

        Returns:
            dewpoint [K], $T_{dew}$
        """
        x = np.log(e_kPa/0.611)/17.2694
        return (-35.86*x+273.16)/(1-x)
    @classmethod
    def saturated_vapor_pressure(cls,temp_K):
        """Calculates the saturated vapor pressure [kPa]

        Args:
            temp_K: Air temperature [K]

        Returns:
            saturated vapor pressure [kPa], $e_{sat}$
        """
        return 0.611 * np.exp(17.2694*(temp_K-273.16)/(temp_K-35.86))
    @classmethod
    def wind_speed(cls,u,v):
        """Caluclates the wind speed from U and V components"""
        return np.sqrt(np.power(u,2)+np.power(v,2))
    @classmethod
    def wind_direction(cls,u,v):
        """Caluclates the wind direction from U and V components"""
        return np.deg2rad(180+180/np.pi*np.arctan2(u,v))
