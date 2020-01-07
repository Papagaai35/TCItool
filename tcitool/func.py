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

class MeteoFuncs(object):
    @classmethod
    def wind_speed(cls,u,v):
        """Caluclates the wind speed from U and V components"""
        return np.sqrt(np.power(u,2)+np.power(v,2))
    @classmethod
    def wind_direction(cls,u,v):
        """Caluclates the wind direction from U and V components"""
        return np.deg2rad(180+180/np.pi*np.arctan2(u,v))
