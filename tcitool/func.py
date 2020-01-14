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
    STEFAN_BOLTZMANN = 5.67e-8

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
    @classmethod
    def wind_at_height_using_fsr(cls,ws10,fsr,height):
        """Calculates the reduced wind speed, at the surface, using the surface roughness

        Args:
            ws10: Windspeed at 10m [m/s]
            fsr: Surface roughness [m]
            height: Height at which the windspeed will be calculated [m]

        Returns:
            Reduced windspeed at height_target [m/s]
        """
        return ws10*(np.log(height/fsr)/np.log(10/fsr))

class ThermodynamicFuncs(object):
    GAS_CONSTANT = 8314.34
    GAS_CONSTANT_AIR = 287.05
    GRAVITATIONAL_AC = 9.81
    MOLAR_MASS_AIR = 28.965
    MOLAR_MASS_H2O = 18.015

    @classmethod
    def atmospheric_emissivity(cls,e_kPa): # ε_atm [-]
        """Calculates the atmospheric emissivity [-]

        Args:
            e_kPa: specific humidity [kPa]

        Returns:
            atmospheric emissivity [-], $\varepsilon_{atm}$

        (Liljegren et al., 2008; Oke T. R., 1978)
        """
        return 0.575 * np.power(e_kPa,0.143)
    @classmethod
    def density(cls,temp_K,P_kPa):
        """Calculates the densisty of the air [kg/m³]

        Args:
            temp_K: Air temperature [K]
            P_kPa: Air pressure [kPa]

        Returns:
            density, $\rho$ [kg/m³]
        """
        return (P_kPa*1e3)/((cls.GAS_CONSTANT/cls.MOLAR_MASS_AIR)*temp_K)
    @classmethod
    def dynamic_viscosity(cls,temp_K):
        """Calculates the dynamic viscosity [Pa s]

        Args:
            temp_K: Air temperature [K]

        Returns:
            dynamic viscosity [Pa s], $\mu$
        """
        return 1.458e-6*(temp_K**1.5) / (temp_K +110.4)  # Pa s (or Ns/m²)
    @classmethod
    def heat_of_evaporation(cls,temp_K):
        """Calculates the heat of evaporation for water

        Args:
            temp_K: Air temperature [K]

        Returns:
            heat of evaporation for water [J/kg], $\Delta H_v$ or $L_v$
        """
        return (313.15-temp_K)/30*-71100+2.4073e6
    @classmethod
    def kinematic_viscosity(cls,temp_K,P_kPa):
        """Calculates the kinematic viscosity [m²/s]

        Args:
            temp_K: Air temperature [K]
            P_kPa: Air pressure [kPa]

        Returns:
            kinematic viscosity [m²/s], $\nu$

        """
        return cls.dynamic_viscosity(temp_K)/cls.density(temp_K,P_kPa)
    @classmethod
    def thermal_capacity(cls,temp_K):
        """Calculates the thermal capacity of the air [J/(kg K)]

        Args:
            temp_K: Air temperature [K]

        Returns:
            thermal capacity of the air, at constant pressure [J/(kg K)], $C_p$
        """
        return 1002.5+275e-6 * np.power(temp_K-200,2)
    @classmethod
    def thermal_conductivity(cls,temp_K):
        """Calculates the thermal conductivity of the air [W/(m K)]

        Args:
            temp_K: Air temperature [K]

        Returns:
            thermal conductivity of the air [W/(m K)], $k$ or $\lambda$
        """
        return 0.02624 * np.power(temp_K/300,0.8646)
    @classmethod
    def thermal_expansion(cls,temp_K):
        """Calculates the thermal expansion of a ideal gas [K^-1]

        Args:
            temp_K: Air temperature [K]

        Returns:
            thermal conductivity of the air [K^-1], $\alpha_V$
        """
         # α_V [K^-1]
        return 1/temp_K
    @classmethod
    def thermal_diffusivity(cls,temp_K,P_kPa):
        """Calculates the thermal diffusivity [m²/s]

        Args:
            temp_K: Air temperature [K]
            P_kPa: Air pressure [kPa]

        Returns:
            thermal diffusivity [m²/s], $\alpha$, $a$, $h$, $\kappa$, $L$ or $D$

        """
        return cls.thermal_conductivity(temp_K) / cls.volumetric_heat_capacity(temp_K,P_kPa)
    @classmethod
    def volumetric_heat_capacity(cls,temp_K,P_kPa):
        """Calculates the volumetric heat capacity [J/(m³ K)]

        Args:
            temp_K: Air temperature [K]
            P_kPa: Air pressure [kPa]

        Returns:
            volumetric heat capacity [J/(m³ K)], $s$

        """
        return cls.density(temp_K,P_kPa) * cls.thermal_capacity(temp_K)

    @classmethod
    def reynolds_number(cls,temp_K,P_kPa,wind_ms,length_m):
        """Calculates the Reynolds number [-]

        Args:
            temp_K: Air temperature [K]
            P_kPa: Air pressure [kPa]
            wind_ms: Wind speed [m/s]
            length_m: Characteristic length [m]

        Returns:
            Reynolds number [-], Re
        """
        return wind_ms*length_m/cls.kinematic_viscosity(temp_K,P_kPa)
    @classmethod
    def prandtl_number(cls,temp_K):
        """Calculates the Prandtl number [-]

        Args:
            temp_K: Air temperature [K]

        Returns:
            Prandtl number [-], Pr
        """
        return cls.thermal_capacity(temp_K) \
            * cls.dynamic_viscosity(temp_K) \
            / cls.thermal_conductivity(temp_K)
    @classmethod
    def schmidt_number(cls,temp_K,P_kPa):
        """Calculates the Schmidt number [-]

            Args:
                temp_K: Air temperature [K]
                P_kPa: Air pressure [kPa]

            Returns:
                Schmidt number [-], Sc
        """
        return cls.dynamic_viscosity(temp_K) / ( \
                    cls.density(temp_K,P_kPa) \
                    * cls.thermal_diffusivity(temp_K,P_kPa) \
            )

u = UnitFuncs
m = MeteoFuncs
td = ThermodynamicFuncs
