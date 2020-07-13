import numpy as np
import xarray as xr

import tcitool

class SolarGenerators(object):
    @classmethod
    def register_generators(cls,gr):
        gr.register(cls.main,
            ['soldist','solhour','solza','solazimuth'],
            ['time','longitude','latitude'])

    @classmethod
    def main(cls,tool):
        """Calculator for many different solar parameters.

        Source for the calculation:
            https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
        """
        ts, lon, lat = cls.extractCoordVars(tool.data.ds)
        solarparam = cls.solarParamNOAA(ts,lon,lat)
        for key, array in solarparam.items():
            tool.data[key] = array
        tool.data.transpose_default()

    @classmethod
    def extractCoordVars(cls,ds):
        def get_coord_var(ds,coord):
            coordxarray_1d = xr.DataArray(np.atleast_1d(ds.coords[coord].values),
                dims=[coord],attrs=ds.coords[coord].attrs)
            _, coordxarray_md = xr.broadcast(ds,coordxarray_1d)
            return coordxarray_md
        ts = (ds['ts'] if 'ts' in ds.data_vars else get_coord_var(ds,'time'))
        lon = (ds['lon'] if 'lon' in ds.data_vars else get_coord_var(ds,'longitude'))
        lat = (ds['lat'] if 'lat' in ds.data_vars else get_coord_var(ds,'latitude'))
        return ts, lon, lat

    @classmethod
    def solarParamNOAA(cls,ts,lon,lat):
        lon_rad = np.deg2rad(lon)
        lon_rad.attrs = {'units':'rad','long_name': 'longitude'}
        lat_rad = np.deg2rad(lat)
        lat_rad.attrs = {'units':'rad','long_name': 'latitude'}
        unix_ts = ts.astype('datetime64[s]').astype('int')
        unix_ts.attrs = {'units':'s','long_name':'seconds since 1970-01-01'}
        julian_day = unix_ts/86400.0 + 2440587.5
        julian_day.attrs = {'units':'days',
            'long_name':'days since 4713BC-01-01 12h UTC'}
        julian_century = (julian_day - 2451545) / 36525
        julian_century.attrs = {'units':'centuries',
            'long_name':'decimal julian centuries since 2000-01-01 12hUTC'}
        time_frac = unix_ts / 86400 % 1
        time_frac.attrs = {'units':'days',
            'long_name':'time of a day represented as decimal day [0-1]'}
        long_sun_deg = ( 280.46646 + julian_century*(
            36000.76983 + 0.0003032*julian_century )) % 360
        long_sun_deg.attrs = {'units':'deg','long_name':'Geom Mean Long Sun'}
        long_sun = np.deg2rad(long_sun_deg)
        long_sun.attrs = {'units':'rad','long_name':'Geom Mean Long Sun'}
        anom_sun_deg = 357.52911 + julian_century*(
            35999.05029 - 0.0001537*julian_century)
        anom_sun_deg.attrs = {'units':'deg','long_name':'Geom Mean Anom Sun'}
        anom_sun = np.deg2rad(anom_sun_deg)
        anom_sun.attrs = {'units':'rad','long_name':'Geom Mean Anom Sun'}
        eccent_earth_orbit = 0.016708634-julian_century*(
            0.000042037+0.0000001267*julian_century)
        eccent_earth_orbit.attrs = {'units':'','long_name':'Eccent Earth Orbit'}
        sun_eqofctr_deg = (
            np.sin(anom_sun) * (1.914602 -
                julian_century * (0.004817 + 0.000014 * julian_century)) +
            np.sin(2 * anom_sun) * (0.019993 - 0.000101 * julian_century) +
            np.sin(3 * anom_sun) * 0.000289
        )
        sun_eqofctr_deg.attrs = {'units':'deg','long_name':'Sun Eq of Ctr'}
        sun_eqofctr = np.deg2rad(sun_eqofctr_deg)
        sun_eqofctr.attrs = {'units':'rad','long_name':'Sun Eq of Ctr'}
        long_sun_true = long_sun + sun_eqofctr
        long_sun_true.attrs = {'units':'rad','long_name':'Sun True Long'}
        anom_sun_true = anom_sun + sun_eqofctr
        anom_sun_true.attrs = {'units':'rad','long_name':'Sun True Anom'}
        sun_rad_vector = (
            ( 1.000001018 * (1 - np.power(eccent_earth_orbit,2))) /
            (1+eccent_earth_orbit*np.cos(anom_sun_true))
        )
        sun_rad_vector.attrs = {'units':'au','long_name':'Sun Rad Vector'}
        sun_app_long = np.deg2rad(
            np.rad2deg(long_sun_true) - 0.00569 - 0.00478 *
            np.sin(np.deg2rad(125.04 - 1934.136 * julian_century)) )
        sun_app_long.attrs = {'units':'rad','long_name':'Sun App Long'}
        mean_obliq = np.deg2rad(23 + (26 + (21.448 - julian_century * (
                46.815 + julian_century * (0.00059 - julian_century * 0.001813)
            )) / 60) / 60)
        mean_obliq.attrs = {'units':'rad','long_name':'Mean Obliq Ecliptic'}
        obliq_corr = np.deg2rad(
            np.rad2deg(mean_obliq) +
            0.00256 * np.cos(np.deg2rad(125.04 - 1934.136 * julian_century)) )
        obliq_corr.attrs = {'units':'rad','long_name':'Obliq Corr'}

        sun_rt_ascen = np.arctan2(
            np.cos(obliq_corr)*np.sin(sun_app_long),np.cos(sun_app_long))
        sun_rt_ascen.attrs = {'units':'rad','long_name':'Sun Rt Ascen'}
        sun_declin = np.arcsin(np.sin(obliq_corr)*np.sin(sun_app_long))
        sun_declin.attrs = {'units':'rad','long_name':'Sun Declin'}
        var_y = np.tan(obliq_corr/2)**2
        eq_of_time = 4 * np.rad2deg(
            var_y * np.sin(2*long_sun)
            - 2 * eccent_earth_orbit * np.sin(anom_sun)
            + 4 * eccent_earth_orbit * var_y * np.sin(anom_sun) * np.cos(2*long_sun)
            - 0.5 * var_y**2 *np.sin(4*long_sun)
            - 1.25 * eccent_earth_orbit**2 * np.sin(2*anom_sun)
        )
        eq_of_time.attrs = {'units':'min','long_name':'Eq of Time'}
        true_solar_time = (
            time_frac * 1440 + eq_of_time + 4 * np.rad2deg(lon_rad) ) % 1440
        true_solar_time.attrs = {'units':'min','long_name':'True Solar Time'}
        hour_angle = np.deg2rad(true_solar_time/4 - 180)
        hour_angle.attrs = {'units':'rad','long_name':'Hour Angle'}

        zenith_uncorr = np.arccos(
            np.sin(lat_rad) * np.sin(sun_declin)
             + np.cos(lat_rad) * np.cos(sun_declin)*np.cos(hour_angle) )
        zenith_uncorr.attrs = {'units':'rad','long_name':'Solar Zenith Angle'}
        elevation = np.deg2rad(90 - np.rad2deg(zenith_uncorr))
        elevation.attrs = {'units':'rad',
            'long_name':'Solar Elevation Angle'}

        atmospheric_refraction = np.deg2rad(
            xr.where(elevation>np.deg2rad(85),0,
            xr.where(elevation>np.deg2rad(5),(
                58.1 / np.tan(elevation)
                - 0.07 / np.power(np.tan(elevation),3)
                + 8.6e-5 / np.power(np.tan(elevation),5) ),
            xr.where(elevation>np.deg2rad(-0.575),(
                1735 + elevation * (
                    -518.2 + elevation * (
                        103.4 + elevation * (
                        -12.79 + elevation*0.711)))),
                -20.772/np.tan(elevation)
            ))) / 3600)
        atmospheric_refraction.attrs = {'units':'rad',
            'long_name':'Approx Atmospheric Refraction'}
        elevation_corr = elevation + atmospheric_refraction
        elevation_corr.attrs = {'units':'rad',
            'long_name':'Solar Elevation corrected for atm refraction'}
        zenith = np.deg2rad(90)-elevation_corr
        zenith.attrs = {'units':'rad',
            'long_name':'Solar Zenith Angle corrected for atm refraction'}

        azimuth_arccos = np.arccos( ( (
                (np.sin(lat_rad) * np.cos(zenith_uncorr) - np.sin(sun_declin))
                / (np.cos(lat_rad) * np.sin(zenith_uncorr))
                + 1)
            % 2) - 1)
        azimuth = (xr.where(
            hour_angle>0,
            np.rad2deg(azimuth_arccos)+180,
            540 - np.rad2deg(azimuth_arccos)
        ) % 360)
        azimuth.attrs = {'units':'deg CW from N',
            'long_name':'Solar Azimuth Angle'}
        return xr.Dataset({
            'soldist': sun_rad_vector,
            'solhour': hour_angle,
            'solza': zenith,
            'solazimuth': azimuth,
        }).squeeze(drop=True)
