# TCItool
A Python3 package to calculate thermal comfort indexes (such as the Wet Bulb Globe Temperature), from ECMWF-data (e.g. ERA5 or IFS).

![Python v3.7](https://img.shields.io/badge/python-v3.7-blue)
![License](https://img.shields.io/github/license/papagaai35/TCItool)
![Last commit](https://img.shields.io/github/last-commit/papagaai35/TCItool)
![Release v0.1](https://img.shields.io/badge/release-v0.1-blue)

## Supported indexes
This package is now able to calculate the following indexes
* Wind-Chill Equivalent Temperature
	* Using the JAG/TI method ([Osczevski & Bluestein 2005](https://doi.org/10.1175/BAMS-86-10-1453))
* Wet Bulb Globe Temperature
	* Using the improved Argonne model (adapted from [Liljegren et al. 2008](https://doi.org/10.1080/15459620802310770)). This method uses less approximations to solve the radiation balance as the Argonne model, using the skin temperature and albedo. These were not available in the data of Liljegren et al. 2008, but are available in the ECMWF datasets.
	* Using the approximations by [ACSM (1984)](https://doi.org/10.5694/j.1326-5377.1984.tb132981.x), [Bernard & Barrow (2013)](https://doi.org/10.2486/indhealth.2012-0160) and [Dimiceli et al. (2013)](https://doi.org/10.1007/978-94-007-4786-9_26)

## Using this package
`Tool` serves as the main entry point for this module. You can create a tool by calling
```
import tcitool
tool = tcitool.Tool()
```

The tool can than be used to load data, calculate indexes, and export data.
```
tool.data.load('./ECMWF_ERA5.nc')
tool.options.update({'radiation_cumulative': False,'radiation_integration_time': 3600})
tool.calculate('windchill_jagti')
```

All data can be accessed using the `tool.data` interface. This interface is based on `xarray` with a few extra functions.
The xarray.Dataset can be accessed using `tool.data.ds`.  
```
tool.data['wcet_jagit']
>>> <xarray.DataArray 'wcet_jagti' (time: 3, latitude: 3, longitude: 3)>
	array([[[16.  , 17.85, 18.3 ],
                [17.83, 14.74, 12.9 ],
                [14.44, 14.36, 12.19]],
               [[16.03, 17.77, 18.05],
                [17.61, 15.21, 15.27],
                [14.39, 14.91, 13.94]],
               [[19.36, 21.33, 20.75],
                [21.22, 25.93, 27.53],
                [26.35, 27.14, 26.68]]])
tool.data['wcet_jagti'].attrs
>>> {'units': 'deg C',
	 'long_name': 'Wind Chill Equivalent Temperature (using JAG/TI calculation method)',
	 'source': 'Eq. 2.3 in https://cdn.knmi.nl/system/downloads/files/000/000/016/original/gevoelstemperatuur.pdf?1433939065'}
```

The Tool interface also provides exporting capabilities.
```
tool.data.save('./ECMWF_ERA5_withWC.nc')
```

For more information, view the documentation using `help(tool)` (or `help(tool.data)` for more info about the data object, for example).


## Development Outline
When diving in the source code of this package, it is usefull to keep the thing below in mind. This is the general setup of the package.
* `tool` provides the main interactions with the user.
* The `func` module is used to  preform simple operations (e.g. unit conversions, density, Reynold number), using only the data submitted to the function. These are used by the Generators and Calculators. Methods in this module do not run on the enitre dataset, but only work with the data given. Classes in this module only contain `classmethod`-s and are a way to namespace the different operations.
* `Generators` are classes that preform a operation to the data set, uppon request. For example, the user, or a `Calculator` may require the wind speed, while this quantity is unavailable in `tool.data`. What is available, is the u and v components of the wind. In that case, the relevant Generator-method will start and use the `tcitool.func.MeteoFuncs.wind_speed`-function to convert between these parameters, and store the result back into `tool.data`.
* `Calculators` are classes that provide the core functionality of this module. These are often longer running processes that will calculate many different parameters to be able to calculate the requested variable. Using `tool.require_data()` it will kick-off the nessesary Generator-methods. It may also use many different `func`-methods in its internals.
* The `data` module is a wrapper for a `xarray.Dataset`, and is used to provide some handy functions to the Dataset.

`Calculator`-classes need to be registerd into the `tcitool.tool.Tool.calculators`-dictionary.  
`Generator`-methods also need to be registerd into the `tcitool.gens.registery`-class to be able to automaticaly start, when a certain dataset is nessesary. `Generators` may not invoke `tool.require_data`.
