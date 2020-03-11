# TEB version 4.0.0 (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL version 2.1.

standard_quantity_names = {
    "Atmospheric pressure": {
        "unit": "hPa",
        "symbol": "p"
    },
    "Dry-bulb air temperature": {
        "unit": "K",
        "symbol": "T"
    },
    "Relative humidity": {
        "unit": "%",
        "symbol": "\\rm RH"
    },
    "Total sky shortwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "S⇊"
    },
    "Sky diffuse shortwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "S⇓"
    },
    "Sky direct shortwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "S↓"
    },
    "Surface reflected shortwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "S↑"
    },
    "Sky longwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "L↓"
    },
    "Surface longwave radiation flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "L↑"
    },
    "Surface net radiative flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "Q^*"
    },
    "Wind speed": {
        "unit": "m\ s^{-1}",
        "symbol": "U"
    },
    "Wind direction": {
        "unit": "rad",
        "symbol": "U_{\\theta}"
    },
    "CO2 concentration": {
        "unit": "g\ m^{-3}",
        "symbol": "CO_{2}"
    },
    "Mass mixing ratio of water vapour": {
        "unit": "g\ kg^{-1}",
        "symbol": "r"
    },
    "Surface turbulent sensible heat flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "H"
    },
    "Surface turbulent latent heat flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "LE"
    },
    "Rain rate": {
        "unit": "mm\ h^{-1}",
        "symbol": "RR"
    },
    "Buildings\' energy demand for cooling": {
        "unit": "W\ m^{-2}",
        "symbol": "E_{\\rm cooling}"
    },
    "Buildings\' energy demand for heating": {
        "unit": "W\ m^{-2}",
        "symbol": "E_{\\rm heating}"
    },
    "Thermal energy production of solar panels on roofs": {
        "unit": "W\ m^{-2}",
        "symbol": "Q_{\\rm solar}"
    },
    "Electrical energy production of solar panels on roofs": {
        "unit": "W\ m^{-2}",
        "symbol": "E_{\\rm solar}"
    },
    "Dry-bulb air temperature at half building height": {
        "unit": "K",
        "symbol": "T_{\\rm canyon}"
    },
    "Mass mixing ratio of water vapour at half building height": {
        "unit": "kg\ kg^{-1}",
        "symbol": "r_{\\rm canyon}"
    },
    "Canyon horizontal wind speed": {
        "unit": "m\ s^{-1}",
        "symbol": "U_{\\rm canyon}"
    },
    "Canyon horizontal wind direction": {
        "unit": "rad",
        "symbol": "U_{\\theta_{canyon}}"
    },
    "Zonal component of wind velocity at half building height": {
        "unit": "m\ s^{-1}",
        "symbol": "u_{\\rm canyon}"
    },
    "Meridional component of wind velocity at half building height": {
        "unit": "m\ s^{-1}",
        "symbol": "v_{\\rm canyon}"
    },
    "Zonal component of wind velocity": {
        "unit": "m\ s^{-1}",
        "symbol": "u"
    },
    "Meridional component of wind velocity": {
        "unit": "m\ s^{-1}",
        "symbol": "v"
    },
    "Surface emissivity": {
        "unit": "1",
        "symbol": "\epsilon"
    },
    "Surface albedo": {
        "unit": "1",
        "symbol": "\\alpha"
    },
    "Surface evaporative moisture flux density": {
        "unit": "kg\ m^{-2}\ s^{-1}",
        "symbol": "E"
    },
    "Surface ground heat flux density": {
        "unit": "W\ m^{-2}",
        "symbol": "G"
    },
    "Surface skin temperature": {
        "unit": "1",
        "symbol": "T_{s}"
    },
    "Specific humidity": {
        "unit": "kg\ kg^{-1}",
        "symbol": "q"
    },
    "Canyon specific humidity": {
        "unit": "kg\ kg^{-1}",
        "symbol": "q_{\\rm canyon}"
    },
    "Surface friction velocity": {
        "unit": "m\ s^{-1}",
        "symbol": "u_{\\star}"
    }
}


teb_quantity_names_inputs = {
    'CO2 concentration': 'CO2',
    'Sky direct shortwave radiation flux density': 'DIR_SW',
    'Wind direction': 'DIR',
    'Sky longwave radiation flux density': 'LW',
    'Atmospheric pressure': 'PS',
    'Specific humidity': 'QA',
    'Rain rate': 'RAIN',
    'Sky diffuse shortwave radiation flux density': 'SCA_SW',
    'Snow rate': 'SNOW',
    'Dry-bulb air temperature': 'TA',
    'Wind speed': 'WIND'
}

teb_quantity_names_output = {
    'Dry-bulb air temperature at half building height': 'T_CANYON',
    'Canyon horizontal wind speed': 'U_CANYON',
    'Canyon specific humidity': 'Q_CANYON',
    'Canyon horizontal wind direction': 'DIR_CANYON',
    'Buildings\' energy demand for cooling': 'HVAC_COOL',
    'Buildings\' energy demand for heating': 'HVAC_HEAT',
    'Surface net radiative flux density': 'RN_TOWN',
    'Surface turbulent sensible heat flux density': 'H_TOWN',
    'Surface turbulent latent heat flux density': 'LE_TOWN',
    'Surface evaporative moisture flux density': 'EVAP_TOWN',
    'Surface ground heat flux density': 'GFLUX_TOWN',
    'Specific humidity': 'Q_TOWN',
    'Surface friction velocity': 'USTAR_TOWN',
    'Surface albedo': 'ALB_TOWN',
    'Surface emissivity': 'EMIS_TOWN',
    'Surface skin temperature': 'TS_TOWN',
    'Thermal energy production of solar panels on roofs': 'THER_PROD_PANEL',
    'Electrical energy production of solar panels on roofs': 'PHOT_PROD_PANEL',
}
