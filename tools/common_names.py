# TEB version 4.0.1 (https://github.com/teb-model/teb).
# Copyright 2018 D. Meyer. Licensed under CeCILL version 2.1.

standard_quantity_names = {
    "Atmospheric pressure": {
        "unit": "\\rm Pa",
        "symbol": "p"
    },
    "Dry-bulb air temperature": {
        "unit": "\\rm K",
        "symbol": "T"
    },
    "Relative humidity": {
        "unit": "\%",
        "symbol": "\\rm RH"
    },
    "Downwelling shortwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "S^{\\downarrow}"
    },
    "Downwelling diffuse shortwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "S^{\\downarrow\\Downarrow}"
    },
    "Downwelling direct shortwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "S^{\\downarrow\\downarrow}"
    },
    "Upwelling shortwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "S^{\\uparrow}"
    },
    "Downwelling longwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "L^{\\downarrow}"
    },
    "Upwelling longwave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "L^{\\uparrow}"
    },
    "Net all-wave radiation flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q^*"
    },
    "Wind speed": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "U"
    },
    "Wind direction": {
        "unit": "\\rm rad",
        "symbol": "U_{\\theta}"
    },
    "CO2 concentration": {
        "unit": "\\rm g\ m^{-3}",
        "symbol": "CO_{2}"
    },
    "Mass mixing ratio of water vapour": {
        "unit": "\\rm g\ kg^{-1}",
        "symbol": "r"
    },
    "Turbulent sensible heat flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_H"
    },
    "Turbulent latent heat flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_E"
    },
    "Rainfall rate": {
        "unit": "\\rm kg\ m^{-2}\ s^{-1}",
        "symbol": "\\rm RR"
    },
    "Snowfall rate": {
        "unit": "\\rm kg\ m^{-2}\ s^{-1}",
        "symbol": "\\rm SR"
    },
    "Buildings\' power demand for cooling": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_{\\rm cooling}"
    },
    "Buildings\' power demand for heating": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_{\\rm heating}"
    },
    "Coefficient of performance": {
        "unit": "1",
        "symbol": "COP"
    },
    "Thermal power production of solar panels on roofs": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_{\\rm thermal}"
    },
    "Electrical power production of solar panels on roofs": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_{\\rm electrical}"
    },
    "Dry-bulb air temperature at half building height": {
        "unit": "\\rm K",
        "symbol": "T_{\\rm canyon}"
    },
    "Mass mixing ratio of water vapour at half building height": {
        "unit": "\\rm kg\ kg^{-1}",
        "symbol": "r_{\\rm canyon}"
    },
    "Canyon horizontal wind speed": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "U_{\\rm canyon}"
    },
    "Canyon horizontal wind direction": {
        "unit": "\\rm rad",
        "symbol": "U_{\\theta_{canyon}}"
    },
    "Zonal component of wind velocity at half building height": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "u_{\\rm canyon}"
    },
    "Meridional component of wind velocity at half building height": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "v_{\\rm canyon}"
    },
    "Zonal component of wind velocity": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "u"
    },
    "Meridional component of wind velocity": {
        "unit": "\\rm m\ s^{-1}",
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
    "Evaporation mass flux density": {
        "unit": "\\rm kg\ m^{-2}\ s^{-1}",
        "symbol": "E"
    },
    "Ground heat flux density": {
        "unit": "\\rm W\ m^{-2}",
        "symbol": "Q_G"
    },
    "Surface (skin) temperature": {
        "unit": "1",
        "symbol": "T_{s}"
    },
    "Specific humidity": {
        "unit": "\\rm kg\ kg^{-1}",
        "symbol": "q"
    },
    "Canyon specific humidity": {
        "unit": "\\rm kg\ kg^{-1}",
        "symbol": "q_{\\rm canyon}"
    },
    "Surface friction velocity": {
        "unit": "\\rm m\ s^{-1}",
        "symbol": "u_{\\star}"
    }
}


teb_quantity_names_inputs = {
    'CO2 concentration': 'CO2',
    'Downwelling direct shortwave radiation flux density': 'DIR_SW',
    'Wind direction': 'DIR',
    'Downwelling longwave radiation flux density': 'LW',
    'Atmospheric pressure': 'PS',
    'Specific humidity': 'QA',
    'Rainfall rate': 'RAIN',
    'Downwelling diffuse shortwave radiation flux density': 'SCA_SW',
    'Snowfall rate': 'SNOW',
    'Dry-bulb air temperature': 'TA',
    'Wind speed': 'WIND'
}

teb_quantity_names_output = {
    'Dry-bulb air temperature at half building height': 'T_CANYON',
    'Canyon horizontal wind speed': 'U_CANYON',
    'Canyon specific humidity': 'Q_CANYON',
    'Canyon horizontal wind direction': 'DIR_CANYON',
    'Buildings\' power demand for cooling': 'HVAC_COOL',
    'Buildings\' power demand for heating': 'HVAC_HEAT',
    'Coefficient of performance': 'COP',
    'Net all-wave radiation flux density': 'RN_TOWN',
    'Turbulent sensible heat flux density': 'H_TOWN',
    'Turbulent latent heat flux density': 'LE_TOWN',
    'Evaporation mass flux density': 'EVAP_TOWN',
    'Ground heat flux density': 'GFLUX_TOWN',
    'Specific humidity': 'Q_TOWN',
    'Surface friction velocity': 'USTAR_TOWN',
    'Surface albedo': 'ALB_TOWN',
    'Surface emissivity': 'EMIS_TOWN',
    'Surface (skin) temperature': 'TS_TOWN',
    'Thermal power production of solar panels on roofs': 'THER_PROD_PANEL',
    'Electrical power production of solar panels on roofs': 'PHOT_PROD_PANEL',
    'Canyon atmospheric pressure': 'P_CANYON',
    'Buildings\' internal dry-bulb air temperature': 'TI_BLD',
    'Road surface (skin) temperature': 'T_ROAD1',
    'Roof surface (skin) temperature': 'T_ROOF1',
    'Wall A surface (skin) temperature': 'T_WALLA1',
    'Wall B surface (skin) temperature': 'T_WALLB1',
}


teb_quantity_names_output_by_keys = {v:k for k,v in teb_quantity_names_output.items()}