# Makefile created by mkmf.pl $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include gfortran_args


.DEFAULT:
	-touch $@
all: tester.exe
abor1_sfx.o: src_driver/abor1_sfx.F90 close_file.o modd_surf_conf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/abor1_sfx.F90
add_forecast_to_date_surf.o: src_driver/add_forecast_to_date_surf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/add_forecast_to_date_surf.F90
avg_urban_fluxes.o: src_teb/avg_urban_fluxes.F90 modd_csts.o mode_thermos.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/avg_urban_fluxes.F90
bem.o: src_teb/bem.F90 modd_csts.o mode_thermos.o mode_psychro.o modi_dx_air_cooling_coil_cv.o modi_floor_layer_e_budget.o modi_mass_layer_e_budget.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/bem.F90
bem_morpho.o: src_teb/bem_morpho.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/bem_morpho.F90
bld_e_budget.o: src_teb/bld_e_budget.F90 modd_csts.o modd_surf_par.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/bld_e_budget.F90
circumsolar_rad.o: src_solar/circumsolar_rad.F90 hook.o modd_csts.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_solar/circumsolar_rad.F90
close_file.o: src_driver/close_file.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/close_file.F90
close_file_asc.o: src_driver/close_file_asc.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/close_file_asc.F90
dx_air_cooling_coil_cv.o: src_teb/dx_air_cooling_coil_cv.F90 mode_thermos.o mode_psychro.o modd_csts.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/dx_air_cooling_coil_cv.F90
facade_e_budget.o: src_teb/facade_e_budget.F90 modd_surf_par.o modd_csts.o modi_wall_layer_e_budget.o modi_window_e_budget.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/facade_e_budget.F90
floor_layer_e_budget.o: src_teb/floor_layer_e_budget.F90 modi_layer_e_budget_get_coef.o modi_layer_e_budget.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/floor_layer_e_budget.F90
flxsurf3bx.o: src_teb/flxsurf3bx.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/flxsurf3bx.F
garden.o: src_proxi_SVAT/garden.F90 modd_csts.o mode_thermos.o modd_type_date_surf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/garden.F90
garden_properties.o: src_proxi_SVAT/garden_properties.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/garden_properties.F90
greenroof.o: src_proxi_SVAT/greenroof.F90 modd_csts.o mode_thermos.o modd_type_date_surf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/greenroof.F90
greenroof_properties.o: src_proxi_SVAT/greenroof_properties.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/greenroof_properties.F90
hook.o: src_teb/hook.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/hook.F90
ini_csts.o: src_teb/ini_csts.F90 modd_csts.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/ini_csts.F90
init_surfconsphy.o: src_teb/init_surfconsphy.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/init_surfconsphy.F
layer_e_budget.o: src_teb/layer_e_budget.F90 modi_tridiag_ground.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/layer_e_budget.F90
layer_e_budget_get_coef.o: src_teb/layer_e_budget_get_coef.F90 hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/layer_e_budget_get_coef.F90
mass_layer_e_budget.o: src_teb/mass_layer_e_budget.F90 modi_layer_e_budget_get_coef.o modi_layer_e_budget.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/mass_layer_e_budget.F90
modd_arch.o: src_driver/modd_arch.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/modd_arch.F90
modd_bem_cst.o: src_teb/modd_bem_cst.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_bem_cst.F90
modd_csts.o: src_teb/modd_csts.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_csts.F90
modd_flood_par.o: src_teb/modd_flood_par.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_flood_par.F90
modd_forc_atm.o: src_driver/modd_forc_atm.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/modd_forc_atm.F90
modd_snow_par.o: src_teb/modd_snow_par.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_snow_par.F90
modd_surf_atm.o: src_teb/modd_surf_atm.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_surf_atm.F90
modd_surf_conf.o: src_driver/modd_surf_conf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/modd_surf_conf.F90
modd_surf_par.o: src_teb/modd_surf_par.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_surf_par.F90
modd_type_date_surf.o: src_teb/modd_type_date_surf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_type_date_surf.F90
modd_water_par.o: src_teb/modd_water_par.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modd_water_par.F90
mode_char2real.o: src_driver/mode_char2real.F90 modd_arch.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/mode_char2real.F90
mode_conv_DOE.o: src_teb/mode_conv_DOE.F90 hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/mode_conv_DOE.F90
mode_psychro.o: src_teb/mode_psychro.F90 hook.o modd_csts.o modd_surf_par.o mode_thermos.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/mode_psychro.F90
mode_surf_snow_frac.o: src_teb/mode_surf_snow_frac.F90 hook.o modd_snow_par.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/mode_surf_snow_frac.F90
mode_thermos.o: src_teb/mode_thermos.F90 hook.o modd_csts.o modd_surf_par.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/mode_thermos.F90
modi_avg_urban_fluxes.o: src_teb/modi_avg_urban_fluxes.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_avg_urban_fluxes.F90
modi_bem.o: src_teb/modi_bem.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_bem.F90
modi_bem_morpho.o: src_teb/modi_bem_morpho.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_bem_morpho.F90
modi_bld_e_budget.o: src_teb/modi_bld_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_bld_e_budget.F90
modi_dx_air_cooling_coil_cv.o: src_teb/modi_dx_air_cooling_coil_cv.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_dx_air_cooling_coil_cv.F90
modi_facade_e_budget.o: src_teb/modi_facade_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_facade_e_budget.F90
modi_floor_layer_e_budget.o: src_teb/modi_floor_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_floor_layer_e_budget.F90
modi_flxsurf3bx.o: src_teb/modi_flxsurf3bx.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_flxsurf3bx.f
modi_garden.o: src_proxi_SVAT/modi_garden.F90 modd_type_date_surf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/modi_garden.F90
modi_garden_properties.o: src_proxi_SVAT/modi_garden_properties.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/modi_garden_properties.F90
modi_greenroof.o: src_proxi_SVAT/modi_greenroof.F90 modd_type_date_surf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/modi_greenroof.F90
modi_greenroof_properties.o: src_proxi_SVAT/modi_greenroof_properties.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_proxi_SVAT/modi_greenroof_properties.F90
modi_ini_csts.o: src_teb/modi_ini_csts.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_ini_csts.F90
modi_init_surfconsphy.o: src_teb/modi_init_surfconsphy.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_init_surfconsphy.f
modi_layer_e_budget.o: src_teb/modi_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_layer_e_budget.F90
modi_layer_e_budget_get_coef.o: src_teb/modi_layer_e_budget_get_coef.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_layer_e_budget_get_coef.F90
modi_mass_layer_e_budget.o: src_teb/modi_mass_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_mass_layer_e_budget.F90
modi_road_layer_e_budget.o: src_teb/modi_road_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_road_layer_e_budget.F90
modi_roof_layer_e_budget.o: src_teb/modi_roof_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_roof_layer_e_budget.F90
modi_snow_cover_1layer.o: src_teb/modi_snow_cover_1layer.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_snow_cover_1layer.F90
modi_surface_aero_cond.o: src_teb/modi_surface_aero_cond.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_surface_aero_cond.F90
modi_surface_cd.o: src_teb/modi_surface_cd.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_surface_cd.F90
modi_surface_ri.o: src_teb/modi_surface_ri.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_surface_ri.F90
modi_teb.o: src_teb/modi_teb.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_teb.F90
modi_teb_garden.o: src_teb/modi_teb_garden.F90 modd_type_date_surf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_teb_garden.F90
modi_tridiag_ground.o: src_teb/modi_tridiag_ground.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_tridiag_ground.F90
modi_urban_drag.o: src_teb/modi_urban_drag.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_drag.F90
modi_urban_exch_coef.o: src_teb/modi_urban_exch_coef.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_exch_coef.F90
modi_urban_fluxes.o: src_teb/modi_urban_fluxes.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_fluxes.F90
modi_urban_hydro.o: src_teb/modi_urban_hydro.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_hydro.F90
modi_urban_lw_coef.o: src_teb/modi_urban_lw_coef.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_lw_coef.F90
modi_urban_snow_evol.o: src_teb/modi_urban_snow_evol.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_snow_evol.F90
modi_urban_solar_abs.o: src_teb/modi_urban_solar_abs.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_urban_solar_abs.F90
modi_wall_layer_e_budget.o: src_teb/modi_wall_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_wall_layer_e_budget.F90
modi_wind_threshold.o: src_teb/modi_wind_threshold.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_wind_threshold.F90
modi_window_data.o: src_teb/modi_window_data.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_window_data.F90
modi_window_e_budget.o: src_teb/modi_window_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_window_e_budget.F90
modi_window_shading.o: src_teb/modi_window_shading.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_window_shading.F90
modi_window_shading_availability.o: src_teb/modi_window_shading_availability.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/modi_window_shading_availability.F90
ol_alloc_atm.o: src_driver/ol_alloc_atm.F90 modd_surf_par.o modd_forc_atm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/ol_alloc_atm.F90
ol_read_atm.o: src_driver/ol_read_atm.F90 ol_read_atm_ascii.o mode_thermos.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/ol_read_atm.F90
ol_read_atm_ascii.o: src_driver/ol_read_atm_ascii.F90 read_surf_atm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/ol_read_atm_ascii.F90
ol_time_interp_atm.o: src_driver/ol_time_interp_atm.F90 modd_csts.o modd_surf_par.o modd_forc_atm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/ol_time_interp_atm.F90
open_close_bin_asc_forc.o: src_driver/open_close_bin_asc_forc.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/open_close_bin_asc_forc.F90
read_surf_atm.o: src_driver/read_surf_atm.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/read_surf_atm.F90
road_layer_e_budget.o: src_teb/road_layer_e_budget.F90 modd_csts.o mode_thermos.o modi_layer_e_budget.o modi_layer_e_budget_get_coef.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/road_layer_e_budget.F90
road_wall_layer_e_budget.o: src_teb/road_wall_layer_e_budget.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/road_wall_layer_e_budget.F90
roof_layer_e_budget.o: src_teb/roof_layer_e_budget.F90 modd_csts.o mode_thermos.o modi_layer_e_budget.o modi_layer_e_budget_get_coef.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/roof_layer_e_budget.F90
snow_cover_1layer.o: src_teb/snow_cover_1layer.F90 modd_csts.o modd_snow_par.o modd_surf_par.o mode_thermos.o modi_surface_ri.o modi_surface_aero_cond.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/snow_cover_1layer.F90
sunpos.o: src_solar/sunpos.F90 modd_csts.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_solar/sunpos.F90
surface_aero_cond.o: src_teb/surface_aero_cond.F90 modd_csts.o modi_wind_threshold.o mode_thermos.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/surface_aero_cond.F90
surface_cd.o: src_teb/surface_cd.F90 modd_csts.o mode_thermos.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/surface_cd.F90
surface_ri.o: src_teb/surface_ri.F90 modd_csts.o modd_surf_atm.o modi_wind_threshold.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/surface_ri.F90
teb.o: src_teb/teb.F90 modd_csts.o modd_surf_par.o modd_snow_par.o mode_thermos.o mode_surf_snow_frac.o modi_snow_cover_1layer.o modi_urban_drag.o modi_urban_snow_evol.o modi_roof_layer_e_budget.o modi_road_layer_e_budget.o modi_facade_e_budget.o modi_urban_fluxes.o modi_urban_hydro.o modi_bld_e_budget.o modi_wind_threshold.o modi_bem.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/teb.F90
teb_garden.o: src_teb/teb_garden.F90 modd_type_date_surf.o modd_csts.o modd_surf_par.o modd_snow_par.o mode_thermos.o mode_surf_snow_frac.o modi_garden_properties.o modi_greenroof_properties.o modi_window_shading_availability.o modi_urban_solar_abs.o modi_urban_lw_coef.o modi_garden.o modi_greenroof.o modi_teb.o modi_avg_urban_fluxes.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/teb_garden.F90
tester.o: src_driver/tester.F90 modd_csts.o modd_surf_atm.o modd_surf_par.o modd_type_date_surf.o mode_thermos.o sunpos.o ol_read_atm.o ol_alloc_atm.o ol_time_interp_atm.o modi_teb_garden.o modi_bem_morpho.o modi_window_data.o circumsolar_rad.o modd_forc_atm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_driver/tester.F90
tridiag_ground.o: src_teb/tridiag_ground.F90 hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/tridiag_ground.F90
urban_drag.o: src_teb/urban_drag.F90 modd_surf_par.o modd_csts.o mode_thermos.o modi_urban_exch_coef.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_drag.F90
urban_exch_coef.o: src_teb/urban_exch_coef.F90 modi_surface_ri.o modi_surface_cd.o modi_surface_aero_cond.o modi_wind_threshold.o modd_csts.o hook.o modi_flxsurf3bx.o modi_init_surfconsphy.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_exch_coef.F90
urban_fluxes.o: src_teb/urban_fluxes.F90 modd_surf_par.o modd_csts.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_fluxes.F90
urban_hydro.o: src_teb/urban_hydro.F90 modd_csts.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_hydro.F90
urban_lw_coef.o: src_teb/urban_lw_coef.F90 modd_csts.o modd_surf_par.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_lw_coef.F90
urban_snow_evol.o: src_teb/urban_snow_evol.F90 modd_snow_par.o modd_csts.o mode_surf_snow_frac.o modi_snow_cover_1layer.o modd_surf_par.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_snow_evol.F90
urban_solar_abs.o: src_teb/urban_solar_abs.F90 modd_csts.o modd_bem_cst.o modd_surf_par.o modi_window_shading.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/urban_solar_abs.F90
vslog.o: src_teb/vslog.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/vslog.f
wall_layer_e_budget.o: src_teb/wall_layer_e_budget.F90 modd_csts.o modi_layer_e_budget_get_coef.o modi_layer_e_budget.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/wall_layer_e_budget.F90
wind_threshold.o: src_teb/wind_threshold.F90 modd_surf_atm.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/wind_threshold.F90
window_data.o: src_teb/window_data.F90 hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/window_data.F90
window_e_budget.o: src_teb/window_e_budget.F90 modd_csts.o mode_conv_DOE.o hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/window_e_budget.F90
window_shading.o: src_teb/window_shading.F90 hook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/window_shading.F90
window_shading_availability.o: src_teb/window_shading_availability.F90 modd_bem_cst.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src_teb/window_shading_availability.F90
./modi_greenroof.F90: src_proxi_SVAT/modi_greenroof.F90
	cp src_proxi_SVAT/modi_greenroof.F90 .
./dx_air_cooling_coil_cv.F90: src_teb/dx_air_cooling_coil_cv.F90
	cp src_teb/dx_air_cooling_coil_cv.F90 .
./urban_snow_evol.F90: src_teb/urban_snow_evol.F90
	cp src_teb/urban_snow_evol.F90 .
./modi_urban_lw_coef.F90: src_teb/modi_urban_lw_coef.F90
	cp src_teb/modi_urban_lw_coef.F90 .
./mass_layer_e_budget.F90: src_teb/mass_layer_e_budget.F90
	cp src_teb/mass_layer_e_budget.F90 .
./modi_urban_exch_coef.F90: src_teb/modi_urban_exch_coef.F90
	cp src_teb/modi_urban_exch_coef.F90 .
./bem.F90: src_teb/bem.F90
	cp src_teb/bem.F90 .
./modd_bem_cst.F90: src_teb/modd_bem_cst.F90
	cp src_teb/modd_bem_cst.F90 .
./urban_fluxes.F90: src_teb/urban_fluxes.F90
	cp src_teb/urban_fluxes.F90 .
./modi_window_data.F90: src_teb/modi_window_data.F90
	cp src_teb/modi_window_data.F90 .
./init_surfconsphy.F: src_teb/init_surfconsphy.F
	cp src_teb/init_surfconsphy.F .
./modi_teb_garden.F90: src_teb/modi_teb_garden.F90
	cp src_teb/modi_teb_garden.F90 .
./modi_layer_e_budget_get_coef.F90: src_teb/modi_layer_e_budget_get_coef.F90
	cp src_teb/modi_layer_e_budget_get_coef.F90 .
./ol_read_atm.F90: src_driver/ol_read_atm.F90
	cp src_driver/ol_read_atm.F90 .
./modi_bld_e_budget.F90: src_teb/modi_bld_e_budget.F90
	cp src_teb/modi_bld_e_budget.F90 .
./mode_psychro.F90: src_teb/mode_psychro.F90
	cp src_teb/mode_psychro.F90 .
./mode_conv_DOE.F90: src_teb/mode_conv_DOE.F90
	cp src_teb/mode_conv_DOE.F90 .
./avg_urban_fluxes.F90: src_teb/avg_urban_fluxes.F90
	cp src_teb/avg_urban_fluxes.F90 .
./circumsolar_rad.F90: src_solar/circumsolar_rad.F90
	cp src_solar/circumsolar_rad.F90 .
./modi_urban_drag.F90: src_teb/modi_urban_drag.F90
	cp src_teb/modi_urban_drag.F90 .
./window_e_budget.F90: src_teb/window_e_budget.F90
	cp src_teb/window_e_budget.F90 .
./greenroof_properties.F90: src_proxi_SVAT/greenroof_properties.F90
	cp src_proxi_SVAT/greenroof_properties.F90 .
./modi_greenroof_properties.F90: src_proxi_SVAT/modi_greenroof_properties.F90
	cp src_proxi_SVAT/modi_greenroof_properties.F90 .
./modi_flxsurf3bx.f: src_teb/modi_flxsurf3bx.f
	cp src_teb/modi_flxsurf3bx.f .
./teb_garden.F90: src_teb/teb_garden.F90
	cp src_teb/teb_garden.F90 .
./road_wall_layer_e_budget.F90: src_teb/road_wall_layer_e_budget.F90
	cp src_teb/road_wall_layer_e_budget.F90 .
./road_layer_e_budget.F90: src_teb/road_layer_e_budget.F90
	cp src_teb/road_layer_e_budget.F90 .
./window_shading_availability.F90: src_teb/window_shading_availability.F90
	cp src_teb/window_shading_availability.F90 .
./modi_snow_cover_1layer.F90: src_teb/modi_snow_cover_1layer.F90
	cp src_teb/modi_snow_cover_1layer.F90 .
./modd_surf_atm.F90: src_teb/modd_surf_atm.F90
	cp src_teb/modd_surf_atm.F90 .
./ol_read_atm_ascii.F90: src_driver/ol_read_atm_ascii.F90
	cp src_driver/ol_read_atm_ascii.F90 .
./window_shading.F90: src_teb/window_shading.F90
	cp src_teb/window_shading.F90 .
./modi_urban_snow_evol.F90: src_teb/modi_urban_snow_evol.F90
	cp src_teb/modi_urban_snow_evol.F90 .
./modi_window_e_budget.F90: src_teb/modi_window_e_budget.F90
	cp src_teb/modi_window_e_budget.F90 .
./modi_surface_cd.F90: src_teb/modi_surface_cd.F90
	cp src_teb/modi_surface_cd.F90 .
./modd_forc_atm.F90: src_driver/modd_forc_atm.F90
	cp src_driver/modd_forc_atm.F90 .
./modd_snow_par.F90: src_teb/modd_snow_par.F90
	cp src_teb/modd_snow_par.F90 .
./close_file.F90: src_driver/close_file.F90
	cp src_driver/close_file.F90 .
./modi_init_surfconsphy.f: src_teb/modi_init_surfconsphy.f
	cp src_teb/modi_init_surfconsphy.f .
./tester.F90: src_driver/tester.F90
	cp src_driver/tester.F90 .
./sunpos.F90: src_solar/sunpos.F90
	cp src_solar/sunpos.F90 .
./modi_ini_csts.F90: src_teb/modi_ini_csts.F90
	cp src_teb/modi_ini_csts.F90 .
./modi_dx_air_cooling_coil_cv.F90: src_teb/modi_dx_air_cooling_coil_cv.F90
	cp src_teb/modi_dx_air_cooling_coil_cv.F90 .
./modd_water_par.F90: src_teb/modd_water_par.F90
	cp src_teb/modd_water_par.F90 .
./modd_flood_par.F90: src_teb/modd_flood_par.F90
	cp src_teb/modd_flood_par.F90 .
./modi_teb.F90: src_teb/modi_teb.F90
	cp src_teb/modi_teb.F90 .
./modi_layer_e_budget.F90: src_teb/modi_layer_e_budget.F90
	cp src_teb/modi_layer_e_budget.F90 .
./modi_window_shading_availability.F90: src_teb/modi_window_shading_availability.F90
	cp src_teb/modi_window_shading_availability.F90 .
./modi_surface_aero_cond.F90: src_teb/modi_surface_aero_cond.F90
	cp src_teb/modi_surface_aero_cond.F90 .
./tridiag_ground.F90: src_teb/tridiag_ground.F90
	cp src_teb/tridiag_ground.F90 .
./modd_csts.F90: src_teb/modd_csts.F90
	cp src_teb/modd_csts.F90 .
./window_data.F90: src_teb/window_data.F90
	cp src_teb/window_data.F90 .
./layer_e_budget.F90: src_teb/layer_e_budget.F90
	cp src_teb/layer_e_budget.F90 .
./mode_thermos.F90: src_teb/mode_thermos.F90
	cp src_teb/mode_thermos.F90 .
./modi_urban_solar_abs.F90: src_teb/modi_urban_solar_abs.F90
	cp src_teb/modi_urban_solar_abs.F90 .
./roof_layer_e_budget.F90: src_teb/roof_layer_e_budget.F90
	cp src_teb/roof_layer_e_budget.F90 .
./modi_surface_ri.F90: src_teb/modi_surface_ri.F90
	cp src_teb/modi_surface_ri.F90 .
./modi_wind_threshold.F90: src_teb/modi_wind_threshold.F90
	cp src_teb/modi_wind_threshold.F90 .
./modi_urban_hydro.F90: src_teb/modi_urban_hydro.F90
	cp src_teb/modi_urban_hydro.F90 .
./ini_csts.F90: src_teb/ini_csts.F90
	cp src_teb/ini_csts.F90 .
./wall_layer_e_budget.F90: src_teb/wall_layer_e_budget.F90
	cp src_teb/wall_layer_e_budget.F90 .
./snow_cover_1layer.F90: src_teb/snow_cover_1layer.F90
	cp src_teb/snow_cover_1layer.F90 .
./hook.F90: src_teb/hook.F90
	cp src_teb/hook.F90 .
./modi_mass_layer_e_budget.F90: src_teb/modi_mass_layer_e_budget.F90
	cp src_teb/modi_mass_layer_e_budget.F90 .
./modd_arch.F90: src_driver/modd_arch.F90
	cp src_driver/modd_arch.F90 .
./modi_facade_e_budget.F90: src_teb/modi_facade_e_budget.F90
	cp src_teb/modi_facade_e_budget.F90 .
./modi_tridiag_ground.F90: src_teb/modi_tridiag_ground.F90
	cp src_teb/modi_tridiag_ground.F90 .
./bem_morpho.F90: src_teb/bem_morpho.F90
	cp src_teb/bem_morpho.F90 .
./urban_lw_coef.F90: src_teb/urban_lw_coef.F90
	cp src_teb/urban_lw_coef.F90 .
./surface_ri.F90: src_teb/surface_ri.F90
	cp src_teb/surface_ri.F90 .
./bld_e_budget.F90: src_teb/bld_e_budget.F90
	cp src_teb/bld_e_budget.F90 .
./facade_e_budget.F90: src_teb/facade_e_budget.F90
	cp src_teb/facade_e_budget.F90 .
./modd_type_date_surf.F90: src_teb/modd_type_date_surf.F90
	cp src_teb/modd_type_date_surf.F90 .
./modi_floor_layer_e_budget.F90: src_teb/modi_floor_layer_e_budget.F90
	cp src_teb/modi_floor_layer_e_budget.F90 .
./open_close_bin_asc_forc.F90: src_driver/open_close_bin_asc_forc.F90
	cp src_driver/open_close_bin_asc_forc.F90 .
./close_file_asc.F90: src_driver/close_file_asc.F90
	cp src_driver/close_file_asc.F90 .
./mode_char2real.F90: src_driver/mode_char2real.F90
	cp src_driver/mode_char2real.F90 .
./modi_window_shading.F90: src_teb/modi_window_shading.F90
	cp src_teb/modi_window_shading.F90 .
./flxsurf3bx.F: src_teb/flxsurf3bx.F
	cp src_teb/flxsurf3bx.F .
./mode_surf_snow_frac.F90: src_teb/mode_surf_snow_frac.F90
	cp src_teb/mode_surf_snow_frac.F90 .
./modi_avg_urban_fluxes.F90: src_teb/modi_avg_urban_fluxes.F90
	cp src_teb/modi_avg_urban_fluxes.F90 .
./floor_layer_e_budget.F90: src_teb/floor_layer_e_budget.F90
	cp src_teb/floor_layer_e_budget.F90 .
./modi_garden_properties.F90: src_proxi_SVAT/modi_garden_properties.F90
	cp src_proxi_SVAT/modi_garden_properties.F90 .
./abor1_sfx.F90: src_driver/abor1_sfx.F90
	cp src_driver/abor1_sfx.F90 .
./urban_hydro.F90: src_teb/urban_hydro.F90
	cp src_teb/urban_hydro.F90 .
./modi_bem_morpho.F90: src_teb/modi_bem_morpho.F90
	cp src_teb/modi_bem_morpho.F90 .
./urban_exch_coef.F90: src_teb/urban_exch_coef.F90
	cp src_teb/urban_exch_coef.F90 .
./modd_surf_conf.F90: src_driver/modd_surf_conf.F90
	cp src_driver/modd_surf_conf.F90 .
./vslog.f: src_teb/vslog.f
	cp src_teb/vslog.f .
./surface_cd.F90: src_teb/surface_cd.F90
	cp src_teb/surface_cd.F90 .
./modi_roof_layer_e_budget.F90: src_teb/modi_roof_layer_e_budget.F90
	cp src_teb/modi_roof_layer_e_budget.F90 .
./urban_solar_abs.F90: src_teb/urban_solar_abs.F90
	cp src_teb/urban_solar_abs.F90 .
./modi_road_layer_e_budget.F90: src_teb/modi_road_layer_e_budget.F90
	cp src_teb/modi_road_layer_e_budget.F90 .
./greenroof.F90: src_proxi_SVAT/greenroof.F90
	cp src_proxi_SVAT/greenroof.F90 .
./layer_e_budget_get_coef.F90: src_teb/layer_e_budget_get_coef.F90
	cp src_teb/layer_e_budget_get_coef.F90 .
./modd_surf_par.F90: src_teb/modd_surf_par.F90
	cp src_teb/modd_surf_par.F90 .
./modi_garden.F90: src_proxi_SVAT/modi_garden.F90
	cp src_proxi_SVAT/modi_garden.F90 .
./add_forecast_to_date_surf.F90: src_driver/add_forecast_to_date_surf.F90
	cp src_driver/add_forecast_to_date_surf.F90 .
./garden_properties.F90: src_proxi_SVAT/garden_properties.F90
	cp src_proxi_SVAT/garden_properties.F90 .
./urban_drag.F90: src_teb/urban_drag.F90
	cp src_teb/urban_drag.F90 .
./surface_aero_cond.F90: src_teb/surface_aero_cond.F90
	cp src_teb/surface_aero_cond.F90 .
./modi_urban_fluxes.F90: src_teb/modi_urban_fluxes.F90
	cp src_teb/modi_urban_fluxes.F90 .
./teb.F90: src_teb/teb.F90
	cp src_teb/teb.F90 .
./ol_time_interp_atm.F90: src_driver/ol_time_interp_atm.F90
	cp src_driver/ol_time_interp_atm.F90 .
./ol_alloc_atm.F90: src_driver/ol_alloc_atm.F90
	cp src_driver/ol_alloc_atm.F90 .
./garden.F90: src_proxi_SVAT/garden.F90
	cp src_proxi_SVAT/garden.F90 .
./modi_wall_layer_e_budget.F90: src_teb/modi_wall_layer_e_budget.F90
	cp src_teb/modi_wall_layer_e_budget.F90 .
./wind_threshold.F90: src_teb/wind_threshold.F90
	cp src_teb/wind_threshold.F90 .
./modi_bem.F90: src_teb/modi_bem.F90
	cp src_teb/modi_bem.F90 .
./read_surf_atm.F90: src_driver/read_surf_atm.F90
	cp src_driver/read_surf_atm.F90 .
SRC = src_teb/modd_water_par.F90 src_teb/vslog.f src_teb/modi_urban_fluxes.F90 src_teb/layer_e_budget.F90 src_driver/open_close_bin_asc_forc.F90 src_driver/close_file_asc.F90 src_teb/window_shading.F90 src_teb/road_wall_layer_e_budget.F90 src_teb/surface_ri.F90 src_teb/modi_urban_hydro.F90 src_teb/modd_bem_cst.F90 src_teb/modi_snow_cover_1layer.F90 src_driver/abor1_sfx.F90 src_driver/add_forecast_to_date_surf.F90 src_teb/bem.F90 src_teb/modi_wall_layer_e_budget.F90 src_teb/road_layer_e_budget.F90 src_teb/modi_facade_e_budget.F90 src_proxi_SVAT/garden.F90 src_teb/urban_exch_coef.F90 src_teb/urban_drag.F90 src_driver/modd_forc_atm.F90 src_proxi_SVAT/modi_greenroof.F90 src_teb/modd_surf_par.F90 src_teb/modi_tridiag_ground.F90 src_teb/modi_window_shading_availability.F90 src_teb/teb.F90 src_teb/modd_flood_par.F90 src_driver/close_file.F90 src_driver/mode_char2real.F90 src_teb/dx_air_cooling_coil_cv.F90 src_teb/modi_floor_layer_e_budget.F90 src_teb/bld_e_budget.F90 src_solar/circumsolar_rad.F90 src_teb/flxsurf3bx.F src_solar/sunpos.F90 src_teb/modi_dx_air_cooling_coil_cv.F90 src_teb/modi_urban_snow_evol.F90 src_teb/tridiag_ground.F90 src_proxi_SVAT/modi_greenroof_properties.F90 src_teb/modi_urban_lw_coef.F90 src_teb/mass_layer_e_budget.F90 src_teb/modi_bld_e_budget.F90 src_driver/ol_alloc_atm.F90 src_teb/modi_init_surfconsphy.f src_driver/read_surf_atm.F90 src_teb/modi_teb_garden.F90 src_teb/wind_threshold.F90 src_teb/window_e_budget.F90 src_teb/mode_conv_DOE.F90 src_teb/mode_surf_snow_frac.F90 src_teb/modd_surf_atm.F90 src_teb/modi_road_layer_e_budget.F90 src_driver/tester.F90 src_teb/mode_psychro.F90 src_teb/modi_roof_layer_e_budget.F90 src_teb/surface_aero_cond.F90 src_teb/window_data.F90 src_teb/modi_flxsurf3bx.f src_teb/init_surfconsphy.F src_teb/avg_urban_fluxes.F90 src_teb/modi_bem_morpho.F90 src_driver/modd_surf_conf.F90 src_teb/urban_fluxes.F90 src_teb/modd_csts.F90 src_teb/modi_layer_e_budget.F90 src_proxi_SVAT/modi_garden.F90 src_teb/snow_cover_1layer.F90 src_teb/urban_snow_evol.F90 src_driver/ol_read_atm_ascii.F90 src_teb/modd_snow_par.F90 src_driver/ol_read_atm.F90 src_proxi_SVAT/modi_garden_properties.F90 src_teb/surface_cd.F90 src_teb/floor_layer_e_budget.F90 src_teb/mode_thermos.F90 src_proxi_SVAT/greenroof.F90 src_teb/wall_layer_e_budget.F90 src_teb/urban_solar_abs.F90 src_teb/modd_type_date_surf.F90 src_teb/modi_teb.F90 src_teb/modi_wind_threshold.F90 src_teb/bem_morpho.F90 src_teb/modi_window_shading.F90 src_teb/modi_urban_exch_coef.F90 src_teb/ini_csts.F90 src_teb/facade_e_budget.F90 src_teb/modi_window_e_budget.F90 src_teb/urban_hydro.F90 src_teb/roof_layer_e_budget.F90 src_teb/modi_ini_csts.F90 src_teb/modi_avg_urban_fluxes.F90 src_teb/urban_lw_coef.F90 src_teb/modi_bem.F90 src_driver/modd_arch.F90 src_teb/layer_e_budget_get_coef.F90 src_teb/modi_mass_layer_e_budget.F90 src_teb/hook.F90 src_teb/modi_layer_e_budget_get_coef.F90 src_teb/modi_urban_drag.F90 src_driver/ol_time_interp_atm.F90 src_teb/teb_garden.F90 src_teb/modi_window_data.F90 src_teb/modi_surface_cd.F90 src_teb/modi_surface_aero_cond.F90 src_proxi_SVAT/greenroof_properties.F90 src_teb/window_shading_availability.F90 src_teb/modi_surface_ri.F90 src_teb/modi_urban_solar_abs.F90 src_proxi_SVAT/garden_properties.F90
OBJ = modd_water_par.o vslog.o modi_urban_fluxes.o layer_e_budget.o open_close_bin_asc_forc.o close_file_asc.o window_shading.o road_wall_layer_e_budget.o surface_ri.o modi_urban_hydro.o modd_bem_cst.o modi_snow_cover_1layer.o abor1_sfx.o add_forecast_to_date_surf.o bem.o modi_wall_layer_e_budget.o road_layer_e_budget.o modi_facade_e_budget.o garden.o urban_exch_coef.o urban_drag.o modd_forc_atm.o modi_greenroof.o modd_surf_par.o modi_tridiag_ground.o modi_window_shading_availability.o teb.o modd_flood_par.o close_file.o mode_char2real.o dx_air_cooling_coil_cv.o modi_floor_layer_e_budget.o bld_e_budget.o circumsolar_rad.o flxsurf3bx.o sunpos.o modi_dx_air_cooling_coil_cv.o modi_urban_snow_evol.o tridiag_ground.o modi_greenroof_properties.o modi_urban_lw_coef.o mass_layer_e_budget.o modi_bld_e_budget.o ol_alloc_atm.o modi_init_surfconsphy.o read_surf_atm.o modi_teb_garden.o wind_threshold.o window_e_budget.o mode_conv_DOE.o mode_surf_snow_frac.o modd_surf_atm.o modi_road_layer_e_budget.o tester.o mode_psychro.o modi_roof_layer_e_budget.o surface_aero_cond.o window_data.o modi_flxsurf3bx.o init_surfconsphy.o avg_urban_fluxes.o modi_bem_morpho.o modd_surf_conf.o urban_fluxes.o modd_csts.o modi_layer_e_budget.o modi_garden.o snow_cover_1layer.o urban_snow_evol.o ol_read_atm_ascii.o modd_snow_par.o ol_read_atm.o modi_garden_properties.o surface_cd.o floor_layer_e_budget.o mode_thermos.o greenroof.o wall_layer_e_budget.o urban_solar_abs.o modd_type_date_surf.o modi_teb.o modi_wind_threshold.o bem_morpho.o modi_window_shading.o modi_urban_exch_coef.o ini_csts.o facade_e_budget.o modi_window_e_budget.o urban_hydro.o roof_layer_e_budget.o modi_ini_csts.o modi_avg_urban_fluxes.o urban_lw_coef.o modi_bem.o modd_arch.o layer_e_budget_get_coef.o modi_mass_layer_e_budget.o hook.o modi_layer_e_budget_get_coef.o modi_urban_drag.o ol_time_interp_atm.o teb_garden.o modi_window_data.o modi_surface_cd.o modi_surface_aero_cond.o greenroof_properties.o window_shading_availability.o modi_surface_ri.o modi_urban_solar_abs.o garden_properties.o
OFF = src_proxi_SVAT/modi_greenroof.F90 src_teb/dx_air_cooling_coil_cv.F90 src_teb/urban_snow_evol.F90 src_teb/modi_urban_lw_coef.F90 src_teb/mass_layer_e_budget.F90 src_teb/modi_urban_exch_coef.F90 src_teb/bem.F90 src_teb/modd_bem_cst.F90 src_teb/urban_fluxes.F90 src_teb/modi_window_data.F90 src_teb/init_surfconsphy.F src_teb/modi_teb_garden.F90 src_teb/modi_layer_e_budget_get_coef.F90 src_driver/ol_read_atm.F90 src_teb/modi_bld_e_budget.F90 src_teb/mode_psychro.F90 src_teb/mode_conv_DOE.F90 src_teb/avg_urban_fluxes.F90 src_solar/circumsolar_rad.F90 src_teb/modi_urban_drag.F90 src_teb/window_e_budget.F90 src_proxi_SVAT/greenroof_properties.F90 src_proxi_SVAT/modi_greenroof_properties.F90 src_teb/modi_flxsurf3bx.f src_teb/teb_garden.F90 src_teb/road_wall_layer_e_budget.F90 src_teb/road_layer_e_budget.F90 src_teb/window_shading_availability.F90 src_teb/modi_snow_cover_1layer.F90 src_teb/modd_surf_atm.F90 src_driver/ol_read_atm_ascii.F90 src_teb/window_shading.F90 src_teb/modi_urban_snow_evol.F90 src_teb/modi_window_e_budget.F90 src_teb/modi_surface_cd.F90 src_driver/modd_forc_atm.F90 src_teb/modd_snow_par.F90 src_driver/close_file.F90 src_teb/modi_init_surfconsphy.f src_driver/tester.F90 src_solar/sunpos.F90 src_teb/modi_ini_csts.F90 src_teb/modi_dx_air_cooling_coil_cv.F90 src_teb/modd_water_par.F90 src_teb/modd_flood_par.F90 src_teb/modi_teb.F90 src_teb/modi_layer_e_budget.F90 src_teb/modi_window_shading_availability.F90 src_teb/modi_surface_aero_cond.F90 src_teb/tridiag_ground.F90 src_teb/modd_csts.F90 src_teb/window_data.F90 src_teb/layer_e_budget.F90 src_teb/mode_thermos.F90 src_teb/modi_urban_solar_abs.F90 src_teb/roof_layer_e_budget.F90 src_teb/modi_surface_ri.F90 src_teb/modi_wind_threshold.F90 src_teb/modi_urban_hydro.F90 src_teb/ini_csts.F90 src_teb/wall_layer_e_budget.F90 src_teb/snow_cover_1layer.F90 src_teb/hook.F90 src_teb/modi_mass_layer_e_budget.F90 src_driver/modd_arch.F90 src_teb/modi_facade_e_budget.F90 src_teb/modi_tridiag_ground.F90 src_teb/bem_morpho.F90 src_teb/urban_lw_coef.F90 src_teb/surface_ri.F90 src_teb/bld_e_budget.F90 src_teb/facade_e_budget.F90 src_teb/modd_type_date_surf.F90 src_teb/modi_floor_layer_e_budget.F90 src_driver/open_close_bin_asc_forc.F90 src_driver/close_file_asc.F90 src_driver/mode_char2real.F90 src_teb/modi_window_shading.F90 src_teb/flxsurf3bx.F src_teb/mode_surf_snow_frac.F90 src_teb/modi_avg_urban_fluxes.F90 src_teb/floor_layer_e_budget.F90 src_proxi_SVAT/modi_garden_properties.F90 src_driver/abor1_sfx.F90 src_teb/urban_hydro.F90 src_teb/modi_bem_morpho.F90 src_teb/urban_exch_coef.F90 src_driver/modd_surf_conf.F90 src_teb/vslog.f src_teb/surface_cd.F90 src_teb/modi_roof_layer_e_budget.F90 src_teb/urban_solar_abs.F90 src_teb/modi_road_layer_e_budget.F90 src_proxi_SVAT/greenroof.F90 src_teb/layer_e_budget_get_coef.F90 src_teb/modd_surf_par.F90 src_proxi_SVAT/modi_garden.F90 src_driver/add_forecast_to_date_surf.F90 src_proxi_SVAT/garden_properties.F90 src_teb/urban_drag.F90 src_teb/surface_aero_cond.F90 src_teb/modi_urban_fluxes.F90 src_teb/teb.F90 src_driver/ol_time_interp_atm.F90 src_driver/ol_alloc_atm.F90 src_proxi_SVAT/garden.F90 src_teb/modi_wall_layer_e_budget.F90 src_teb/wind_threshold.F90 src_teb/modi_bem.F90 src_driver/read_surf_atm.F90
clean: neat
	-rm -f .cppdefs $(OBJ) tester.exe
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
tester.exe: $(OBJ) 
	$(LD) $(OBJ) -o tester.exe  $(LDFLAGS)
