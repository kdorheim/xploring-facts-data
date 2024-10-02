# Get separate climate variable files

library(ncdf4)

# Get needed vars from 2000 sample hector run file
nc_all <- nc_open(file.path("data","OUTPUTS","hector_results","hector.tlm.global","hector.tlm.global.temperature.hector.temperature_climate.nc"))
oceantemp <- ncvar_get(nc_all,"ssp585/deep_ocean_temperature")
ohc <- ncvar_get(nc_all,"ssp585/ocean_heat_content")
sfc_temp <- ncvar_get(nc_all,"ssp585/surface_temperature")
years <- ncvar_get(nc_all,"ssp585/years")
samples <- ncvar_get(nc_all,"ssp585/samples")
nc_close(nc_all)

# Replace values in FaIR files with Hector
nc_gsat <- nc_open(file.path("data","OUTPUTS","hector_results","hector.tlm.global","hector.tlm.global.temperature.hector.temperature_gsat.nc"),write=TRUE)
ncvar_put(nc_gsat,varid="surface_temperature",vals=sfc_temp)
nc_close(nc_gsat)

nc_oceantemp <- nc_open(file.path("data","OUTPUTS","hector_results","hector.tlm.global","hector.tlm.global.temperature.hector.temperature_oceantemp.nc"),write=TRUE)
ncvar_put(nc_oceantemp,varid="deep_ocean_temperature",vals=oceantemp)
nc_close(nc_oceantemp)

nc_ohc <- nc_open(file.path("data","OUTPUTS","hector_results","hector.tlm.global","hector.tlm.global.temperature.hector.temperature_ohc.nc"),write=TRUE)
ncvar_put(nc_ohc,varid="ocean_heat_content",vals=ohc)
nc_close(nc_ohc)
