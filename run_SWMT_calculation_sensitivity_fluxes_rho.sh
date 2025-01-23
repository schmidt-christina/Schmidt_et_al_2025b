#!/bin/bash

## loop over count, submit job to gadi with count that gets communicated to python

for y in {2156..2159}; do
   qsub -v year=$y,expt_flux='access-om2-01_ryf_wind_50_down_zonal',expt_rho='01deg_jra55v13_ryf9091_rerun_for_easterlies' SWMT_calculation_sensitivity_fluxes_rho.sh
done

for y in {2156..2159}; do
   qsub -v year=$y,expt_flux='01deg_jra55v13_ryf9091_rerun_for_easterlies',expt_rho='access-om2-01_ryf_wind_50_down_zonal' SWMT_calculation_sensitivity_fluxes_rho.sh
done


for y in {2156..2159}; do
   qsub -v year=$y,expt_flux='access-om2-01_ryf_meltwater_50_down',expt_rho='01deg_jra55v13_ryf9091_rerun_for_easterlies' SWMT_calculation_sensitivity_fluxes_rho.sh
done

for y in {2156..2159}; do
   qsub -v year=$y,expt_flux='01deg_jra55v13_ryf9091_rerun_for_easterlies',expt_rho='access-om2-01_ryf_meltwater_50_down' SWMT_calculation_sensitivity_fluxes_rho.sh
done

