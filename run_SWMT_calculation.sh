#!/bin/bash

## loop over count, submit job to gadi with count that gets communicated to python

for y in {2159..2159}; do
   qsub -v year=$y,expt='01deg_jra55v13_ryf9091_rerun_for_easterlies' SWMT_calculation.sh
done

# for y in {2156..2159}; do
#    qsub -v year=$y,expt='access-om2-01_ryf_wind_50_down_zonal' SWMT_calculation.sh
# done


# for y in {2156..2159}; do
#    qsub -v year=$y,expt='access-om2-01_ryf_meltwater_50_down' SWMT_calculation.sh
# done

# for y in {2150..2153}; do
#    qsub -v year=$y,expt='access-om2-01_ryf_meltwater_78_down_excl_Getz' SWMT_calculation.sh
# done

# for y in {2153..2155}; do
#    qsub -v year=$y,expt='access-om2-01_ryf_meltwater_50_up' SWMT_calculation.sh
# done

# for y in {2150..2151}; do
#    qsub -v year=$y,expt='access-om2-01_ryf_wind_zonal_mw_50_down' SWMT_calculation.sh
# done


