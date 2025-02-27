#!/bin/bash
#PBS -P e14
#PBS -l ncpus=48
#PBS -l mem=188GB
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l storage="gdata/hh5+gdata/ik11+gdata/v45+gdata/e14+gdata/cj50+scratch/v45+scratch/x77"
#PBS -l wd
#PBS -o calculation_SWMT.out
#PBS -j oe

module use /g/data/hh5/public/modules
#module load conda/analysis3-unstable

python3 SWMT_calculation_sensitivity_fluxes_rho.py ${year} ${expt_flux} ${expt_rho} &>> SWMT_calculation_${expt_flux}_${expt_rho}_${year}.txt