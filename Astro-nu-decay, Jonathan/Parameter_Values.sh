#!/bin/bash

index_block=0
num_blocks=1

mnu1_min=0
mnu1_max=0.0281
mnu1_npts=10

gamma_min=2
gamma_max=3
gamma_npts=10

log10g21_min=-8
log10g21_max=0
log10g21_npts=10

log10g31_min=-8
log10g31_max=-0
log10g31_npts=10

log10g32_min=-8
log10g32_max=-0
log10g32_npts=10

s12sq_min=0.304
s12sq_max=0.304
s12sq_npts=1

s23sq_min=0.573
s23sq_max=0.573
s23sq_npts=1

s13sq=0.02219
deltacp=3.4382986264

F_e_s=0.3333333333
F_mu_s=0.6666666667

output_path=data
verbose=1

make clean
make

for ((index_block = 0 ; index_block <= $num_blocks ; index_block++)); 
do
  echo $index_block $num_blocks $mnu1_min $mnu1_max $mnu1_npts $gamma_min $gamma_max $gamma_npts $log10g21_min $log10g21_max $log10g21_npts $log10g31_min $log10g31_max $log10g31_npts $log10g32_min $log10g32_max $log10g32_npts $s12sq_min $s12sq_max $s12sq_npts $s23sq_min $s23sq_max $s23sq_npts $s13sq $deltacp $F_e_s $F_mu_s $output_path $verbose
  time ./main $index_block $num_blocks $mnu1_min $mnu1_max $mnu1_npts $gamma_min $gamma_max $gamma_npts $log10g21_min $log10g21_max $log10g21_npts $log10g31_min $log10g31_max $log10g31_npts $log10g32_min $log10g32_max $log10g32_npts $s12sq_min $s12sq_max $s12sq_npts $s23sq_min $s23sq_max $s23sq_npts $s13sq $deltacp $F_e_s $F_mu_s $output_path $verbose
done