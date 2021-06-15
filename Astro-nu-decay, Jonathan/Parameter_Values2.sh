#!/bin/bash
index_block=2
num_blocks=3

mnu1_min=0.01
mnu1_max=0.01
mnu1_npts=1

gamma_min=2
gamma_max=2
gamma_npts=1

log10g21_min=-10
log10g21_max=-10
log10g21_npts=1

log10g31_min=-10
log10g31_max=-10
log10g31_npts=1

log10g32_min=-8
log10g32_max=0
log10g32_npts=150

s12sq_min=0.304
s12sq_max=0.304
s12sq_npts=1

s23sq_min=0.573
s23sq_max=0.573
s23sq_npts=1

s13sq=0.02219
deltacp=3.4382986264

F_e_s=0.3333333333333333333
F_mu_s=0.6666666666666666667

output_path=py
verbose=1
ykbh=0

make clean
make
echo $index_block $num_blocks $mnu1_min $mnu1_max $mnu1_npts $gamma_min $gamma_max $gamma_npts $log10g21_min $log10g21_max $log10g21_npts $log10g31_min $log10g31_max $log10g31_npts $log10g32_min $log10g32_max $log10g32_npts $s12sq_min $s12sq_max $s12sq_npts $s23sq_min $s23sq_max $s23sq_npts $s13sq $deltacp $F_e_s $F_mu_s $output_path $verbose $ykbh
time ./main $index_block $num_blocks $mnu1_min $mnu1_max $mnu1_npts $gamma_min $gamma_max $gamma_npts $log10g21_min $log10g21_max $log10g21_npts $log10g31_min $log10g31_max $log10g31_npts $log10g32_min $log10g32_max $log10g32_npts $s12sq_min $s12sq_max $s12sq_npts $s23sq_min $s23sq_max $s23sq_npts $s13sq $deltacp $F_e_s $F_mu_s $output_path $verbose $ykbh