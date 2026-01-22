#!/bin/sh

# Download ergo --> directory ergo-3.8.2
# In ergo-3.8.2 Run ./bootstrap.h && make
topdir=`pwd`
ergo_dirname=ergo-3.8.2
ergo_output_dirname=aspirin_mtx
rm -rf $ergo_output_dirname ; mkdir $ergo_output_dirname ;
cd $ergo_output_dirname

echo Running Aspirin HF/STO-3G
"$topdir"/"$ergo_dirname"/source/ergo -m "$topdir"/aspirin.xyz <<EOINPUT > /dev/null
basis = "STO-3G"
scf.calculation_identifier = "Aspirin"
run "HF"
EOINPUT
mv ergoscf.out ergoscf_sto3g.out

echo Running Aspirin HF/6-31Gss
"$topdir"/"$ergo_dirname"/source/ergo -m "$topdir"/aspirin.xyz <<EOINPUT > /dev/null
initial_density="density.bin"
basis = "6-31Gss"
scf.calculation_identifier = "Aspirin"
run "HF"
EOINPUT
mv ergoscf.out ergoscf_6-31gss.out

echo Running Aspirin HF/6-311++Gss
"$topdir"/"$ergo_dirname"/source/ergo -m "$topdir"/aspirin.xyz <<EOINPUT > /dev/null
initial_density="density.bin"
basis = "6-311++Gss"
scf.calculation_identifier = "Aspirin"
scf.create_mtx_file_S = 1
scf.create_mtx_files_F = 1
scf.create_mtx_files_D = 1
run "HF"
EOINPUT
