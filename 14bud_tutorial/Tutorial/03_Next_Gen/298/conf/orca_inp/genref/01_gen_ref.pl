#!/usr/bin/perl
# the input files are *.xyz files
foreach $file (@ARGV)
{
$basename=substr $file,0,-4;
@base_path=split(/\//, $basename);
$name=@base_path[@base_path-1];
$d3name=$name.".d3.xyz.d3grad";
$gradname=$name.".orca.engrad";
$cmd1="chunk 2 9999 $file | grep -v ' MW' | ref_gen_step1_cord molinfo | ref_upd_orca_grad $gradname | ref_upd_net | ref_upd_d3_grad $d3name | ref_upd_net | xyz_add_msite Ow Hw Hw M | ref_fix_msite M | xyz_add_msite Omm Hmm Hmm Mmm | ref_fix_msite Mmm | xyz_fix_linenu > ref_files/$name.ref";
system("$cmd1");
}
