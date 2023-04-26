#!/usr/bin/perl
# the inputs are pxyz files generated from 01_runme.sh
foreach $file (@ARGV)
{
$basename=substr $file,0,-5;
@base_path=split(/\//, $basename);
$d3xyz=@base_path[@base_path-1].".d3.xyz";
$cmd1="pxyz_select_n $file 4 3 2 | grep -v MW | xyz_add_linenu | pxyz_pqs_gen_xyz nucinfo.d3 | xyz_add_linenu > $d3xyz";
$cmd2="getd3force b3-lyp $d3xyz";
$cmd3="mv $d3xyz.d3grad genref && rm $d3xyz.eng $d3xyz";
system("$cmd1");
system("$cmd2");
system("$cmd3");
}
