#!/usr/bin/perl
# the inputs are the pxyz files generated from 01_runme.sh
foreach $file (@ARGV)
{
$basename=substr $file,0,-5;
@base_path=split(/\//, $basename);
$orca_input=@base_path[@base_path-1].".orca.inp";
$dirname=@base_path[@base_path-1];

### Create QM Region Orca File
$cmd1="pxyz_select_n $file 4 3 2 | grep -v 'MW' | xyz_add_linenu | pxyz_orca_upd_xyz MyMol.templ.inp nam_translation > hold";
### Update Charges in Orca File
$cmd2="pxyz_select 1 $file | grep -v 'OW' | xyz_add_linenu | pxyz_orca_upd_chg hold chginfo > $orca_input && rm hold";
$cmd3="mkdir $dirname";
$cmd4="mv $orca_input $dirname";
$cmd5="cp -a runme.pinnacle $dirname";
system("$cmd1");
system("$cmd2");
system("$cmd3");
system("$cmd4");
system("$cmd5");
}
