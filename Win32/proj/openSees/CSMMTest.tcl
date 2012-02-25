wipe;
model basic -ndm 2 -ndf 2;
set pi 3.141592654; set PlaneStressMatTag 20;
set wfc 14.3; set wfyv 300; set wfyh1 210; set wE 1.96e5;
set rouv 0.0246; set rouh1 0.01246; # #10@70
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
#uniaxialMaterial ConcreteL02 14 [expr -$wfc] -0.0025
#uniaxialMaterial ConcreteL02 15 [expr -$wfc] -0.0025
uniaxialMaterial ManderConcreteC01 14 [expr -$wfc] -0.002 30000. 0.40 0.60 [expr 0.1*$wfc] 8.0e-5 0.5 2.9 -spall 2.1
uniaxialMaterial ManderConcreteC01 15 [expr -$wfc] -0.002 30000. 0.40 0.60 [expr 0.1*$wfc] 8.0e-5 0.5 2.9 -spall 2.1
nDMaterial CSMMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002
node 1001 0.0 0.0;
node 1002 1000.0 0.0;
node 1003 1000.0 1000.0;
node 1004 0.0 1000.0;
fix 1001 1 1 ;
fix 1002 0 1 ;
#element quad 2001 1001 1002 1003 1004 150.0  "PlaneStress" $PlaneStressMatTag;
element SSPquad 2001 1001 1002 1003 1004 $PlaneStressMatTag "PlaneStress" 150. 0. 0.;
recorder Node -file disp.out -time -node 1003 1004 -dof 1 2 disp;
set numElem 1;
recorder Element -eleRange 1 $numElem -time -file stress.out  stress;
recorder Element -eleRange 1 $numElem -time -file strain.out  strain;
recorder display g3 10 10 800 600 -wipe;
prp 0 0 1000;  vrp 0 0 0;  vup 0 1 0;  vpn 0 0 1;  viewWindow -400 1400 -200 1200;
port -1 1 -1 1; projection 0; fill 0; display 1 -1 1000;
set axialLoad -5.e3;
pattern Plain 3001 "Constant" {    load 1003   0.0 $axialLoad;    load 1004   0.0 $axialLoad;};
integrator LoadControl 0.1
system UmfPack;	#SparseGEN; # Overkill, but may need the pivoting!SparseGeneral -piv
test EnergyIncr  1.0e-4 30; #test NormDispIncr 1.0e-3 250
numberer Plain
constraints Plain
algorithm Newton
analysis Static
analyze 10 0.1
#nodeDisp 1003;nodeDisp 1004;
loadConst 0.0
pattern Plain 3002 "Linear" {    load 1003  1.0 0.0;    load 1004  1.0 0.0;};
set IDctrlNode 1003;
set IDctrlDOF 1;
set Dmax 10;
set Dincr 0.1;
set numIncr 1; #100
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr;
set maxNumIterStatic 6;
set algorithmTypeStatic Newton;
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr $Dmax $Dincr $Dincr
analyze 5
#global LunitTXT;
#if {  [info exists LunitTXT] != 1} {set LunitTXT "Length"};
#
#set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Target=%.4f /%s";
#
#if {$ok != 0} {
#	set Dstep 0.0;
#	set ok 0
#	while {$Dstep <= 1.0 && $ok == 0} {
#		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
#		set Dstep [expr $controlDisp/$Dmax]
#		set ok [analyze 1];
#		if {$ok != 0} {;
#			set Nk 4;
#			set DincrReduced [expr $Dincr/$Nk];
#			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
#			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
#				set ok [analyze 1];
#				if {$ok != 0} {
#					puts "Trying Newton with Initial Tangent .."
#					test NormDispIncr   $TolStatic      2000 0
#					algorithm Newton -initial
#					set ok [analyze 1]
#					test $testTypeStatic $TolStatic      $maxNumIterStatic    0
#					algorithm $algorithmTypeStatic
#				}
#				if {$ok != 0} {
#					puts "Trying Broyden .."
#					algorithm Broyden 8
#					set ok [analyze 1 ]
#					algorithm $algorithmTypeStatic
#				}
#				if {$ok != 0} {
#					puts "Trying NewtonWithLineSearch .."
#					algorithm NewtonLineSearch 0.8
#					set ok [analyze 1]
#					algorithm $algorithmTypeStatic
#				}
#				if {$ok != 0} {;
#					puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
#					return -1
#				};
#			};
#			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;
#		};
#	};
#};
#
#if {$ok != 0 } {
#	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
#} else {
#	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
#}