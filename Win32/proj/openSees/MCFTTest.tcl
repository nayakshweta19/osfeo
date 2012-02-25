wipe; model BasicBuilder -ndm 2 -ndf 2; set dataDir Data; file mkdir $dataDir; source libUnits.tcl;
set PlaneStressMatTag 20;
set pi [expr 2.*asin(1.0)];
set wfc 30.0; set wfyv 550.0; set wfyh1 500.0; set wE 2.e5;
set epsc0 0.0025; set Ec [expr 2.0*$wfc/$epsc0]; set fpcu  46.0; set epsU  0.025;
set lambda 0.1; set ft  5.0; #[expr 0.65*pow($wfc,0.33)];
set Ets [expr $ft/0.002];
set rou1 0; set rou2 0; set db1 20; set db2 20; set b 0.015; set rouv 0.03785; set rouh1 0.02713;
set rn_pre 0.4; set rn_post 0.6; set ecr 0.00008; set rp 0.5; set xp_cr 2.9; set xn_cr 2.1;
uniaxialMaterial MCFTSteel03 11 $wfyv $wE $b
uniaxialMaterial MCFTSteel03 12 $wfyh1 $wE $b
uniaxialMaterial MCFTConcrete03 14 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets
uniaxialMaterial MCFTConcrete03 15 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets
#uniaxialMaterial ManderConcreteM01 14 [expr -$wfc] -$epsc0 $Ec $rn_pre $rn_post $ft $ecr $rp $xp_cr -spall $xn_cr;
#uniaxialMaterial ManderConcreteM01 15 [expr -$wfc] -$epsc0 $Ec $rn_pre $rn_post $ft $ecr $rp $xp_cr -spall $xn_cr;
nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial DSFMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
node 1001 0.0 0.0; node 1002 1000.0 0.0; node 1003 1000.0 1000.0; node 1004 0.0 1000.0;
fix 1001 1 1 ; fix 1002 0 1 ; fix 1003 0 0 ; fix 1004 0 0 ; 
#element quad 2001 1001 1002 1003 1004 150.0  "PlaneStress" $PlaneStressMatTag;
element SSPquad 2001 1001 1002 1003 1004 $PlaneStressMatTag "PlaneStress" 150. 0. 0.;
recorder Node -file disp.out -time -node 1003 1004 -dof 1 2 disp;
set numElem 1;
recorder Element -eleRange 1 $numElem -time -file stress.out  stress;
recorder Element -eleRange 1 $numElem -time -file strain.out  strain;
recorder display g3 10 10 800 600 -wipe;
prp 0 0 1000;  vrp 0 0 0;  vup 0 1 0;  vpn 0 0 1;  viewWindow -400 1400 -200 1200;
port -1 1 -1 1; projection 0; fill 0; display 1 -1 100;
set axialLoad -3.e6;
pattern Plain 3001 "Constant" {    load 1003   0.0 $axialLoad;    load 1004   0.0 $axialLoad;};
integrator LoadControl 0.1
system SparseGEN;	#UmfPack; # Overkill, but may need the pivoting!SparseGeneral -piv
test EnergyIncr  1.0e-9 200 4 ; #test NormDispIncr 1.0e-3 250
numberer Plain; constraints Plain; algorithm Newton; analysis Static;
analyze 10 0.1; #nodeDisp 1003; nodeDisp 1004;

loadConst 0.0
pattern Plain 3002 "Linear" {    load 1003  1.0 0.0;    load 1004  1.0 0.0;};
set IDctrlNode 1003;
set IDctrlDOF 1;
set Dmax 10;
set Dincr 1;
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr $Dmax $Dincr $Dincr
analyze 5
set Dmax -10;
set Dincr -0.1;
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr $Dmax $Dincr $Dincr
analyze 5

wipe; model BasicBuilder -ndm 2 -ndf 2; set dataDir Data; file mkdir $dataDir; source libUnits.tcl;
set PlaneStressMatTag 20; set pi [expr 2.*asin(1.0)];
set k       27777.78; #---bulk modulus
set G       9259.26;  #---shear modulus
set sig0    300.0;      #---yield stress
set sigInf  600.0;      #---final saturation yield stress
set delta   0.1;      #---exponential hardening parameter eta>= 0 for rate independent case
set H       0.05;     #---linear hardening parameter
nDMaterial J2Plasticity 1 $k $G $sig0 $sigInf $delta $H;
#set k       27777.78; #---bulk modulus
#set G       9259.26;  #---shear modulus
#set sigY    5.0;      #---yield stress
#set rho     0.398;    #---failure surface and associativity
#set rhoBar  0.398;    #
#set Kinf    0.0;      #---isotropic hardening
#set Ko 	    0.0
#set delta1  0.0
#set H 	    0.0;      #---kinematic hardening
#set theta   1.0
#set delta2  0.0;      #---tension softening
#nDMaterial DruckerPrager 1 7.8 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta
nDMaterial PlaneStress $PlaneStressMatTag 1;
node 1001 0.0 0.0; node 1002 1000.0 0.0; node 1003 1000.0 1000.0; node 1004 0.0 1000.0;
fix 1001 1 1 ; fix 1002 1 1 ; fix 1003 0 0 ; fix 1004 0 0 ; 
#element quad 2001 1001 1002 1003 1004 100.0  "PlaneStress" $PlaneStressMatTag;
element SSPquad 2001 1001 1002 1003 1004 $PlaneStressMatTag "PlaneStress" 100. 0. 0.;
recorder Node -file disp.out -time -node 1003 1004 -dof 1 2 disp;
set numElem 1;
recorder Element -eleRange 1 $numElem -time -file stress.out  stress;
recorder Element -eleRange 1 $numElem -time -file strain.out  strain;
#recorder display g3 10 10 800 600 -wipe;
#prp 0 0 1000;  vrp 0 0 0;  vup 0 1 0;  vpn 0 0 1;  viewWindow -400 1400 -200 1200;
#port -1 1 -1 1; projection 0; fill 0; display 1 -1 100;
set axialLoad -5.e7;
pattern Plain 3001 "Constant" {    load 1003   0.0 $axialLoad;    load 1004   0.0 $axialLoad;};
integrator LoadControl 0.1
system UmfPack;	#; SparseGEN# Overkill, but may need the pivoting!SparseGeneral -piv
test EnergyIncr  1.0e-4 200 4 ; #test NormDispIncr 1.0e-3 250
numberer Plain; constraints Plain; algorithm ModifiedNewton; analysis Static;
analyze 10 0.1; #nodeDisp 1003;nodeDisp 1004;
#set numIncr 1; #100
#set TolStatic 1.e-9;
#set testTypeStatic EnergyIncr;
#set maxNumIterStatic 6;
#set algorithmTypeStatic Newton;
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