wipe;
model BasicBuilder -ndm 2 -ndf 2;
set dataDir Data;
file mkdir $dataDir;
source libUnits.tcl;

set pi 3.141592654;
set nDMatTag 6;
set PlaneStressMatTag 20;
set testMatTag 101;
set testSecTag 201;
#nDMaterial ElasticIsotropic $matTag $E $v
nDMaterial ElasticIsotropic 1 1.0e4 0.499;
nDMaterial PlaneStress 2 1;
nDMaterial PlaneStrain 3 1;

##nDMaterial PlaneStress $PlaneStressMatTag $nDMatTag

set wfc 19.0; #30.0;
set wfyv 458.0;
set wfyh1 300.0;
set wE 2.e5;
set epsc0 0.00215;
set fpcu  46.0;
set epsU  0.025;
set lambda 0.1;
set ft  1.72; #[expr 0.14*$wfc];
set Ets [expr $ft/0.002];
set rou1 0;
set rou2 0;
set db1 16;
set db2 18;
set b 0.02
set rouv 0.01785;
set rouh1 0.00713;

uniaxialMaterial MCFTSteel03 11 $wfyv $wE $b
uniaxialMaterial MCFTSteel03 12 $wfyh1 $wE $b

uniaxialMaterial MCFTConcrete03 14 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets
uniaxialMaterial MCFTConcrete03 15 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets

nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial DSFMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0;
node 1002 1000.0 0.0;
node 1003 1000.0 1000.0
node 1004 0.0 1000.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 ; #1
fix 1002 0 1 ; #1

set PlaneStressMatTag 20
# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $testSecTag
element quad 2001 1001 1002 1003 1004 100.0  "PlaneStress" $PlaneStressMatTag
#element SSPquad 
#element SSPquad 2001 1001 1002 1003 1004 $PlaneStressMatTag "PlaneStress" 1. 0. 0.

# Create recorder
recorder Node -file disp.out -time -node 1003 1004 -dof 1 2 disp;
set numElem 1;
recorder Element -eleRange 1 $numElem -time -file stress.out  stress
recorder Element -eleRange 1 $numElem -time -file strain.out  strain

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp 0 0 1000
  vrp 0 0 0
  vup 0 1 0
  vpn 0 0 1
  viewWindow -400 1400 -200 1200
}
port -1 1 -1 1
projection 0
fill 0
display 1 -1 10


# Define constant axial load
set axialLoad -10;
pattern Plain 3001 "Constant" {
    #load 1002 $axialLoad 0.0 0.0
    load 1003   0.0 $axialLoad 
    load 1004   0.0 $axialLoad 
}

# Define analysis parameters
integrator LoadControl 0 1 0 0
system SparseGEN; #UmfPack;	# Overkill, but may need the pivoting!SparseGeneral -piv
test EnergyIncr  1.0e-4 30
#test NormDispIncr 1.0e-3 250
numberer Plain
constraints Plain
algorithm Newton
analysis Static

# Do one analysis for constant axial load
analyze 10 0.1
loadConst 0.0

pattern Plain 3002 "Linear" {
    load 1003  1.0 0.0
    load 1004  1.0 0.0
}
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

# Do the section analysis
set ok [analyze $numIncr]

global LunitTXT;
if {  [info exists LunitTXT] != 1} {set LunitTXT "Length"};
	
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Target=%.4f /%s";

if {$ok != 0} {  
	set Dstep 0.0;
	set ok 0
	while {$Dstep <= 1.0 && $ok == 0} {	
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr $controlDisp/$Dmax]
		set ok [analyze 1];
		if {$ok != 0} {;
			set Nk 4;
			set DincrReduced [expr $Dincr/$Nk];
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
				set ok [analyze 1];
				if {$ok != 0} {  
					puts "Trying Newton with Initial Tangent .."
					test NormDispIncr   $TolStatic      2000 0
					algorithm Newton -initial
					set ok [analyze 1]
					test $testTypeStatic $TolStatic      $maxNumIterStatic    0
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					puts "Trying Broyden .."
					algorithm Broyden 8
					set ok [analyze 1 ]
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					puts "Trying NewtonWithLineSearch .."
					algorithm NewtonLineSearch 0.8 
					set ok [analyze 1]
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {;
					puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					return -1
				};
			};
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;
		};
	};
};

if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}