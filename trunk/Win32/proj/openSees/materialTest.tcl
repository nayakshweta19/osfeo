# Build OpenSees Model
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 2;	# Define the model builder
set RigidFloor off

# Define nodes
node  1   0   0
node  2   0   1
# Single point constraints -- Boundary Conditions
fix   1   1  1 
fix   2   1  0 

set IDctrlNode 2
set IDctrlDOF 2

set MaterialTag 99
set ElasticMaterialTag 999

# set fc fy E
#set wfc 57.0;
#set wfy 419.2;
#set wE  187544.0;
set wfc 30.0;
set wfy 235.0;
set wE 2.06e5;
set rou1 0.0023;
set rou2 0.0023;

# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
#uniaxialMaterial    SteelZ01  $MaterialTag   $wfy    $wE  $wfc  $rou1
 
# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteZ01  $MaterialTag [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  $MaterialTag [expr -$wfc] -0.0025  
#uniaxialMaterial TendonL01    $MaterialTag $wfy $wE $wfc 0.015 0.02
#uniaxialMaterial Concrete01 $MaterialTag  -4.136854E+004  -4.000000E-003  -3.447379E+004  -1.400000E-002  +5.000000E-001  +4.136854E+003  +3.447379E+006 
# Material "Concrete02":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
#uniaxialMaterial  Concrete02   $MaterialTag  -4.136854E+004  -4.000000E-003  -3.447379E+004  -1.400000E-002  +5.000000E-001  +4.136854E+003  +3.447379E+006 
uniaxialMaterial MCFTConcrete01 $MaterialTag 30.0 0.0022

set Eelastic 1.0e-9
uniaxialMaterial Elastic $ElasticMaterialTag $Eelastic
set Area 1.
element truss 1 1 2 $Area $MaterialTag ;			# defined section behavior, same as above, since A=1
#element truss 2 1 2 $Area $ElasticMaterialTag ;			# defined section behavior, same as above, since A=1
set RecorderFilename test_Response.out
#recorder Node -file $RecorderFilename -time -node $IDctrlNode -dof $IDctrlDOF disp
#set RecorderFilename test_ElementForces.out
#recorder Element -file $RecorderFilename -ele 1 axialForce
#set RecorderFilename test_ElementDeformations.out
recorder Element -file $RecorderFilename -time -ele 1 deformations
#set fileName matTestData;
#set windowTitle "Concrete Response"
#set xLoc 10; set yLoc 10;
#set xPixels 400; set yPixels 400;
#recorder plot $fileName $windowTitle $xLoc $yLoc $xPixels $yPixels -columns 2 1
# -------------------------------------------------define load pattern
pattern Plain 1 Linear {
   #   node  FX          FY
   load  2  0.0 1.0 
} 
# Perform the pushover analysis
# characteristics of pushover analysis
set Dincr [expr 0.00001 ];	# displacement increment. you want this to be small, but not too small to slow analysis

# characteristics of cyclic analysis
set iDmax { 0.0005 0.0010 0.0015 0.002 0.0025  0.0030}; # 0.03 0.04 0.05 0.06 0.075 0.1 0.3 0.4 0.5 0.75 1.0 };		# vector of displacement-cycle peaks, in terms of storey drift ratio
#set iDmax { 0.01 0.02 0.03 };
set CycleType Full;				# you can do Full / Push / HalfCycle  with the proc
set Ncycles 1;				# specify the number of cycles at each peak
set fmt1 "%s Cyclic analysis : CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis

# Plain Constraints#Lagrange Multipliers#  Penalty Method#Transformation Method
set constraintsTypeStatic Plain;		# default;
constraints $constraintsTypeStatic
#   Plain#   RCM
set numbererTypeStatic Plain
numberer $numbererTypeStatic
#          ProfileSPD #          BandGeneral#          BandSPD#          SparseGeneral#          SparseSPD -- D
#          UmfPack#          FullGeneral
set systemTypeStatic FullGeneral;		# try UmfPack for large model
system $systemTypeStatic
set Tol 1.e-6;
set testType EnergyIncr
set maxNumIter 6
#set algorithmType Newton
set algorithmType KrylovNewton
set printFlag 0

test $testType $Tol	$maxNumIter $printFlag

#algorithm Newton -initial
algorithm KrylovNewton
integrator DisplacementControl $IDctrlNode $IDctrlDOF 0.0;
set analysisTypeStatic Static
analysis $analysisTypeStatic 

#set ok [analyze 1]

#  ---------------------------------    perform Static Cyclic Displacements Analysis
# ----------- set up analysis parameters

source GeneratePeaks.tcl;
foreach Dmax $iDmax {
	#set iDstep [GeneratePeaks $Dmax $Dincr $CycleType];	# this proc is defined above
	source tmpDsteps.tcl;
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep
			set Dincr [expr $D1 - $D0]
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			analysis Static
			# ----------------------------------------------first analyze command------------------------
			set ok [analyze 1]
			# ----------------------------------------------if convergence failure-------------------------
			# if analysis fails, we try some other stuff
			# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent .."
				test NormDispIncr $Tol 2000 0;
				algorithm Newton -initial;
				set ok [analyze 1];
				test $testType $Tol $maxNumIter 0;
				algorithm $algorithmType;
				puts "Trying Newton with Initial Tangent failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Broyden ..";
				algorithm Broyden 8
				set ok [analyze 1 ]
				algorithm $algorithmType
				puts "Trying Broyden failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Newton With LineSearch ..";
				algorithm NewtonLineSearch 0.8 
				set ok [analyze 1]
				algorithm $algorithmType
				puts "Trying Newton With LineSearch failed to converge...";
			}
			if {$ok != 0} {
				set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] m]
				puts $putout
				return -1
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
			puts "$D0";
			#puts $f1 "[expr [nodeCoord 34 1]+[nodeDisp 34 1]] [expr [nodeCoord 34 2]+[nodeDisp 34 2]]";
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
	puts [format $fmt1 "$CycleType"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] "m error."];
} else {
	puts [format $fmt1 "$CycleType"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] "m done"];
}
remove recorders;
