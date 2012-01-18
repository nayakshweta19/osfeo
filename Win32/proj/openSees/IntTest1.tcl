wipe;
model basic -ndm 2 -ndf 3 
set tolAx 1.0e-3;
set iterAx 100;
set tolLatNew 1.0e-6;
set iterLatNew 100;
set tolLatIni 1.0e-5;
set iterLatIni 1000; 
set dUi 5.0; # Displacement increment
set maxU 1600; # Max. Displacement

#concrete parameters
set fc -25.0;
set eo -0.003;
set r 1.88; 
set k 1.0;
set alphaC 0.32;
set fcr 2.14;
set ecr 0.00008;
set b 0.4;
set alphaT 0.08;
set fcc -18.9;
set ecc -0.0038;
set rcc 1.88; 

#steel parameters
set fyRed 345.0;
set by 0.02;
set Ry 1.564;
set roy 0.00327;
set fyBRed 325.0;
set byB 0.02; 
set RyB 25.0;
set royB 0.0293;
set fxRed 325.0;
set bx 0.02;
set Rx 1.564;
set rox 0.00327;

#concrete (confined and unconfined)
#uniaxialMaterial Concrete06 1 $fcc $ecc $rcc $k $alphaC $fcr $ecr $b $alphaT;
#uniaxialMaterial Concrete06 2 $fc $eo $r $k $alphaC $fcr $ecr $b $alphaT;
uniaxialMaterial Concrete01 1 $fcc  $ecc [expr 0.2*$fcc] -0.008
uniaxialMaterial Concrete01 2 $fc $eo [expr 0.2*$fcc] -0.008

# steel (bound., web and horiz. reinforcement)
set E 2.0e5; 
uniaxialMaterial Steel02 1003 $fyBRed $E $byB $RyB 0.9 0.1 0 0.1 0 0.1 ; 
uniaxialMaterial Steel02 1004 $fyRed $E $by $Ry 0.9 0.1 0 0.1 0 0.1 ; 
uniaxialMaterial Steel02 1005 $fxRed $E $bx $Rx 0.9 0.1 0 0.1 0 0.1 ;

# Define cross-section 
set t1 360;
set NStrip1 1; # thickness 1
set t2 120;
set NStrip2 4; # thickness 2
set t3 360;
set NStrip3 1; # thickness 3
geomTransf LinearInt 1
set np 4; # int. points
set C 0.4; # center of rotation 

#section definition 
section FiberInt 2 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { 
#vertical fibers: fiber yLoc zLoc area matTag
fiber -240    0 43200 1;   fiber -240    0 503   1003;
fiber  -90  135 21600 1;   fiber  -90  135 251.5 1003;
fiber   90  135 21600 1;   fiber   90  135 251.5 1003;
fiber  -90 -135 21600 1;   fiber  -90 -135 251.5 1003;
fiber   90 -135 21600 1;   fiber   90 -135 251.5 1003;
fiber  240    0 43200 1;   fiber  240    0 503   1003;
#horiz. reinf.: Hfiber yLoc zLoc area matTag
Hfiber 0 0 1680 1005;
}

#nodes
set n 2;
node 1 0 0.0;
node 2 0 1200.0;

mass 1 1e6 1e6 0;
mass 2 1e6 1e6 0;

#element definition
#element Timoshenko 1 1 2 $np 2 1;
element dispBeamColumnInt 1 1 2 $np 2 1 $C
#element dispBeamColumnInt 3 4 5 $np 2 1 $C 
#element dispBeamColumnInt 4 5 6 $np 2 1 $C
#element dispBeamColumnInt 5 6 7 $np 2 1 $C
#element dispBeamColumnInt 6 7 8 $np 2 1 $C
#element dispBeamColumnInt 7 8 9 $np 2 1 $C
#element dispBeamColumnInt 8 9 2 $np 2 1 $C

fix 1 1 1 1;
fix 2 0 0 0;

# Set axial load 
pattern Plain 1 Constant {;
	load 2 0.0 -85.0 0.0
}

initialize;
integrator LoadControl 0 1 0 0;
system BandGeneral;
test NormUnbalance $tolAx $iterAx 0;
numberer Plain;
constraints Plain;
algorithm ModifiedNewton -initial;
analysis Static;

# perform the gravity load analysis, 
analyze [expr 1]

## Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0
remove recorders

# set lateral load 
pattern Plain 2 Linear { 
	load 2 1.0 0.0 0.0;
} 
system BandGeneral;
constraints Transformation;
numberer Plain;

# Create a recorder to monitor nodal displacement and element forces
recorder Node -file nodeTop.out -node 2 -dof 1 disp
recorder Element -file elesX.out -time -ele 1 globalForce

# recorder for element1 section1 steel stress/strain and section force-def.
#recorder Element -file Sect_FandD.out -ele 1 section 1 forceAndDeformation
#recorder Element -file Sect_eX.out -ele 1 section 1 eX
#recorder Element -file Sect_eY.out -ele 1 section 1 eY
#recorder Element -file Sect_sX.out -ele 1 section 1 sX
#recorder Element -file Sect_sY.out -ele 1 section 1 sY
set dU $dUi;
set numSteps [expr int($maxU/$dU)];
integrator DisplacementControl 2 1 $dU 1 $dU $dU
test NormDispIncr $tolLatNew $iterLatNew; 
algorithm KrylovNewton
analysis Static

# characteristics of pushover analysis
set Dincr [expr 1.0 ];	# displacement increment. you want this to be small, but not too small to slow analysis
set IDctrlNode 2;
set IDctrlDOF 1;
source LibAnalysisStaticParameters.tcl;
# initialize in case we need to do an initial stiffness iteration
initialize

# perform the analysis
#analyze 10

# Perform the Cycle pushover analysis
# characteristics of cyclic analysis	
set iDmax { 20 25 30 35 40 45 50 60 70 80 100 };  #}120 140 160 240  ; ## vector of displacement-cycle peaks, in terms of storey drift ratio
set CycleType Full;				# you can do Full / Push / Half cycles with the proc
set Ncycles 1;				# specify the number of cycles at each peak
#  ---------------------------------    perform Static Cyclic Displacements Analysis
source GeneratePeaks.tcl

# -- STATIC PUSHOVER/CYCLIC ANALYSIS

foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep            
			puts "D1=$D1"
			set Dincr [expr $D1 - $D0]
			#puts "Dincr = $Dincr"
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			puts "Displacement Node: $IDctrlNode. Dof: $IDctrlDOF. Incr: $Dincr"
			#integrator ArcLength 1.0 0.1
			# ------------------first analyze command------------------------
			set ok [analyze 1]
			# -----------------if convergence failure-------------------------
			# max no. of iterations performed before "failure to converge" is ret'd
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent .."
				test NormDispIncr $Tol 2000 0;
				algorithm Newton -initial;
				set ok [analyze 1];
				test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;
				algorithm $algorithmTypeStatic;
				puts "Trying Newton with Initial Tangent failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Broyden ..";
				algorithm Broyden 40
				set ok [analyze 1 ]
				algorithm $algorithmTypeStatic
				puts "Trying Broyden failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Newton With LineSearch ..";
				algorithm NewtonLineSearch 0.8 
				set ok [analyze 1]
				algorithm $algorithmTypeStatic
				puts "Trying Newton With LineSearch failed to converge...";
			}
			#if {$ok != 0} {
			#	puts "Trying ArcLength ..";
			#	algorithm $algorithmTypeStatic 
			#	integrator ArcLength 1.0 0.1
			#	set ok [analyze 1]
			#	puts "Trying ArcLength With KrylovNewton failed to converge...";
			#}
			if {$ok != 0} {
				set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";
				set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] m]
				puts $putout
				return -1
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
			puts "$D0";
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl