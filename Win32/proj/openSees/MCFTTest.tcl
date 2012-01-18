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
nDMaterial ElasticIsotropic 1 1.0e4 0.2;
nDMaterial PlaneStress 2 1;
nDMaterial PlaneStrain 3 1;

set k       27777.78; #---bulk modulus
set G       9259.26; #---shear modulus
set sig0    30.0; #---yield stress
set sigInf   60.5; #---final saturation yield stress
set delta    0.1; #---exponential hardening parameter eta>= 0 for rate independent case
set H        0.05; #---linear hardening parameter
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

set k       27777.78; #---bulk modulus
set G       9259.26; #---shear modulus
set sigY    5.0; #---yield stress
set rho     0.398; #---failure surface and associativity
set rhoBar  0.398; 
set Kinf    0.0; #---isotropic hardening
set Ko 	    0.0
set delta1  0.0
set H 	    0.0; #---kinematic hardening
set theta   1.0
set delta2  0.0; #---tension softening
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
#nDMaterial DruckerPrager  $nDMatTag 7.8e-6 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

##nDMaterial PlaneStress $PlaneStressMatTag $nDMatTag

set wfc 30.0;
set wfyv 500.0;
set wfyh1 450.0;
set wE 1.96e5;
set epsc0 0.002;
set fpcu  46.0;
set epsU  -0.04;
set lambda 0.1;
set ft  [expr 0.33*pow($wfc,0.5)];
set Ets [expr $ft/0.002];
set rou1 0;
set rou2 0;
set db1 16;
set db2 18;
set b 0.04
set rouv 0.0152;
set rouh1 0.0267;

uniaxialMaterial MCFTSteel03 11 $wfyv $wE $b
uniaxialMaterial MCFTSteel03 12 $wfyh1 $wE $b

#uniaxialMaterial MCFTConcrete01 14 [expr -$wfc] -0.0022
#uniaxialMaterial MCFTConcrete01 15 [expr -$wfc] -0.0022
uniaxialMaterial MCFTConcrete03 14 [expr -$wfc] -$epsc0 $fpcu $epsU $lambda $ft $Ets
uniaxialMaterial MCFTConcrete03 15 [expr -$wfc] -$epsc0 $fpcu $epsU $lambda $ft $Ets

nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.002 10 75 75.
#nDMaterial DSFMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.002 10 75 75.

nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

#set np 2; # int. points
##set C 0.4; # center of rotation 
#
#section Timoshenko $testSecTag { ;
#	layer2d straight $testMatTag 5 5376.33 -116 0 116 0
#	#fiber2d -100 0 5376 $testMatTag;
#	#fiber2d  -50 0 5376 $testMatTag;
#	#fiber2d    0 0 5376 $testMatTag;
#	#fiber2d   50 0 5376 $testMatTag;
#	#fiber2d  100 0 5376 $testMatTag;
#};

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0;
node 1002 1000.0 0.0;
node 1003 1000.0 1000.0
node 1004 0.0 1000.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 ; #1
fix 1002 0 1 ; #1

# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $testSecTag
element quad 2001 1001 1002 1003 1004 100.0  "PlaneStress" $PlaneStressMatTag

# Create recorder
#recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)
recorder Node -file data/Mphi.out -time -node 1003 1004 -dof 1 2 disp;
# Define constant axial load
set axialLoad -100;
pattern Plain 3001 "Constant" {
    #load 1002 $axialLoad 0.0 0.0
    load 1003  $axialLoad 0.0   
    load 1004  $axialLoad 0.0  
}

# Define analysis parameters
integrator LoadControl 0 1 0 0
system SparseGEN;	# Overkill, but may need the pivoting!SparseGeneral -piv
#test EnergyIncr  1.0e-4 10
test NormDispIncr 1.0e-6 250
numberer Plain
constraints Plain
algorithm Newton
analysis Static

# Do one analysis for constant axial load
analyze 1
loadConst 0.0

# Define reference moment
#pattern Plain 3002 "Linear" {
#    load 1003  1.0  0.0 
#    load 1004  1.0  0.0 
#}
pattern Plain 3002 "Linear" {
    load 1003  0.0 1.0
    load 1004  0.0 1.0
}
# Compute curvature increment
set maxK 10.0;
set numIncr 100;
set dK [expr $maxK/$numIncr];

# Use displacement control at node 1002 for section analysis, dof 3
#integrator DisplacementControl 1003 1 $dK 1 $dK $dK
integrator DisplacementControl 1003 1 $dK 1 $dK $dK
# Do the section analysis
set ok [analyze $numIncr]

# ----------------------------------------------if convergence failure-------------------------
set IDctrlNode 1002
set IDctrlDOF 3
set Dmax $maxK
set Dincr $dK
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr  
set maxNumIterStatic 6
set algorithmTypeStatic Newton

global LunitTXT;					# load time-unit text
if {  [info exists LunitTXT] != 1} {set LunitTXT "Length"};		# set blank if it has not been defined previously
	
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Curv=%.4f /%s";	# format for screen/file output of DONE/PROBLEM analysis

if {$ok != 0} {  
	# if analysis fails, we try some other stuff, performance is slower inside this loop
	set Dstep 0.0;
	set ok 0
	while {$Dstep <= 1.0 && $ok == 0} {	
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr $controlDisp/$Dmax]
		set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
		if {$ok != 0} {;				# reduce step size if still fails to converge
			set Nk 4;			# reduce step size
			set DincrReduced [expr $Dincr/$Nk];
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
				set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
				if {$ok != 0} {  
					# if analysis fails, we try some other stuff
					# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
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
				if {$ok != 0} {;				# stop if still fails to converge
					puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					return -1
				}; # end if
			}; # end for
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;	# bring back to original increment
		}; # end if
	};	# end while loop
};      # end if ok !0
# -----------------------------------------------------------------------------------------------------
    
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}