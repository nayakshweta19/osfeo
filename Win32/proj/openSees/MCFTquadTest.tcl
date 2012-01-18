# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 2;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

set pi 3.141592654;
set nDMatTag 6;
set PlaneStressMatTag 20;
set testMatTag 101;
set testSecTag 201;
#nDMaterial ElasticIsotropic $matTag $E $v
nDMaterial ElasticIsotropic 1 10e3 0.167;
nDMaterial PlaneStress 2 1;
nDMaterial PlaneStrain 3 1;

#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sig0    30.0
#---final saturation yield stress
        set sigInf   60.5
#---exponential hardening parameter eta>= 0 for rate independent case
        set delta    0.1
#---linear hardening parameter
        set H        0.05
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sigY    5.0
#---failure surface and associativity
	set rho     0.398
	set rhoBar  0.398
#---isotropic hardening
	set Kinf    0.0
	set Ko 	    0.0
	set delta1  0.0
#---kinematic hardening
	set H 	    0.0
	set theta   1.0
#---tension softening
	set delta2  0.0
 
#--material models
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
#nDMaterial DruckerPrager  $nDMatTag 7.8e-6 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

##nDMaterial PlaneStress $PlaneStressMatTag $nDMatTag

# set fc fy E
set wfc 76.0;
set wfyv 500.0;
set wfyh1 450.0;
set wE 1.96e5;
set rou1 0;
set rou2 0;
set rouv 0.0352;
set rouh1 0.0367;

# tag fy E0 fpc rou
#uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
#uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
uniaxialMaterial MCFTSteel01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial MCFTSteel01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
#uniaxialMaterial ConcreteZ01 14 [expr -$wfc] -0.0022
#uniaxialMaterial ConcreteZ01 15 [expr -$wfc] -0.0022
uniaxialMaterial MCFTConcrete01 14 [expr -$wfc] -0.0022
uniaxialMaterial MCFTConcrete01 15 [expr -$wfc] -0.0022

# NDMaterial: FAFourSteelPCPlaneStress
#                                                  tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 19:52 2011-10-30rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelPCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteetPCPlaneStress           matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteetPCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAFourSteelRCPlaneStress  (100)
#                   $PlaneStressMatTag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelRCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 $rouv $rouh1 $wfc 0.0 $wfyv $wE 0.003

#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
#nDMaterial CSMMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $wfc $wfyv $wE -0.002

nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002 10 140 140.

nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

set np 2; # int. points
#set C 0.4; # center of rotation 

section Timoshenko $testSecTag { ;
	layer2d straight $testMatTag 5 5376.33 -116 0 116 0
	#fiber2d -100 0 5376 $testMatTag;
	#fiber2d  -50 0 5376 $testMatTag;
	#fiber2d    0 0 5376 $testMatTag;
	#fiber2d   50 0 5376 $testMatTag;
	#fiber2d  100 0 5376 $testMatTag;
};

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0
node 1002 1.0 0.0
node 1003 1.0 1.0
node 1004 0.0 1.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 ; #1
fix 1002 1 1 ; #1

# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $testSecTag
element quad 2001 1001 1002 1003 1004 1000.0  "PlaneStress"  $PlaneStressMatTag

# Create recorder
#recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)
recorder Node -file data/Mphi.out -time -node 1003 1004 -dof 1 2 disp;
# Define constant axial load
set axialLoad 100;
pattern Plain 3001 "Constant" {
    #load 1002 $axialLoad 0.0 0.0
    load 1003 $axialLoad 0.0 
    load 1003 $axialLoad 0.0 
}

# Define analysis parameters
integrator LoadControl 0 1 0 0
system SparseGEN;	# Overkill, but may need the pivoting!SparseGeneral -piv
test EnergyIncr  1.0e-9 10
numberer Plain
constraints Plain
algorithm Newton
analysis Static

# Do one analysis for constant axial load
analyze 1
loadConst 0.0

# Define reference moment
pattern Plain 3002 "Linear" {
	load 1002 0.0 0.0 1.0 
}
# Compute curvature increment
set maxK 10.0;
set numIncr 100;
set dK [expr $maxK/$numIncr];

# Use displacement control at node 1002 for section analysis, dof 3
integrator DisplacementControl 1002 3 $dK 1 $dK $dK

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