# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

# MATERIAL parameters -------------------------------------------------------------------

set nDMatTag 1;
# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D $nDMatTag 5.e5 0.3 0.0
#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sig0    3.0
#---final saturation yield stress
        set sigInf   6.5
#---exponential hardening parameter eta>= 0 for rate independent case
        set delta    0.1
#---linear hardening parameter
        set H        0.05
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

set nDBeamTag 11;

#--define material parameters for the model
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
#nDMaterial DruckerPrager $nDMatTag 7.8 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

#nDmaterial for beam fiber 3D
#nDMaterial BeamFiber $nDBeamTag $nDMatTag;
#                Dodd_Restrepo $tag $Fy $Fsu $ESH $ESU $Youngs $ESHI $FSHI <$OmegaFac $Conv>
#uniaxialMaterial Dodd_Restrepo $nDBeamTag 200. 275. 0.002 0.015 2.e5 0.005 245.
#                         $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
#uniaxialMaterial Steel01 $nDBeamTag 500 0.002 0.03;
	
set beamSecTag 1;
#section definition 
#section Fiber $beamSecTag { ;
#	fiber -10 0 1 $nDBeamTag;
#	fiber  10 0 1 $nDBeamTag;
#};
set uDBeamTag 2;
nDMaterial BeamFiber2d $uDBeamTag $nDMatTag;
#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
#section WSection2d $beamSecTag $nDBeamTag  300 10 120 20 16 12
section Timoshenko $beamSecTag { ;
	fiber2d -10 0 100 $uDBeamTag;
	fiber2d  -5 0 100 $uDBeamTag;
	fiber2d   0 0 100 $uDBeamTag;
	fiber2d   5 0 100 $uDBeamTag;
	fiber2d  10 0 100 $uDBeamTag;
};

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0
node 1002 0.0 0.0

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 1
fix 1002 0 0 0

# Define element
#                         tag ndI ndJ  secTag
element zeroLengthSection  2001   1001   1002  $beamSecTag

# Create recorder
recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)

# Define constant axial load
set axialLoad 100;
pattern Plain 3001 "Constant" {
	load 1002 $axialLoad 0.0 0.0
}

# Define analysis parameters
integrator LoadControl 0 1 0 0
system SparseGeneral -piv;	# Overkill, but may need the pivoting!
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