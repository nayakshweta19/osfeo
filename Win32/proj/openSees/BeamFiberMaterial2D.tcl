# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

# MATERIAL parameters -------------------------------------------------------------------

set nDMatTag 1;
# Material "Steel":    matTag    E    v    rho 
#nDMaterial  ElasticIsotropic3D $nDMatTag 2.0e5 0.3 20

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
nDMaterial DruckerPrager $nDMatTag 7.8 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

#nDmaterial for beam fiber 3D
nDMaterial BeamFiber $nDBeamTag $nDMatTag;

set beamSecTag 21;
#section definition 
#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
section WSection2d $beamSecTag $nDBeamTag  300 10 120 20 16 12
	
#section CSMMFiber2d $SecTag -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 21;
#	fiber2d  10 0 1 21;
#};

#section FiberInt 4 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 2;
#	fiber2d  10 0 1 2;
#};

# --------------------------------------------------------------------------------------------------

# Define two nodes
node 1001 0.0 0.0;
node 1002 0.0 1000.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 1
#fix 1002 0 0 0

set transfTag 1;
geomTransf Linear $transfTag
set np 5; # int. points
# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $secTag
#forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag
element nonlinearBeamColumn 2001 1001 1002 $np $beamSecTag $transfTag;
	
# Create recorder
recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)
recorder Node -file data/Reaction.out -time -node 1001 -dof 1 2 3 reaction;

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -4750 -7800 8000
  #vrp 0 -500 250
  vup 0 1 0
  #vpn -1 -1 0.5
  viewWindow -400 400 -400 400
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 1

## Define constant axial load
#set axialLoad 1000;
#pattern Plain 3001 "Constant" {;
#	load 1002 0.0 0.0 $axialLoad 0.0 0.0 0.0;
#};
#
## Define analysis parameters
integrator LoadControl 0 1 0 0;
system UmfPack;	# Overkill, but may need the pivoting!
test EnergyIncr  1.0e-9 10;
numberer Plain;
constraints Plain;
algorithm Newton;
analysis Static;
#analyze 1
## Do one analysis for constant axial load
# Compute curvature increment
set maxK 1.0;
set numIncr 100;
set dK [expr $maxK/$numIncr];

set IDctrlNode 1002;
set IDctrlDOF 3;
set Dmax $maxK;
set Dincr $dK;
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr;
set maxNumIterStatic 6;
set algorithmTypeStatic Newton;

# Define reference moment
pattern Plain 3002 "Linear" {;
	load $IDctrlNode 0.0 0.0 1.0;
};

# Use displacement control at node 1002 for section analysis, dof 3
integrator DisplacementControl $IDctrlNode $IDctrlDOF $dK 1 $dK $dK

# Do the section analysis
set ok [analyze $numIncr]

# ----------------------------------------------if convergence failure-------------------------

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
