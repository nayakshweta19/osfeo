# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

# MATERIAL parameters -------------------------------------------------------------------
set SecTagFlex 2;			# assign a tag number to the column flexural behavior
set SecTagAxial 3;			# assign a tag number to the column axial behavior	
set SecTag 1;			# assign a tag number to the column section tag

# COLUMN section
# calculated stiffness parameters
set EASec $Ubig;				# assign large value to axial stiffness
set MySec [expr 130000*$kip*$in];		# yield moment
set PhiYSec [expr 0.65e-4/$in];		# yield curvature
set EICrack [expr $MySec/$PhiYSec];		# cracked section inertia
set b 0.01 ;				# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $SecTagFlex $MySec $EICrack $b; 		# bilinear behavior for flexural moment-curvature
uniaxialMaterial Elastic $SecTagAxial $EASec;			# this is not used as a material, this is an axial-force-strain response
section Aggregator $SecTag $SecTagAxial P $SecTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)

set nDMatTag 6;
# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D $nDMatTag 2e4 0.167  0 

set nDBeamTag 10;
#nDmaterial for beam fiber 3D
nDMaterial BeamFiber $nDBeamTag $nDMatTag;

set beamSecTag 21;
#section definition 
#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
section WSection2d $beamSecTag $nDBeamTag  300 10 120 10 10 6
#section CSMMFiber2d $SecTag -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 21;
#	fiber2d  10 0 1 21;
#};
#section FiberInt 4 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 2;
#	fiber2d  10 0 1 2;
#};

# --------------------------------------------------------------------------------------------------

# Define two nodes at (0,0)
node 1001 0.0 0.0 0.0
node 1002 0.0 0.0 1000.0

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 1 1 1 1
#fix 1002 0 0 0

set transfTag 1;
geomTransf Linear $transfTag 0 1 0
set np 5; # int. points
# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $secTag
#forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag
element nonlinearBeamColumn 2001 1001 1002 $np $beamSecTag $transfTag;
	
# Create recorder
recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 4 5 6 disp;	# output moment (col 1) & curvature (col 2)
recorder Node -file data/Reaction.out -time -node 1001 -dof 1 2 3 4 5 6 reaction;

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -1750 -5800 6000
  #vrp 0 -500 250
  vup 0 0 1
  #vpn -1 -1 0.5
  viewWindow -400 400 -400 400
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 1
	
## Define constant axial load
#set axialLoad 1000;
#pattern Plain 3001 "Constant" {
#	load 1002 0.0 0.0 $axialLoad 0.0 0.0 0.0
#}
#
## Define analysis parameters
integrator LoadControl 0 1 0 0
system SparseGeneral -piv;	# Overkill, but may need the pivoting!
test EnergyIncr  1.0e-9 10
numberer Plain
constraints Plain
algorithm Newton
analysis Static

## Do one analysis for constant axial load
#analyze 1
set IDctrlNode 1002
set IDctrlDOF 4
set Dmax $maxK
set Dincr $dK
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr  
set maxNumIterStatic 6
set algorithmTypeStatic Newton

# Define reference moment
pattern Plain 3002 "Linear" {
	load $IDctrlNode 0.0 0.0 0.0 1.0 0.0 0.0
}

# Compute curvature increment
set maxK 1.
set numIncr 100;
set dK [expr $maxK/$numIncr]

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
