# Input File for Specimen C8C
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
# ------------------------------
# Start of model generation
# ------------------------------
# 4.1
# -----------------------------------------------------------------
# Create ModelBuilder for moment carrying walls of hollow bridge columns (with two-dimensions and ?? DOF/node)
# -----------------------------------------------------------------
wipe;
model basic -ndm 3 -ndf 6
set fileModel [open "smodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  6 \n"
# ------------------------------------------------------------------
# Create nodes needed for moment carrying walls of hollow bridge columns (?? ? dof)
# ------------------------------------------------------------------
# set dimension of the wall and mesh
set L 540;
set B 300;
set H 1200;
set deltaL 270; #Y 
set deltaB 300; #X
set deltaH 300; #Z
set nL [expr $L/$deltaL]; #1
set nB [expr $B/$deltaB]; #1
set nH [expr $H/$deltaH]; #4
set t 120;
#  ^
#  |Y  (L)
#  
#  ---------
#  |       |
#  |       |    nL(k)   
#  |       |
#  ---------     -> X  (B)
#      nB(j)
# Creating nodes for moment carrying walls of hollow bridge columns
# tag X Y Z
set nodeID 1;
for {set i 0} {$i <= [expr $nH] } {incr i 1} {
	set Z [expr $i*$deltaH];
	for {set j 0} {$j <=[expr $nB]} {incr j 1} {
		set X [expr $j*$deltaB];
		set Y 0.;
		node $nodeID $X $Y $Z;
		puts $fileModel "node $nodeID $X $Y $Z";
		incr nodeID;
		set Y [expr $L];
		node $nodeID $X $Y $Z;
		puts $fileModel "node $nodeID $X $Y $Z";
		incr nodeID;
	}
	for {set k 1} {$k <= [expr $nL-1]} {incr k 1} {
		set Y [expr $k*$deltaL];
		set X 0.;
		node $nodeID $X $Y $Z;
		puts $fileModel "node $nodeID $X $Y $Z";
		incr nodeID;
		set X [expr $B];
		node $nodeID $X $Y $Z;
		puts $fileModel "node $nodeID $X $Y $Z";
		incr nodeID;
	}
}
# Fix supports at base of columns
# tag DX DY RZ
fixZ 0.0 1 1 1 1 1 1;
puts $fileModel "fixZ 0.0 1 1 1 1 1 1";

# 4.2
# ------------------------------------------------------
# Create ModelBuilder for 2D elements representing shear carrying walls of the column (with two-dimensions and 2 DOF/node)
# ------------------------------------------------------
#model basic -ndm 2 -ndf 2
# -----------------------------------------------------
# Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------
# set fc fy E
set wfc 57.0;
set wfy 419.2;
set wE  187544.0;
set rou1 0.014;
set rou2 0.035;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
 
# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
 
set pi 3.141592654
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
nDMaterial FAReinforcedConcretePlaneStress 15  6.7e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

# create the material
nDMaterial ElasticIsotropic   1   4.0e4   0.22  6.75e-6;

nDMaterial FAReinforcedConcretePlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

# ---------------------------------------
# Define Section tag
# ---------------------------------------

set elastEleSecTag 2;
set thick 120;
section  PlateFiber  $elastEleSecTag   1  $thick

set eleSecTag 3;
section  PlateFiber  $eleSecTag   2  $thick  

# --------------------------------------------------
# Define 2D ReinforceConcretePlaneStress element
# --------------------------------------------------
set eleID 1;
set N [expr ($nB+$nL)*2];
for {set i 0} {$i < [expr $nH-1]} {incr i 1} {
	for {set j 0} {$j < [expr $nB]} {incr j 1} {; # along X axial one side by one side
		element shell02 $eleID [expr $i*$N+$j*2+1] [expr $i*$N+$j*2+3] [expr ($i+1)*$N+$j*2+3] [expr ($i+1)*$N+$j*2+1] $eleSecTag; #$t PlaneStress 15
		puts $fileModel "element quad $eleID [expr $i*$N+$j*2+1] [expr $i*$N+$j*2+3] [expr ($i+1)*$N+$j*2+3] [expr ($i+1)*$N+$j*2+1] $t PlaneStress 15"
		incr eleID 1;
		element shell02 $eleID [expr $i*$N+$j*2+2] [expr $i*$N+$j*2+4] [expr ($i+1)*$N+$j*2+4] [expr ($i+1)*$N+$j*2+2] $eleSecTag; #$t PlaneStress 15
		puts $fileModel "element quad $eleID [expr $i*$N+$j*2+2] [expr $i*$N+$j*2+4] [expr ($i+1)*$N+$j*2+4] [expr ($i+1)*$N+$j*2+2] $t PlaneStress 15"
		incr eleID 1;
	}
	# along Y axial one side by one side
	# left
	element shell02 $eleID [expr $i*$N+2*($nB+1)+1] [expr $i*$N+1] [expr ($i+1)*$N+1] [expr ($i+1)*$N+2*($nB+1)+1] $eleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+1] [expr $i*$N+1] [expr ($i+1)*$N+1] [expr ($i+1)*$N+2*($nB+1)+1] $t PlaneStress 15"
	incr eleID 1;
	element shell02 $eleID [expr $i*$N+2*($nB+1)+2] [expr $i*$N+$nB*2+1] [expr ($i+1)*$N+$nB*2+1] [expr ($i+1)*$N+2*($nB+1)+2] $eleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+2] [expr $i*$N+$nB*2+1] [expr ($i+1)*$N+$nB*2+1] [expr ($i+1)*$N+2*($nB+1)+2] $t PlaneStress 15"
	incr eleID 1;
	#middle
	for {set k 1} {$k < [expr $nL-1]} {incr k 1} {
	    element shell02 $eleID [expr $i*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+$k*2+1] [expr $i*$N+2*($nB+1)+$k*2+1] $eleSecTag; #$t PlaneStress 15
		puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+$k*2+1] [expr $i*$N+2*($nB+1)+$k*2+1] $t PlaneStress 15"
		incr eleID 1;
		element shell02 $eleID [expr $i*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2+2] [expr $i*$N+2*($nB+1)+$k*2+2] $eleSecTag; #$t PlaneStress 15
		puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2+2] [expr $i*$N+2*($nB+1)+$k*2+2] $t PlaneStress 15"
		incr eleID 1;
	}
	#right
	element shell02 $eleID [expr $i*$N+2] [expr $i*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2] $eleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2] [expr $i*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2] $t PlaneStress 15"
	incr eleID 1;
	element shell02 $eleID [expr $i*$N+2*($nB+1)] [expr $i*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)] $eleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)] [expr $i*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)] $t PlaneStress 15"
	incr eleID 1;
}
set i [expr $nH-1];
for {set j 0} {$j < [expr $nB]} {incr j 1} {; # along X axial one side by one side
	element shell02 $eleID [expr $i*$N+$j*2+1] [expr $i*$N+$j*2+3] [expr ($i+1)*$N+$j*2+3] [expr ($i+1)*$N+$j*2+1] $elastEleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+$j*2+1] [expr $i*$N+$j*2+3] [expr ($i+1)*$N+$j*2+3] [expr ($i+1)*$N+$j*2+1] $t PlaneStress 15"
	incr eleID 1;
	element shell02 $eleID [expr $i*$N+$j*2+2] [expr $i*$N+$j*2+4] [expr ($i+1)*$N+$j*2+4] [expr ($i+1)*$N+$j*2+2] $elastEleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+$j*2+2] [expr $i*$N+$j*2+4] [expr ($i+1)*$N+$j*2+4] [expr ($i+1)*$N+$j*2+2] $t PlaneStress 15"
	incr eleID 1;
}
# along Y axial one side by one side
# left
element shell02 $eleID [expr $i*$N+2*($nB+1)+1] [expr $i*$N+1] [expr ($i+1)*$N+1] [expr ($i+1)*$N+2*($nB+1)+1] $elastEleSecTag; #$t PlaneStress 15
puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+1] [expr $i*$N+1] [expr ($i+1)*$N+1] [expr ($i+1)*$N+2*($nB+1)+1] $t PlaneStress 15"
incr eleID 1;
element shell02 $eleID [expr $i*$N+2*($nB+1)+2] [expr $i*$N+$nB*2+1] [expr ($i+1)*$N+$nB*2+1] [expr ($i+1)*$N+2*($nB+1)+2] $elastEleSecTag; #$t PlaneStress 15
puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+2] [expr $i*$N+$nB*2+1] [expr ($i+1)*$N+$nB*2+1] [expr ($i+1)*$N+2*($nB+1)+2] $t PlaneStress 15"
incr eleID 1;
#middle
for {set k 1} {$k < [expr $nL-1]} {incr k 1} {
    element shell02 $eleID [expr $i*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+$k*2+1] [expr $i*$N+2*($nB+1)+$k*2+1] $elastEleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+($k-1)*2+1] [expr ($i+1)*$N+2*($nB+1)+$k*2+1] [expr $i*$N+2*($nB+1)+$k*2+1] $t PlaneStress 15"
	incr eleID 1;
	element shell02 $eleID [expr $i*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2+2] [expr $i*$N+2*($nB+1)+$k*2+2] $elastEleSecTag; #$t PlaneStress 15
	puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2] [expr ($i+1)*$N+2*($nB+1)+$k*2+2] [expr $i*$N+2*($nB+1)+$k*2+2] $t PlaneStress 15"
	incr eleID 1;
}
#right
element shell02 $eleID [expr $i*$N+2] [expr $i*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2] $elastEleSecTag; #$t PlaneStress 15
puts $fileModel "element quad $eleID [expr $i*$N+2] [expr $i*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2-1] [expr ($i+1)*$N+2] $t PlaneStress 15"
incr eleID 1;
element shell02 $eleID [expr $i*$N+2*($nB+1)] [expr $i*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)] $elastEleSecTag; #$t PlaneStress 15
puts $fileModel "element quad $eleID [expr $i*$N+2*($nB+1)] [expr $i*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)+($nL-1)*2] [expr ($i+1)*$N+2*($nB+1)] $t PlaneStress 15"
incr eleID 1;

# 4.3
# -------------------------------------
# Define prestress and gravity loads
# -------------------------------------
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
	# Create the nodal load - command: load nodeID xForce yForce
	for {set i 0} {$i <[expr $N]} {incr i 1} {
		load [expr $nodeID-$i-1] 0.0 0.0 -4e3 0.0 0.0 0.0;
		puts $fileModel "load [expr $nodeID-$i] 0.0 0.0 -4e3 0.0 0.0 0.0";
	}
}

close $fileModel;
# ------------------------------
# End of model generation
# ------------------------------
# 4.4
# ------------------------------
# Start of analysis generation
# ------------------------------
# Creating the system of equation, a sparse solver with partial pivoting
system BandGeneral
# Creating the constraint handler
constraints Plain
# Creating the DOF numberer
numberer Plain
# Creating the convergence test
test NormDispIncr 1.0e-3 100 5
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
integrator LoadControl 0.1
# Creating the analysis object
analysis Static
# ------------------------------
# End of analysis generation
# ------------------------------
# perform the analysis
analyze 10
# Print out the state of nodes
#print node 25 26 27 44 45 46
# 4.5
# Set the gravity and prestress loads to be constant & reset the time in the domain
loadConst -time 0.0
# ------------------------------------------------------
# End of Model Generation & Initial Gravity and Prestress Load Analysis
# ------------------------------------------------------
# ----------------------------------------------------
# Start of additional modeling for lateral loads
# ----------------------------------------------------
# ---------------------------
# Define horizontal loads
# ---------------------------
set P 1000.0;
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 2 "Linear" {
	# Create the nodal load - command: load nodeID xForce yForce
	for {set i 0} {$i <[expr $N]} {incr i 1} {
		load [expr $nodeID-$i-1] 1.0 0.0 0.0 0.0 0.0 0.0;
	}
}
# ------------------------------
# End of model generation
# ------------------------------
# create the display
recorder display nodeNum 10 10 800 600 -wipe
# next three commmands define viewing system, all values in global coords
prp -800 -100 400;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 0 1;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -800 800 -800 800; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 -1 10; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# 4.6
# ------------------------------
# Start of analysis generation
# ------------------------------
# Creating the system of equation, a sparse solver with partial pivoting
system BandGeneral
# Creating the constraint handler
constraints Plain
# Creating the DOF numberer
numberer Plain
# Creating the convergence test
#test NormDispIncrVaryIter 0.001 45 5 numStep 10 20 20 20 25 30 30 30 35 40 40 40 50 60 60 30 70 80 80 80 100 120 120 120 140 160 160 160 31 15 240 240 240 140 80 40 80 45 50 50 50 53 4 48 30 -numIter 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
test NormDispIncr 1.0e-3 100 5
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
#integrator DisplacementPath 26 1 45 numStep 10 20 20 20 25 30 30 30 35 40 40 40 50 60 60 30 70 80 80 80 100 120 120 120 140 160 160 160 31 15 240 240 240 140 80 40 80 45 50 50 50 53 4 48 30 increment 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -2.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 5.0 3.0 -1.0 1.0 -1.0 2.0 -4.0 8.0 -4.0 8.0 -8.0 8.0 -8.0 8.0 4.0 -10.0 8.0 
integrator DisplacementControl [expr $nodeID-1] 1 0.05
# Creating the analysis object
analysis Static
# initialize in case we need to do an initial stiffness iteration
initialize
# ------------------------------
# End of analysis generation
# ------------------------------
set IDctrlNode [expr $nodeID-1];
set IDctrlDOF 1;
set Dincr 0.1;
# perform the analysis    /////////////////////////////////////////////////////////////
#                0.08   -0.08  0.12   -0.12   0.3   -0.3   0.4   -0.4   0.6   -0.6   1.1   -1.1
#set numSteps     {   508   1016  1270    1524  2667   3810  4445   5080  6350   7620 10795  13970 }
set numIters     {     0    300   400     400   400    400   400    400   400    400   400    400 }
#set increments   { 0.004 -0.004 0.004  -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 }
set numSteps     { 20   40   50   60  106  152 178  204 254  304 432  560 }
set increments   { 0.1 -0.1 0.1  -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 }
	
for { set i 0 } { $i<12 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set Dincr [lindex $increments $i]
    integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
    if {$numIter == 0} {  set numIter 100 };
    test NormDispIncr 1.0e-3 $numIter 5
    analyze $numStep

    set j 0;
    while { $j < $numStep } {
      set ok [analyze 1]
      # ----------------------------------------------if convergence failure-------------------------
      # if analysis fails, we try some other stuff
      # performance is slower inside this loop  global maxNumIterStatic;      # max no. of iterations performed before "failure to converge" is ret'd
      if {$ok != 0} {
        puts "Trying Newton with Initial Tangent .."
        algorithm Newton -initial
        set ok [analyze 1]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying Broyden .."
        algorithm Broyden 8
        set ok [analyze 1 ]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying NewtonWithLineSearch .."
        algorithm NewtonLineSearch 0.8 
        set ok [analyze 1]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        set putout "PROBLEM: $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF]";
        puts $putout
        set ok [analyze 1]
      }; # end if
      incr j 1;
      #record;
      #puts $end
    }
    puts $i
}
#remove recorders;
