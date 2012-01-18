# N_FSW13_3D.tcl
wipe;
# Apply axial load first
 
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
 
# ------------------------------
# Start of model generation
# ------------------------------
 
#4.1
# ------------------------------------------------------------------------------
# Create ModelBuilder for beams and columns (with two-dimensions and 3 DOF/node)
# ------------------------------------------------------------------------------
  model basic -ndm 3 -ndf 6
 
# ---------------------------------------------------------
# Create nodes needed for beams and columns (3dof)
# ---------------------------------------------------------
 
# set dimension of the wall and mesh
 
set L 1068;
set H 1068;
set deltaL 356;
set deltaH 356;
set nL [expr $L/$deltaL];
set nH [expr $H/$deltaH];
set t 76.2;
 
set beamType nonlinearBeamColumn;  # dispBeamColumn

# Create nodes for beams and columns
# tag X Y
 
 
set nodeStartID [expr ($nH+1)*($nL+1)+1];   
set eleStartID  [expr $nH*$nL+1]
 
 
set j 0;
 
while {$j < [expr $nH+1]} {
    if { $j < [expr $nH] } {
        node [expr $nodeStartID+$j*2]    0.0       [expr $deltaH*$j]    0.0
        node [expr $nodeStartID+$j*2+1]  [expr $L] [expr $deltaH*$j]    0.0
    } else {
        set i 0;
        while {$i < [expr $nL+1]} {
            node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaL] [expr $H]   0.0
            set i [expr $i+1]
        }  
    } 
    set j [expr $j+1]   
}
 
  node 27 356 0.0 0.0
  node 28 712 0.0 0.0
 
# --------------------------------------------------
# Define materials for nonlinear columns and beams
# --------------------------------------------------
 
# CONCRETE tag f'c ec0 f'cu ecu
# Core concrete (confined)
 
  uniaxialMaterial Concrete01 1 -64.7 -0.0024 -13.0 -0.006
 
# Cover concrete (unconfined)
 
  uniaxialMaterial Concrete01 2 -57.0 -0.002 -0.0 -0.005
 
# STEEL
# Reinforcing steel 
 
set fy 370.0; # Yield stress for #7 bar
set E 216082.0; # Young's modulus
 
# tag fy E0 b
 
    uniaxialMaterial Steel01 3 $fy $E 0.023
 #  uniaxialMaterial SteelZ01  3  $fy $E  49.75 0.033
 
# -----------------------------------------------
# Define cross-section for nonlinear columns
# -----------------------------------------------
# set some paramaters
 
  set colWidth 152.4
  set colDepth 152.4
  set cover 20.0
set As 126.7; # area of no. 4 bars
 
# some variables derived from the parameters
 
set cy1 [expr $colDepth/2.0]
set cz1 [expr $colWidth/2.0]
 
section Fiber 1 {
 
  # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover]
 
  # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1
    patch rect 2 10 1 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1]
    patch rect 2 2 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover]
    patch rect 2 2 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover]
 
  # Create the reinforcing fibers (4 layers)
    layer straight 3 3 $As [expr $cy1-$cover] [expr $cz1-$cover] [expr $cy1-$cover] [expr $cover-$cz1]
    layer straight 3 3 $As [expr $cover-$cy1] [expr $cz1-$cover] [expr $cover-$cy1] [expr $cover-$cz1]
} 
 
 
# -----------------------------------------------
# Define cross-section for nonlinear beams
# -----------------------------------------------
 
# set some paramaters
 
  set beamWidth 152.4
  set beamDepth 152.4
 
# some variables derived from the parameters
 
set by1 [expr $beamDepth/2.0]
 
set bz1 [expr $beamWidth/2.0]
 
section Fiber 2 {
 
  # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$by1] [expr $cover-$bz1] [expr $by1-$cover] [expr $bz1-$cover]
 
  # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1 [expr -$by1] [expr $bz1-$cover] $by1 $bz1
    patch rect 2 10 1 [expr -$by1] [expr -$bz1] $by1 [expr $cover-$bz1]
    patch rect 2 2 1 [expr -$by1] [expr $cover-$bz1] [expr $cover-$by1] [expr $bz1-$cover]
    patch rect 2 2 1 [expr $by1-$cover] [expr $cover-$bz1] $by1 [expr $bz1-$cover]
 
  # Create the reinforcing fibers (2 layers)
    layer straight 3 3 $As [expr $by1-$cover] [expr $bz1-$cover] [expr $by1-$cover] [expr $cover-$bz1]  
    layer straight 3 3 $As [expr $cover-$by1] [expr $bz1-$cover] [expr $cover-$by1] [expr $cover-$bz1]
 
} 


# -------------------------------------------------------
# Define column elements
# -------------------------------------------------------

geomTransf Linear 1 0 0 1

set np 3
set iterNum 40
set iterTol 1e-3

set j 0;  
while {$j < [expr $nH]} {
    if {$j < [expr $nH-1]} {
         # define the columns elements  
	element $beamType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
	element $beamType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol
    } else {
	element $beamType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
	element $beamType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol
    }
    set j [expr $j+1]
}
 
 
# ----------------------------------------------------------
# Define beam elments
# ---------------------------------------------------------
 
  geomTransf Linear 2 0 0 1
 
set j [expr $nH];
 
set bA 48000;
set bE 2e15;
set bI 1e10;  
 
set i 0;
 
while {$i < [expr $nL]} {
 
    element $beamType [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 2 2 -iter $iterNum $iterTol
    #   element  elasticBeamColumn [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $bA $bE $bI 2
 
    set i [expr $i+1]
}  
 
element $beamType 19 17 27 $np 2 2 -iter $iterNum $iterTol
element $beamType 20 27 28 $np 2 2 -iter $iterNum $iterTol
element $beamType 21 28 18 $np 2 2 -iter $iterNum $iterTol
 
# 4.2
# -------------------------------------------------------------------------
# Create ModelBuilder for 2D element (with two-dimensions and 2 DOF/node)
# -------------------------------------------------------------------------
 
model basic -ndm 3 -ndf 3
 
 
# Create nodes & add to Domain - command: node nodeId xCrd yCrd
 
set j 0;
while {$j < [expr $nH+1]} {
    set i 0;
    while {$i < [expr $nL+1]} {
	node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH] 0.0
	set i [expr $i+1]
    }
    set j [expr $j+1]
}
 
# fix one end as a pin, the other end as a roller
fix 1 1 1 1
fix 4 0 1 1
 
# tie nodes between beam, column and 2D elements
 
equalDOF 1 17 1 2 3
equalDOF 4 18 1 2 3
equalDOF 5 19 1 2 3
equalDOF 8 20 1 2 3
equalDOF 9 21 1 2 3
equalDOF 12 22 1 2 3
equalDOF 13 23 1 2 3
equalDOF 14 24 1 2 3
equalDOF 15 25 1 2 3
equalDOF 16 26 1 2 3
 
equalDOF 2 27 1 2 3
equalDOF 3 28 1 2 3
 
 
# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 57.0;
set wfy 419.2;
set wE  187544.0;
set rou1 0.0023;
set rou2 0.0023;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
  uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
  uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
 
 
 
# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
 
set pi 3.141592654
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
nDMaterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
 
 
# -----------------------------------------------------------------
#  Define 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
set j 0;
while {$j < $nH} {
    set i 0;
    while {$i < $nL} {
     # Create quad elements - command:
     # element quad eleID node1 node2 node3 node4 thick  type         matID
 
	element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15
 
	set i [expr $i+1]
    }
    set j [expr $j+1]
}
 
#4.3
# ---------------------------
# Define horizontal loads
# ---------------------------
 
 
set N 89000.0;
 
 
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
 
    # Create the nodal load - command: load nodeID xForce yForce
    load 13 0 [expr  -$N]  0
    load 16 0 [expr  -$N]  0
}
 
#4.4
# ------------------------------
# End of model generation
# ------------------------------
 
# ------------------------------
# Start of analysis generation
# ------------------------------
 
 
# Create the system of equation, a sparse solver with partial pivoting
  system BandGeneral
 
 
# Create the constraint handler
  constraints Plain
 
# Create the DOF numberer
  numberer Plain
 
 
# Create the convergence test
  test NormDispIncr 1.0e-3 100 5
 
 
# Create the solution algorithm
  algorithm KrylovNewton
 
# Create the integration scheme, the DisplacementControl scheme
  integrator LoadControl 0.1
 
# Create the analysis object
  analysis Static
 
# initialize in case we need to do an initial stiffness iteration
#  initialize
 
# ------------------------------
# End of analysis generation
# ------------------------------
 
# Create a recorder to monitor nodal displacements
#  recorder Node -file N_FSW13.out -time -node 15 -dof 1 2 3 disp
 
# perform the analysis
   analyze 10
 
 
# Print out the state of nodes, if wanted
  print node 13 14 15 16 26 27 28 18
 
# Print out the state of elements, if wanted
 print ele 4 
 
# 4.5
# Set the gravity loads to be constant & reset the time in the domain
 
loadConst -time 0.0
 
# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------
 
# ----------------------------------------------------
# Start of additional modeling for lateral loads
# ----------------------------------------------------
 
 
# -------------------------------
# Define horizontal reference load
# -------------------------------
 
set P 1000.0;
 
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 2 "Linear" {
 
    # Create the nodal load - command: load nodeID xForce yForce
    load 13 [expr   $P/2] 0 0
    load 16 [expr   $P/2] 0 0
    load 4  [expr  -$P/2] 0 0
}
 
# ------------------------------
# End of model generation
# ------------------------------
 
# ------------------------------
# Start of analysis generation
# ------------------------------
 
 
# Create the system of equation, a sparse solver with partial pivoting
  system BandGeneral
 
# Create the constraint handler
  constraints Plain
 
# Create the DOF numberer
  numberer Plain
 
# Create the convergence test
#  test NormDispIncrVaryIter 0.01 14 5 numStep 100 400 300 400 600 800 1200 1600 2050 2500 2960 3470 1000 10 numIter 100 0 0 0 0 100 0 0 0 0 0 0 0 0 
  test NormDispIncr 0.01 100 5
 
# Create the solution algorithm
  algorithm KrylovNewton
 
# Create the integration scheme, the DisplacementControl scheme
# integrator DisplacementPath 13 1 14 numStep 100 400 300 400 600 800 1200 1600 2050 2500 2960 3470 1000 10 increment 0.01 -0.005 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 
  integrator DisplacementControl 13 1 0.01
 
  analysis Static
 
# initialize in case we need to do an initial stiffness iteration
  initialize
 
# ------------------------------
# End of analysis generation
# ------------------------------
 
# Create a recorder to monitor nodal displacements
  recorder Node -file N_FSW13_n13.out -time -node 13 -dof 1 disp
#  recorder Node -file a.out -time -node 13 -dof 1 disp
 
#  MODIFICATION TO MEET WHAT WAS ORIGINALLY IN THE .TCL FILE IN BOOK
#  test NormDispIncrVaryIter 0.01 14 5 numStep 100 400 300 400 600 800 1200 1600 2050 2500 2960 3470 1000 10 numIter 100 0 0 0 0 100 0 0 0 0 0 0 0 0 
#  integrator DisplacementPath 13 1 14 numStep 100 400 300 400 600 800 1200 1600 2050 2500 2960 3470 1000 10 increment 0.01 -0.005 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 
 
# perform the analysis
# analyze 17380
set numSteps   { 100    400  300   400  600   800 1200  1600  2050  2500 2960  3470 1000    10}
set numIters   { 100      0    0     0    0   100    0     0     0     0  500   500  500   300}
set increments {0.01 -0.005 0.01 -0.01 0.01 -0.01 0.01 -0.01  0.01 -0.01 0.01 -0.01 0.01 -0.01}
 
for {set i 0} {$i<14} {incr i 1} {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set increment [lindex $increments $i]
 
    integrator DisplacementControl 13 1 $increment
    if {$numIter == 0} {
	set numIter 1
    } 
    test NormDispIncr 1e-2 $numIter 5
    analyze $numStep
}
 
# Print out the state of elements, if wanted
#  print ele 1 2 3 4 5 6 7 8 9 
 
# Print out the state of nodes, if wanted
#  print node 27 28 18 25 13 14 15 16