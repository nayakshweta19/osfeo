wipe;
# Input File for Specimen C8C
set modelFile [open "model.ops" "w"]
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
# ------------------------------
# Start of model generation
# ------------------------------
# 4.1
# -----------------------------------------------------------------
# Create ModelBuilder for moment carrying walls of hollow bridge columns (with two-dimensions and 3 DOF/node)
# -----------------------------------------------------------------
model basic -ndm 2 -ndf 3
puts $modelFile "model basic -ndm 2 -ndf 3"
# ------------------------------------------------------------------
# Create nodes needed for moment carrying walls of hollow bridge columns (3dof)
# ------------------------------------------------------------------
# set dimension of the wall and mesh
set L 330;
set H 500;
#set deltaL 330;
#set deltaH 500;
set deltaL 330;
set deltaH 500
set nL [expr $L/$deltaL];
set nH [expr $H/$deltaH];
set t 200;
# Creating nodes for moment carrying walls of hollow bridge columns
# tag X Y
set nodeStartID [expr ($nH+1)*($nL+1)+1];
set eleStartID [expr $nH*$nL+1];
set j 0;
while {$j < [expr $nH+1]} {
if { $j < [expr $nH] } {
node [expr $nodeStartID+$j*2] 0.0 [expr $deltaH*$j]
node [expr $nodeStartID+$j*2+1] [expr $L] [expr $deltaH*$j]
puts $modelFile "node [expr $nodeStartID+$j*2] 0.0 [expr $deltaH*$j]"
puts $modelFile "node [expr $nodeStartID+$j*2+1] [expr $L] [expr $deltaH*$j]"
} else {
set i 0;
while {$i < [expr $nL+1]} {
node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaL] [expr $H]
puts $modelFile "node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaL] [expr $H]"
set i [expr $i+1]
}
}
set j [expr $j+1]
}
# ---------------------------------
# Create nodes for tendons
# ---------------------------------
set tendonnodeStartID [expr ($nH+1)*($nL+1)+$nH*2+($nL+1) +1];
set tendoneleStartID [expr $nH*$nL+$nH*2+$nL+1];
set j 0;
while {$j < [expr $nH+1]} {
node [expr $tendonnodeStartID+$j] [expr $deltaL] [expr $deltaH*$j]
puts $modelFile "node [expr $tendonnodeStartID+$j] [expr $deltaL] [expr $deltaH*$j]"
set j [expr $j+1]
}
# Fix supports at base of columns
# tag DX DY RZ
fix [expr $nodeStartID] 1 1 1
fix [expr $nodeStartID+1] 1 1 1
fix [expr $tendonnodeStartID] 1 1 1
puts $modelFile "fix [expr $nodeStartID] 1 1 1"
puts $modelFile "fix [expr $nodeStartID+1] 1 1 1"
puts $modelFile "fix [expr $tendonnodeStartID] 1 1 1"
# -----------------------------------------------------------
# Define nonlinear materials for moment carrying walls and tendons
# -----------------------------------------------------------
# CONCRETE tag f¡¯c ec0 f¡¯cu ecu
# Cover concrete (unconfined)
uniaxialMaterial Concrete01 1 -30 -0.003 -6 -0.01
# Core concrete (confined)
uniaxialMaterial Concrete01 2 -45 -0.003 -9 -0.03
# STEEL
# Longitudinal Reinforcing steel
set fy 434; # Yield stress for bare bar
set E 200000.0; # Young¡¯s modulus
# tag fy E0 fpc rou
uniaxialMaterial Steel01 4 $fy $E 0.01
# tag fy E0 fpc rou epsi
uniaxialMaterial TendonL01 5 1670 $E 1860 0.001 0.0036
# -----------------------------------------------
# Define cross-section for nonlinear columns
# -----------------------------------------------
# set some parameters
set colWidth 860.0
set colDepth 200.0
set cover 50.0
set As 430;
# some variables derived from the parameters
set cy1 [expr $colDepth/2.0]
set cz1 [expr $colWidth/2.0]
# the section 1 is for stirrup confinement is #4@80
section Fiber 1 {
# Create the concrete core fibers
# mat num num
patch rect 2 10 1 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover]
# Create the concrete cover fibers (top, bottom, left, right)
patch rect 1 10 1 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1
patch rect 1 10 1 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1]
patch rect 1 2 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover]
patch rect 1 2 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover]
# Create the reinforcing fibers (2 layers)
layer straight 4 2 $As [expr $cy1-$cover] [expr $cz1 -$cover] [expr $cy1-$cover] [expr $cover-$cz1]
layer straight 4 2 $As [expr $cover-$cy1] [expr $cz1 -$cover] [expr $cover-$cy1] [expr $cover-$cz1]
}
# -----------------------------------------------
# Define cross-section for top beam
# -----------------------------------------------
# E A I
section Elastic 3 2e5 1e6 1e14
# -----------------------------------------------
# Define cross-section for tendons
# -----------------------------------------------
section Fiber 4 {
# Create the concrete core fibers
# Create the reinforcing fibers (2 layers)
layer straight 5 2 280 -185 0 185 0
layer straight 5 2 280 0 -185 0 185
}
# -------------------------------------------------------
# Define column elements
# -------------------------------------------------------
geomTransf Linear 1
set np 5
set iterNum 10
set iterTol 1e-3
set j 0;
while {$j < [expr $nH]} {
if {$j < [expr $nH-1]} {
element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol"
} else {
element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol"
}
set j [expr $j+1]
}
# -------------------------------------------------------
# Define tendon elements
# -------------------------------------------------------
# geomTransf Linear 1
set i 0;
while {$i < [expr $nH]} {
element nonlinearBeamColumn [expr $tendoneleStartID+$i] [expr $tendonnodeStartID+$i] [expr $tendonnodeStartID+1+$i] $np 4 1 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $tendoneleStartID+$i] [expr $tendonnodeStartID+$i] [expr $tendonnodeStartID+1+$i] $np 4 1 -iter $iterNum $iterTol"
set i [expr $i+1]
}
# -------------------------------------------------------
# Define beam elements
# ------------------------------------------------------
geomTransf Linear 2
set j [expr $nH];
set i 0;
while {$i < [expr $nL]} {
element nonlinearBeamColumn [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 3 2 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 3 2 -iter $iterNum $iterTol"
set i [expr $i+1]
}
# 4.2
# ------------------------------------------------------
# Create ModelBuilder for 2D elements representing shear carrying walls of the column (with two-dimensions and 2 DOF/node)
# ------------------------------------------------------
model basic -ndm 2 -ndf 2
puts $modelFile "model basic -ndm 2 -ndf 2"
# Create nodes & add to Domain - command: node nodeId xCrd yCrd
set j 0;
while {$j < [expr $nH+1]} {
set i 0;
while {$i < [expr $nL+1]} {
node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]
puts $modelFile "node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]"
set i [expr $i+1]
}
set j [expr $j+1]
}
# Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
set i 0;
while {$i < [expr $nL+1]} {
fix [expr $i+1] 1 1
puts $modelFile "fix [expr $i+1] 1 1"
set i [expr $i+1]
}
# tying nodes between moment carrying walls and 2D elements representing shear carrying walls of the column
equalDOF 2 7 1 2
equalDOF 4 8 1 2
equalDOF 4 10 1 2
# -----------------------------------------------------
# Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------
# set fc fy E
set wfc 30.0;
set wfyv 434;
set wfyh1 413;
set wE 200000.0;
set rou1 0.001;
set rou2 0.001;
set rouv 0.01;
set rouh1 0.0246; # #10@70
# UniaxialMaterial: steelZ01
# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteL01 14 [expr -$wfc] -0.003
uniaxialMaterial ConcreteL01 15 [expr -$wfc] -0.003
set pi 3.141592654
# NDMaterial: FAFourSteelPCPlaneStress
# tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0
#nDMaterial FAFourSteelPCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003
nDMaterial FAFourSteelPCPlaneStress 21 2.3 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003
# --------------------------------------------------
# Define 2D ReinforceConcretePlaneStress element
# --------------------------------------------------
set j 0;
while {$j < $nH} {
set i 0;
while {$i < $nL} {
# Create quad elements - command:
# element quad eleID node1 node2 node3 node4 thick type matID
element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 21
puts $modelFile "element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 21"
set i [expr $i+1]
}
set j [expr $j+1]
}
close $modelFile
# 4.3
# -------------------------------------
# Define prestress and gravity loads
# -------------------------------------
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
# Create the nodal load - command: load nodeID xForce yForce
load 2 0 [expr -346e3]
load 4 0 [expr -346e3]
#load 6 0 [expr 346e3]
#load 22 0 [expr -346e3]
#load 23 0 [expr -346e3]
#load 24 0 [expr -346e3]
#load 25 0 [expr -486e3]
#load 26 0 [expr -486e3]
#load 27 0 [expr -486e3]
}
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
print node 2 4 7 8
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
load 2 [expr $P/3] 0
load 4 [expr $P/3] 0
}
# ------------------------------
# End of model generation
# ------------------------------
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
test NormDispIncr 1.0e-3 20 0
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
#integrator DisplacementPath 26 1 45 numStep 10 20 20 20 25 30 30 30 35 40 40 40 50 60 60 30 70 80 80 80 100 120 120 120 140 160 160 160 31 15 240 240 240 140 80 40 80 45 50 50 50 53 4 48 30 increment 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -2.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 5.0 3.0 -1.0 1.0 -1.0 2.0 -4.0 8.0 -4.0 8.0 -8.0 8.0 -8.0 8.0 4.0 -10.0 8.0 
integrator DisplacementControl 2 1 1.0
# Creating the analysis object
analysis Static
# initialize in case we need to do an initial stiffness iteration
initialize
# ------------------------------
# End of analysis generation
# ------------------------------
# Creating a recorder to monitor nodal displacements
recorder Node -file C8C_m1_disp_Output.out -time -node 2 -dof 1 disp
# perform the analysis
analyze 10
# Print out the state of nodes
print node 2 4