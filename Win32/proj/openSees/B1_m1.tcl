# Input file for beam B1
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
set modelFile [open "model.ops" "w"]
# ------------------------------
# Start of model generation
# ------------------------------
# ---------------------------------------------------------
# Create ModelBuilder for top and bottom flanges (with twodimensions and 3 DOF/node)
# ---------------------------------------------------------
model basic -ndm 2 -ndf 3
puts $modelFile "model basic -ndm 2 -ndf 3"
# ---------------------------------------------------------
# Create nodes needed for top and bottom flanges (3dof)
# ---------------------------------------------------------
# set dimension of the beam and mesh
set L 457; #24 feet
set H 546; #21.5 inches
set deltaL 457;
set deltaH 546;
set nL [expr $L/$deltaL]; #16 elements along length
set nH [expr $H/$deltaH]; # 1 element along depth of web
set topH 546;
set botH 0;
# Create nodes for top and bottom flanges
# tag X Y
set nodeStartID [expr ($nH+1)*($nL+1)+1]; # 35
set eleStartID [expr $nH*$nL+1]; # 17
set j 0;
while {$j < [expr $nL+1]} {
node [expr $nodeStartID+$j] [expr $j*$deltaL] [expr $botH]
node [expr $nodeStartID+$nL+1+$j] [expr $j*$deltaL] [expr $topH]
puts $modelFile "node [expr $nodeStartID+$j] [expr $j*$deltaL] [expr $botH]"
puts $modelFile "node [expr $nodeStartID+$nL+1+$j] [expr $j*$deltaL] [expr $topH]"
set j [expr $j+1]
}
# Provide roller and hinge supports at two ends of beams
# tag DX DY RZ
fix [expr $nodeStartID] 1 1 0
fix [expr $nodeStartID+$nL] 0 1 0
puts $modelFile "fix [expr $nodeStartID] 1 1 0"
puts $modelFile "fix [expr $nodeStartID+$nL] 0 1 0"
# --------------------------------------------------
# Define nonlinear materials for flanges
# --------------------------------------------------
# CONCRETE tag f¡¯c ec0 f¡¯cu ecu
uniaxialMaterial Concrete01 1 -72.4 -0.0024 -14.0 -0.008
# STEEL
# Reinforcing steel
set E 200000.0; # Young¡¯s modulus
# tag fy E0 b
uniaxialMaterial Steel01 3 413.7 $E 0.001
# for strands fu=270ksi=1861.65MPa
# for strands fy=0.7fu = 1303.0 MPa
uniaxialMaterial Steel01 4 1303.0 $E 0.001
# -----------------------------------------------
# Define cross-section for top and bottom flanges
# -----------------------------------------------
# set some paramaters
# #5 bar
set As1 200
# 1/2 inch seven-wire strands
set As2 99
# for top flange
section Fiber 1 {
patch rect 1 10 4 -152 -70 152 70
# Creating the reinforcing fibers
layer straight 3 2 $As1 -127.0 40.0 127.0 40.0
}
# for bottom flange
section Fiber 2 {
patch rect 1 10 4 -203 -95 203 95
# Creating the prestressing strands
layer straight 4 2 $As2 -153 -49 153 -49
}
# ---------------------------------------------------
# Define flange elements
# ---------------------------------------------------
geomTransf Linear 2
set np 3;
set iterNum 10;
set iterTol 1e-2;
set i 0;
while {$i < [expr $nL]} {
# bottom flange elements
element nonlinearBeamColumn [expr $eleStartID+$i] [expr $nodeStartID+$i] [expr $nodeStartID+1+$i] $np 2 2 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$i] [expr $nodeStartID+$i] [expr $nodeStartID+1+$i] $np 2 2 -iter $iterNum $iterTol"
# top flange elements
element nonlinearBeamColumn [expr $eleStartID+$nL+$i] [expr $nodeStartID+$nL+1+$i] [expr $nodeStartID+$nL+2+$i] $np 1 2 -iter $iterNum $iterTol
puts $modelFile "element nonlinearBeamColumn [expr $eleStartID+$nL+$i] [expr $nodeStartID+$nL+1+$i] [expr $nodeStartID+$nL+2+$i] $np 1 2 -iter $iterNum $iterTol"
set i [expr $i+1]
}
# ------------------------------------------------------
# Create ModelBuilder for 2D web elements (with twodimensions and 2 DOF/node)
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
# Provide roller and hinge supports at two ends of beams
fix 1 1 1
fix 17 0 1
puts $modelFile "fix 1 1 1"
puts $modelFile "fix 17 0 1"
# tie nodes between flange elements and web elements
# tying nodes along bottom flange
equalDOF 2 36 1 2
equalDOF 3 37 1 2
equalDOF 4 38 1 2
equalDOF 5 39 1 2
equalDOF 6 40 1 2
equalDOF 7 41 1 2
equalDOF 8 42 1 2
equalDOF 9 43 1 2
equalDOF 10 44 1 2
equalDOF 11 45 1 2
equalDOF 12 46 1 2
equalDOF 13 47 1 2
equalDOF 14 48 1 2
equalDOF 15 49 1 2
equalDOF 16 50 1 2
# tying nodes along top flange
equalDOF 18 52 1 2
equalDOF 19 53 1 2
equalDOF 20 54 1 2
equalDOF 21 55 1 2
equalDOF 22 56 1 2
equalDOF 23 57 1 2
equalDOF 24 58 1 2
equalDOF 25 59 1 2
equalDOF 26 60 1 2
equalDOF 27 61 1 2
equalDOF 28 62 1 2
equalDOF 29 63 1 2
equalDOF 30 64 1 2
equalDOF 31 65 1 2
equalDOF 32 66 1 2
equalDOF 33 67 1 2
equalDOF 34 68 1 2
# ------------------------------------------------------
# Define materials for 2D PrestressConcretePlaneStress element
# ------------------------------------------------------
# set fc fy E
set wfc 72.4;
set wfpu 1862;
set wfy 413.7;
set wE 200000.0;
set rou1 0.0055;
set rou2 0.00164;
set ec 0.002;
set t 152.4;
# UniaxialMaterial: steelZ01
# tag fy E0 fpu rou epsi
uniaxialMaterial TendonL01 11 [expr 0.7*$wfpu] $wE $wfpu $rou1 0.006
# tag fy E0 fpu rou
uniaxialMaterial SteelZ01 12 $wfy $wE $wfc $rou2
# UniaxialMaterial: concreteL01
# ConcreteL01 tag f¡¯c ec0
uniaxialMaterial ConcreteL01 13 [expr -$wfc] [expr -$ec]
uniaxialMaterial ConcreteL01 14 [expr -$wfc] [expr -$ec]
set pi 3.141592654
# NDMaterial: FAPrestressConcretePlaneStress
# tag rho s1 s2 c1 c2 angle1 angle2 rou1 rou2 fpc fpy fy E0 ec
nDMaterial FAPrestressedConcretePlaneStress 15 0.0 11 12 13 14 [expr 1.0*$pi] [expr 0.5*$pi] $rou1 $rou2 -0.0 $wfc [expr 0.7*$wfpu] $wfy $wE $ec
nDMaterial FAPrestressedConcretePlaneStress 16 0.0 11 12 13 14 [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 -0.0 $wfc [expr 0.7*$wfpu] $wfy $wE $ec
# ------------------------------------------------------
# Define 2D ReinforceConcretePlaneStress element
# ------------------------------------------------------
set j 0;
while {$j < $nH} {
set i 0;
while {$i < [expr $nL]} {
# Create quad elements - command:
# element quad eleID node1 node2 node3 node4 thick type matID } 
element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15
#element quad [expr $j*$nL+($i+$nL/2)+1] [expr $j*($nL+1)+($i+$nL/2)+1] [expr $j*($nL+1)+($i+$nL/2)+2] [expr ($j+1)*($nL+1)+($i+$nL/2)+2] [expr ($j+1)*($nL+1)+($i+$nL/2)+1] $t PlaneStress 16
puts $modelFile "element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15"
set i [expr $i+1]
}
set j [expr $j+1]
}

close $modelFile

model basic -ndm 2 -ndf 3
#puts $modelFile "model basic -ndm 2 -ndf 2"
# ---------------------------
# Define prestress loads
# ---------------------------
# set pForce 1.654e6; # 372kips=1.654e6N
set pForce 1.654e6;
# Creating a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
# Create the nodal load - command: load nodeID xForce yForce
load 36 [expr $pForce/2] 0 0;
load 37 [expr $pForce/2] 0 0;
load 49 [expr -$pForce/2] 0 0;
load 50 [expr -$pForce/2] 0 0;
}
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
# algorithm Newton
# algorithm NewtonLineSearch 0.8
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
integrator LoadControl 0.1
# Creating the analysis object
analysis Static
# initialize in case we need to do an initial stiffness iteration
initialize
# ------------------------------
# End of analysis generation
# ------------------------------
# perform the analysis
analyze 10
# Print out the state of nodes
puts "After Prestress Force"
print node 3 9 15
print node 37 43 49
# Set the prestress loads to be constant & reset the time in the domain
loadConst -time 0.0
# ----------------------------------------------------
# Start of modeling for vertical loads to 200kips
# ----------------------------------------------------
set P1 1000;
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 2 "Linear" {
# Create the nodal load - command: load nodeID xForce yForce
load 53 0 [expr -$P1/3] 0;
load 54 0 [expr -$P1/3] 0;
load 55 0 [expr -$P1/3] 0;
load 65 0 [expr -$P1/3] 0;
load 66 0 [expr -$P1/3] 0;
load 67 0 [expr -$P1/3] 0;
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
# Creating the convergence test
test NormDispIncrVaryIter 0.1 1 5 numStep 2000 numIter 100
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
integrator DisplacementPath 66 2 1 numStep 2000 increment -0.01
# Creating the analysis object analysis Static
# puts "analysis performed"
# initialize in case we need to do an initial stiffness iteration
initialize
# ------------------------------
# End of analysis generation
# ------------------------------
# Creating a recorder to monitor nodal displacements
recorder Node -file B1_m1_disp_Output.out -time -node 3 9 15 -dof 2 disp
# perform the analysis
analyze 2000
# Print out the state of nodes
print node 3 9 15 54
