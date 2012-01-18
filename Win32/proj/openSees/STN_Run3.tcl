# STN Run3.tcl
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
# ------------------------------
# Start of model generation
# ------------------------------
# --------------------------------------------------------
# Create ModelBuilder for beams and columns
# --------------------------------------------------------
model basic -ndm 2 -ndf 3
# ---------------------------------------------------------
# Create nodes for beams and columns
# ---------------------------------------------------------
# set dimension of the wall and mesh
set L 1158;
set H 699;
set deltaL 193;
set deltaH 233;
set nL [expr $L/$deltaL];
set nH [expr $H/$deltaH];
set t 60;
# Create nodes for beams and columns
# tag X Y
set LB 1278;
set deltaLB 213;
set nodeStartID [expr ($nH+1)*($nL+1)+1];
set eleStartID [expr $nH*$nL+1]
set j 0;
while {$j < [expr $nH+1]} {
if {$j < [expr $nH] } {
node [expr $nodeStartID+$j*2] 0.0 [expr $deltaH*$j]
node [expr $nodeStartID+$j*2+1] [expr $LB] [expr $deltaH*$j]
} else {
set i 0;
while {$i < [expr $nL+1]} {
node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaLB] [expr $H]
set i [expr $i+1]
}
}
set j [expr $j+1]
}
# --------------------------------------------------
# Define materials for nonlinear columns and beams
# --------------------------------------------------
set ec 0.003;
# CONCRETE tag f¡¯c ec0 f¡¯cu ecu
# Core concrete (confined)
uniaxialMaterial Concrete01 1 -34.0 [expr -$ec] -6.8 -0.008
# Cover concrete (unconfined)
uniaxialMaterial Concrete01 2 -34.0 [expr -$ec] -0.0 -0.006
# STEEL
# Reinforcing steel
set fy 360.0; # Smeared yield stress for steel bars
set E 200000.0; # Young¡¯s modulus
# tag fy E0 b
uniaxialMaterial Steel01 3 $fy $E 0.02
# -----------------------------------------------
# Define cross-section for columns
# -----------------------------------------------
# set some parameters
set colWidth 117.0
set colDepth 117.0
set cover 20.0
set As1 197.9;
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
layer straight 3 2 $As1 [expr $cy1-$cover] [expr $cz1-$cover] [expr $cy1-$cover] [expr $cover-$cz1]
layer straight 3 2 $As1 [expr $cover-$cy1] [expr $cz1-$cover] [expr $cover-$cy1] [expr $cover-$cz1]
}
# -----------------------------------------------
# Define cross-section for beams
# -----------------------------------------------
# set some parameters
set beamWidth 480
set beamDepth 150
set As2 387.9;
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
layer straight 3 3 $As2 [expr $by1-$cover] [expr $bz1-$cover] [expr $by1-$cover] [expr $cover-$bz1]
layer straight 3 3 $As2 [expr $cover-$by1] [expr $bz1-$cover] [expr $cover-$by1] [expr $cover-$bz1]
}
# -------------------------------------------------------
# Define column elements
# -------------------------------------------------------
geomTransf Linear 1
set np 2
set iterNum 10
set iterTol 1e-3
set j 0;
while {$j < [expr $nH]} {
if {$j < [expr $nH-1]} {
# define the columns elements
element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol
} else {
element nonlinearBeamColumn [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
element nonlinearBeamColumn [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol
}
set j [expr $j+1]
}
# -------------------------------------------------------
# Define beam elements
# ------------------------------------------------------
geomTransf Linear 2
set j [expr $nH];
set bA 48000;
set bE 2e15;
set bI 1e10;
set i 0;
while {$i < [expr $nL]} {
element nonlinearBeamColumn [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 2 2 -iter $iterNum $iterTol
set i [expr $i+1]
}
# ------------------------------------------
# Define cantilever beams on both sides
# ------------------------------------------
node 42 -500.0 699.0
node 43 -250.0 699.0
node 44 1530.0 699.0
node 45 1780.0 699.0
element nonlinearBeamColumn 31 42 43 $np 2 2 -iter $iterNum $iterTol
element nonlinearBeamColumn 32 43 35 $np 2 2 -iter $iterNum $iterTol
element nonlinearBeamColumn 33 41 44 $np 2 2 -iter $iterNum $iterTol
element nonlinearBeamColumn 34 44 45 $np 2 2 -iter $iterNum $iterTol
# ------------------------------------------------------
# Create ModelBuilder for wall elements
# ------------------------------------------------------
model basic -ndm 2 -ndf 2
# Create nodes & add to Domain - command: node nodeId xCrd yCrd
set j 0;
while {$j < [expr $nH+1]} {
set i 0;
while {$i < [expr $nL+1]} {
node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]
set i [expr $i+1]
}
set j [expr $j+1]
}
# Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
set i 0;
while {$i < [expr $nL+1]} {
fix [expr $i+1] 1 1
set i [expr $i+1]
}
# tie nodes between beams, columns and wall elements
equalDOF 1 29 1 2
equalDOF 7 30 1 2
equalDOF 8 31 1 2
equalDOF 14 32 1 2
equalDOF 15 33 1 2
equalDOF 21 34 1 2
equalDOF 22 35 1 2
equalDOF 23 36 1 2
equalDOF 24 37 1 2
equalDOF 25 38 1 2
equalDOF 26 39 1 2
equalDOF 27 40 1 2
equalDOF 28 41 1 2
# -----------------------------------------------------
# Define materials for wall elements
# -----------------------------------------------------
# set fc fy E
set wfc 34.0;
set wfy 371.2;
set wE 200000.0;
set rou1 0.005;
set rou2 0.005;
# UniaxialMaterial: steelZ01
# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfy $wE $wfc $rou1
uniaxialMaterial SteelZ01 12 $wfy $wE $wfc $rou2
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteZ01 13 [expr -$wfc] [expr -$ec]
uniaxialMaterial ConcreteZ01 14 [expr -$wfc] [expr -$ec]
set pi 3.141592654
# nDMaterial FAReinforceConcretePlaneStress
# tag rho s1 s2 c1 c2 angle1 angle2 rou1 rou2 fpc fy E0 ec
nDMaterial FAReinforceConcretePlaneStress 15 0.0 11 12 13 14 [expr 0.25*$pi] [expr 0.75*$pi] $rou1 $rou2 $wfc $wfy $wE $ec
# --------------------------------------------------------
# Define wall elements
# --------------------------------------------------------
set j 0;
while {$j < $nH} {
set i 0;
while {$i < $nL} {
# Create quad elements - command:
# element quad eleID node1 node2 node3 node4 thick type matID
element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15
set i [expr $i+1]
}
set j [expr $j+1]
}
# -----------------------------------------------------
# Define dynamic analysis parameters
# -----------------------------------------------------
model basic -ndm 2 -ndf 3
# define masses
set m 1.86;
mass 35 $m $m 0
mass 36 $m $m 0
mass 37 $m $m 0
mass 38 $m $m 0
mass 39 $m $m 0
mass 40 $m $m 0
mass 41 $m $m 0
mass 42 $m $m 0
mass 43 $m $m 0
mass 44 $m $m 0
mass 45 $m $m 0
# define earthquake excitation
set inputAccel "Path -filePath STN_run3_input_accel.TXT -dt 0.005 -factor [expr -30000.0]"
pattern UniformExcitation 2 1 -accel $inputAccel
# define damping
set dampingratio 0.04;
set a1 [expr $dampingratio*0.0073];
rayleigh 0.0 0.0 0.0 $a1
# ------------------------------
# End of model generation
# ------------------------------
# ------------------------------
# Start of analysis generation
# ------------------------------
# Create the system of equation
system BandGeneral
# Create the constraint handler
constraints Plain
# Create the DOF numberer
numberer Plain
# Create the convergence test
test NormDispIncr 0.01 100 5
# Create the solution algorithm
algorithm KrylovNewton
# Create the integration scheme
integrator WilsonTheta 1.42
# Create the analysis object
analysis Transient
# ------------------------------
# End of analysis generation
# ------------------------------
# ------------------------------------
# Create recorder and perform analysis
# ------------------------------------
# Create a recorder to record lateral displacement on the top of the central wall
recorder Node -file STN_run3_disp.out -time -node 38 -dof 1 disp
# Perform analysis
analyze 15000 0.005