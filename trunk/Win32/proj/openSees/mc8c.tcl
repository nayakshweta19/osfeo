# Input File for Specimen C8C
# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------
# ------------------------------
# Start of model generation
# ------------------------------
# 4.1
wipe;
# -----------------------------------------------------------------
# Create ModelBuilder for moment carrying walls of hollow bridge columns (with two-dimensions and 3 DOF/node)
# -----------------------------------------------------------------
model basic -ndm 2 -ndf 3
# ------------------------------------------------------------------
# Create nodes needed for moment carrying walls of hollow bridge columns (3dof)
# ------------------------------------------------------------------
# set dimension of the wall and mesh
set L 300;
set H 1240;
set deltaL 100;
set deltaH 155;
set nL [expr $L/$deltaL];
set nH [expr $H/$deltaH];
set t 240;
# Creating nodes for moment carrying walls of hollow bridge columns
# tag X Y
set nodeStartID [expr ($nH+1)*($nL+1)+1];
set eleStartID [expr $nH*$nL*2+1];
set j 0;
while {$j < [expr $nH+1]} {
  if { $j < [expr $nH] } {
    node [expr $nodeStartID+$j*2]    0.0       [expr $deltaH*$j]
    node [expr $nodeStartID+$j*2+1]  [expr $L] [expr $deltaH*$j]
    puts "node [expr $nodeStartID+$j*2]    0.0       [expr $deltaH*$j]"
    puts "node [expr $nodeStartID+$j*2+1]  [expr $L] [expr $deltaH*$j]"
  } else {
    set i 0;
    while {$i < [expr $nL+1]} {
      node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaL] [expr $H]
      puts "node [expr $nodeStartID+$j*2+$i] [expr $i*$deltaL] [expr $H]"
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
  #node [expr $tendonnodeStartID+$j] [expr $deltaL] [expr $deltaH*$j]
  set j [expr $j+1]
}
# Fix supports at base of columns
# tag DX DY RZ
fix [expr $nodeStartID] 1 1 1
fix [expr $nodeStartID+1] 1 1 1
#fix [expr $tendonnodeStartID] 1 1 1

puts "fix [expr $nodeStartID] 1 1 1"
puts "fix [expr $nodeStartID+1] 1 1 1"
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
set fy 300; # Yield stress for bare bar
set E 200000.0; # Young¡¯s modulus
# tag fy E0 fpc rou
uniaxialMaterial Steel01 4 $fy $E 0.01
# tag fy E0 fpc rou epsi
uniaxialMaterial TendonL01 5 1670 $E 1860 0.001 0.0036
# -----------------------------------------------
# Define cross-section for nonlinear columns
# -----------------------------------------------
# set some parameters
set colWidth 360.0
set colDepth 120.0
set cover 25.0
set As 50.8;
# some variables derived from the parameters
set cy1 [expr $colDepth/2.0]
set cz1 [expr $colWidth/2.0]
# the section 1 is for stirrup confinement is #4@80
section Fiber 1 {
  # Create the concrete core fibers
  # mat num num
  patch rect 2 12 4 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover]
  # Create the concrete cover fibers (top, bottom, left, right)
  patch rect 1 1 6 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1
  patch rect 1 1 6 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1]
  patch rect 1 12 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover]
  patch rect 1 12 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover]
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
set colEleIDs "";
set columnType "dispBeamColumn";
set np 2
set iterNum 55
set iterTol 1e-3
set j 0;
while {$j < [expr $nH]} {
  if { $j < [expr $nH-1] } {
  #  element $columnType [expr $eleStartID+$j*2]   [expr $nodeStartID+$j*2]   [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
  #  element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol
    puts "element $columnType [expr $eleStartID+$j*2]   [expr $nodeStartID+$j*2]   [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
    puts "element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol"
    lappend colEleIDs "[expr $eleStartID+$j*2] [expr $eleStartID+$j*2+1]"
  } else {
  #  element $columnType [expr $eleStartID+$j*2]   [expr $nodeStartID+$j*2]   [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
  #  element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+$nL+2] $np 1 1 -iter $iterNum $iterTol
    puts "element $columnType [expr $eleStartID+$j*2]   [expr $nodeStartID+$j*2]   [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
    puts "element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol"
    lappend colEleIDs "[expr $eleStartID+$j*2] [expr $eleStartID+$j*2+1]"
  }
  incr j 1;
}
# -------------------------------------------------------
# Define tendon elements
# -------------------------------------------------------
# geomTransf Linear 1
##set i 0;
##while {$i < [expr $nH]} {
##	#element nonlinearBeamColumn [expr $tendoneleStartID+$i] [expr $tendonnodeStartID+$i] [expr $tendonnodeStartID+1+$i] $np 4 1 -iter $iterNum $iterTol
##	incr i 1
##}
# -------------------------------------------------------
# Define beam elements
# ------------------------------------------------------
geomTransf Linear 2
set j [expr $nH];
set i 0;
while {$i < [expr $nL]} {
  #element $columnType [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 3 2 -iter $iterNum $iterTol
  puts "element nonlinearBeamColumn [expr $eleStartID+$j*2+$i] [expr $nodeStartID+$j*2+$i] [expr $nodeStartID+$j*2+1+$i] $np 3 2 -iter $iterNum $iterTol"
  
  incr i 1
}
# 4.2
# ------------------------------------------------------
# Create ModelBuilder for 2D elements representing shear carrying walls of the column (with two-dimensions and 2 DOF/node)
# ------------------------------------------------------
model basic -ndm 2 -ndf 2
# Create nodes & add to Domain - command: node nodeId xCrd yCrd
set j 0;
while {$j < [expr $nH+1]} {
  set i 0;
  while {$i < [expr $nL+1]} {
    node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]
    puts "node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]"
    if { $j > 0} {
    # tying nodes between moment carrying walls and 2D elements representing shear carrying walls of the column
      equalDOF [expr $j*($nL+1)+$i+1] [expr $nodeStartID+$j*($nL+1)+$i] 1 2;
      puts "equalDOF [expr $j*($nL+1)+$i+1] [expr $nodeStartID+$j*($nL+1)+$i] 1 2";
    }
    incr i 1
  }
  incr j 1
}
# Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
set i 0;
while {$i < [expr $nL+1]} {
  fix [expr $i+1] 1 1
  puts "fix [expr $i+1] 1 1"
  incr i 1
}

# -----------------------------------------------------
# Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------
# set fc fy E
set wfc 30.0;
set wfyv 434;
set wfyh1 433;
set wE 200000.0;
set rou1 0;
set rou2 0;
set rouv 0.025;
set rouh1 0.0246; # #10@70
# UniaxialMaterial: steelZ01
# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteL02 14 [expr -$wfc] -0.003
uniaxialMaterial ConcreteL02 15 [expr -$wfc] -0.003
set pi 3.141592654
# NDMaterial: FAFourSteelPCPlaneStress
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelPCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteetPCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteetPCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAFourSteelRCPlaneStress  (100)
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.0*$pi]  [expr 0.5*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
nDMaterial CSMMRCPlaneStress 21 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.0*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.003

# --------------------------------------------------
# Define 2D ReinforceConcretePlaneStress element
# --------------------------------------------------
set j 0;
set quadEleIDs "";
while {$j < $nH} {
	set i 0;
	while {$i < $nL} {
		# Create quad elements - command:
		# element quad eleID node1 node2 node3 node4 thick type matID
		element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 21
		#element quad [expr $nH*$nL+$j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 21
		puts "element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1]"
		#puts "element quad [expr $nH*$nL+$j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 21"
		
		lappend quadEleIDs "[expr $j*$nL+$i+1]"
		incr i 1
	}
	incr j 1
}
# 4.3
# -------------------------------------
# Define prestress and gravity loads
# -------------------------------------
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
	# Create the nodal load - command: load nodeID xForce yForce
  set i 0;
  while {$i < [expr $nL+1]} {
	load [expr $nH*($nL+1)+$i+1] 0 [expr -518400./($nL+1)]
	puts "load [expr $nH*($nL+1)+$i+1] 0 [expr -518400./($nL+1)]"
    incr i 1;
  }
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
#test NormDispIncr 1.0e-3 1000 5
test EnergyIncr  1.0e-12    10         0
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
integrator LoadControl 0.1
# Creating the analysis object
analysis Static
# ------------------------------
# End of analysis generation
# ------------------------------
# create the display
#remove recorders;
recorder display nodeNum 10 10 600 420 -wipe
# next three commmands define viewing system, all values in global coords
prp -8000 -1000 4000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 1 0;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -1400 1400 -1400 1400; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 -1 10; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# create the display
#remove recorders
recorder display Deform 620 10 600 420 -wipe
# next three commmands define viewing system, all values in global coords
prp -1000 -10000 7000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 1 0;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -1200 1200 -1200 1200; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 -1 10; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# perform the analysis
analyze 10
puts "Graivty analysis passed!"
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
  set i 0;
  while {$i < [expr $nL+1]} {
    load [expr $nH*($nL+1)+$i+1] 1.0 0.0
    puts "load [expr $nH*($nL+1)+$i+1] 1.0 0.0"
    incr i 1
  }
}
set nodeIDs [getNodeTags];
set eleIDs [getEleTags];
# Creating a recorder to monitor nodal displacements

set cmd "recorder Element -file mC8C_Col_S_E.dat -time -ele $colEleIDs section fiber -$cy1 0.0 1 stressStrain";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p1_S.dat -time -ele $quadEleIDs integrPoint 1 stresses";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p1_E.dat -time -ele $quadEleIDs integrPoint 1 strains ";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p2_S.dat -time -ele $quadEleIDs integrPoint 2 stresses";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p2_E.dat -time -ele $quadEleIDs integrPoint 2 strains ";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p3_S.dat -time -ele $quadEleIDs integrPoint 3 stresses";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p3_E.dat -time -ele $quadEleIDs integrPoint 3 strains ";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p4_S.dat -time -ele $quadEleIDs integrPoint 4 stresses";
eval $cmd;
set cmd "recorder Element -file mC8C_ele_p4_E.dat -time -ele $quadEleIDs integrPoint 4 strains ";
eval $cmd;

set cmd "recorder Node -file mC8C_disp_Output.dat -time -node $nodeIDs -dof 1 2 disp";
eval $cmd;

# ------------------------------
# End of model generation
# ------------------------------
# 4.6
# ------------------------------
# Start of analysis generation
# ------------------------------
# ---------------------------------- 
# perform the analysis
# ---------------------------------- 
set IDctrlNode [expr ($nH+1)*($nL+1)];
set IDctrlDOF 1;
set Dincr 0.01;
# Constraint Handler Lagrange Transformation
constraints  Plain
# Convergence Test 
test  NormDispIncr  +1.000000E-003    250     0
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
puts "integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr"
# DOF Numberer 
numberer  Plain 
# System of Equations SparseGeneralBandGeneralUmfPack
system BandGeneral
# Analysis Type 
analysis  Static 

set rcdcmd "recorder  Node  -file  PushoverCase_KeyNode_Dsp.out  -time -node $IDctrlNode -dof $IDctrlDOF  disp "
eval $rcdcmd;

#  perform Static Cyclic Displacements Analysis   /////////////////////////////////////////////////////////////
#                     5     -5     10    -10    +15    -15    20   -20    25     -25   30    -30     35    -35   40   -40     45     -45   50   -50      55    -55    60    -60     65    -65   70   -70     0
set numSteps     {    50    100   150     200   250    300   300    400   450    500   550    600    650   700    750   800    850   900    950  1000   1050   1100   1150  1200   1250  1300   1350  1400   700      }
set numIters     {   100    500   400     400   600    400   200    200   200    200   500    500    200   200    200   500    500   200    200   500    500    200    200   500    500   200    200   500   500      }
set increments   {   0.1   -0.1   0.1    -0.1   0.1  -0.01  0.01  -0.01  0.01  -0.01  0.01  -0.01   0.01 -0.01   0.01 -0.01   0.01 -0.01   0.01 -0.01   0.01  -0.01   0.01 -0.01   0.01 -0.01   0.01 -0.01  0.01      }
 
for { set i 0 } { $i<29 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set Dincr [lindex $increments $i]
    integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
    if {$numIter == 0} {
	set numIter 100
    }
    test NormDispIncr 1.0e-3 $numIter 5
    #analyze $numStep

    set j 0;
    while { $j < $numStep } {
      set ok [analyze 1]
      # ----------------------------------------------if convergence failure-------------------------
      # if analysis fails, we try some other stuff
      # performance is slower inside this loop  global maxNumIterStatic;      # max no. of iterations performed before "failure to converge" is ret'd
      if {$ok != 0} {
        puts "Trying Newton with Initial Tangent .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm Newton -initial
        set ok [analyze 1]
        test NormDispIncr 1e-3 [expr $numIter*2] 5
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying Broyden .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm Broyden 8
        set ok [analyze 1 ]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying NewtonWithLineSearch .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm NewtonLineSearch 0.8 
        set ok [analyze 1]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        set putout "PROBLEM: $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF]";
        puts $putout
        test NormDispIncr 1e-3 $numIter 5
        set ok [analyze 1]
      }; # end if
      incr j 1;
      #record;
      #puts $end
    }

    puts $i
}
#remove recorders;
