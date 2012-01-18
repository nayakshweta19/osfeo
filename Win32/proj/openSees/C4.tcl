wipe;
source units&constants.tcl
set modelFile [open "model.ops" "w"]
# Input File for Specimen C8C
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
set L 375;
set H 1600;
set deltaL 375;
set deltaH 400;
set nL [expr $L/$deltaL];
set nH [expr $H/$deltaH];
set t 150;
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

# Fix supports at base of columns
# tag DX DY RZ
fix [expr $nodeStartID] 1 1 1
fix [expr $nodeStartID+1] 1 1 1
puts $modelFile "fix [expr $nodeStartID] 1 1 1"
puts $modelFile "fix [expr $nodeStartID+1] 1 1 1"

# -----------------------------------------------------------
# Define nonlinear materials for moment carrying walls and tendons
# -----------------------------------------------------------
# CONCRETE tag f¡¯c ec0 f¡¯cu ecu
# Cover concrete (unconfined)
uniaxialMaterial Concrete01 1 -24.800000000000004 -0.002345151796691333 -24.640000000000004 -0.04 
# Core concrete (confined)
uniaxialMaterial Concrete02 2 -28.4 -0.00254944878790797 -7.28 -0.04 0.1 3.64 1820.0

# STEEL
# Longitudinal Reinforcing steel
set fy 560; # Yield stress for bare bar
set E 190000.0; # Young¡¯s modulus
# tag fy E0 fpc rou
uniaxialMaterial ReinforcingSteel 3 430.0 494.49999999999994 199948.00558280002 19994.80055828 0.008 0.03  -GABuck 6.0 2.0 0.1 0.25  -CMFatigue 0.2 0.5 0.4  -IsoHard 4.3 0.01  -MPCurveParams 0.333 18.0 4.0

uniaxialMaterial Steel02 4 $fy $E 0.05 18.0 0.925 0.15 0.0 1.0 0.0 1.0 0.0

# -----------------------------------------------
# Define cross-section for nonlinear columns
# -----------------------------------------------
# set some parameters
set colWidth 825.0
set colDepth 75.0
set cover 10.0
set As [expr 53*13];
# some variables derived from the parameters
set cy1 [expr $colDepth/2.0]
set cz1 [expr $colWidth/2.0]
# the section 1 is for stirrup confinement is #4@80
section Fiber 1 -GJ 1189799672666.6667 { 
layer straight 4 13 50.3 29.5 442.0 29.5 -442.0
patch quad 1 16 4 -37.5 450.0 -37.5 -450.0 -29.5 -442.0 -29.5 442.0
layer straight 4 0 50.3 29.5 442.0 29.5 -442.0
patch quad 1 4 16 -29.5 -442.0 -37.5 -450.0 37.5 -450.0 29.5 -442.0
patch quad 2 4 4 -29.5 442.0 -29.5 -442.0 29.5 -442.0 29.5 442.0
layer straight 4 13 50.3 -29.5 442.0 -29.5 -442.0
patch quad 1 16 4 37.5 450.0 29.5 442.0 29.5 -442.0 37.5 -450.0
patch quad 1 4 16 -37.5 450.0 -29.5 442.0 29.5 442.0 37.5 450.0
};
# -----------------------------------------------
# Define cross-section for top beam
# -----------------------------------------------
# E A I
section Elastic 3 2e5 1e6 1e14

# -------------------------------------------------------
# Define column elements
# -------------------------------------------------------
geomTransf Corotational 1
set columnType nonlinearBeamColumn ;  # dispBeamColumn
set np 2
set iterNum 150
set iterTol 1e-3
set j 0;
while {$j < [expr $nH]} {
	if {$j < [expr $nH-1]} {
		element $columnType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
		element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol
		puts $modelFile "element $columnType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
		puts $modelFile "element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$j*2+3] $np 1 1 -iter $iterNum $iterTol"
	} else {
		element $columnType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol
		element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol
		puts $modelFile "element $columnType [expr $eleStartID+$j*2] [expr $nodeStartID+$j*2] [expr $nodeStartID+$j*2+2] $np 1 1 -iter $iterNum $iterTol"
		puts $modelFile "element $columnType [expr $eleStartID+$j*2+1] [expr $nodeStartID+$j*2+1] [expr $nodeStartID+$nH*2+$nL] $np 1 1 -iter $iterNum $iterTol"
	}
	set j [expr $j+1]
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
equalDOF 3 13 1 2
equalDOF 4 14 1 2
equalDOF 5 15 1 2
equalDOF 6 16 1 2
equalDOF 7 17 1 2
equalDOF 8 18 1 2
equalDOF 9 19 1 2
equalDOF 10 20 1 2

# -----------------------------------------------------
# Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------
# set fc fy E
set wfc 28.0;
set wfyv 560;
set wfyh1 450;
set wE 190000.0;
set rou1 0;
set rou2 0;
set rouv 0.0246;
set rouh1 0.00246; # #10@70
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
#nDMaterial FAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 $rouv $rouh1 $wfc 0.0 $wfyv $wE 0.003

#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
nDMaterial CSMMRCPlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi] [expr 0.0*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002

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
		#element SSPquad eleTag? iNode? jNode? kNode? lNode? matTag? type? thickness? <b1? b2?>?
		#element SSPquad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] 21 PlaneStress $t 
		puts $modelFile "element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1]"
		set i [expr $i+1]
	}
	set j [expr $j+1]
}
close $modelFile;
# 4.3
# -------------------------------------
# Define prestress and gravity loads
# -------------------------------------
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
	# Create the nodal load - command: load nodeID xForce yForce
	load 9  0 [expr -125e3]
	load 10 0 [expr -125e3]

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
test NormDispIncr 1.0e-3 1000 5
#test  EnergyIncr  +1.000000E-004    25     0     2 
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
integrator LoadControl 0.1
# Creating the analysis object
analysis Static
# ------------------------------
# End of analysis generation
# ------------------------------
recorder Node -file C4_Grav_disp_Output.dat -time -nodeRange 1 20 -dof 1 2 disp
# perform the analysis
analyze 10
# Print out the state of nodes
print node 9 10
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
set P 1.0;
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 2 "Linear" {
	# Create the nodal load - command: load nodeID xForce yForce
	load 9 [expr $P/2] 0
	load 10 [expr $P/2] 0
}
# Creating a recorder to monitor nodal displacements
recorder Node -file C4_disp_Output.dat -time -nodeRange 1 20 -dof 1 2 disp

recorder Element -file C4_Col_S_E.dat -time -eleRange 5 12 section fiber -$cy1 0.0 1 stressStrain

recorder Element -file C4_ele_p1_S.dat -time -ele 1 2 3 4 integrPoint 1 stresses
recorder Element -file C4_ele_p1_E.dat -time -ele 1 2 3 4 integrPoint 1 strains
recorder Element -file C4_ele_p2_S.dat -time -ele 1 2 3 4 integrPoint 2 stresses
recorder Element -file C4_ele_p2_E.dat -time -ele 1 2 3 4 integrPoint 2 strains
recorder Element -file C4_ele_p3_S.dat -time -ele 1 2 3 4 integrPoint 3 stresses
recorder Element -file C4_ele_p3_E.dat -time -ele 1 2 3 4 integrPoint 3 strains
recorder Element -file C4_ele_p4_S.dat -time -ele 1 2 3 4 integrPoint 4 stresses
recorder Element -file C4_ele_p4_E.dat -time -ele 1 2 3 4 integrPoint 4 strains

# create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -1750 -5800 6000
  #vrp 0 -500 250
  vup 0 1 0
  #vpn -1 -1 0.5
  viewWindow -4000 4000 -4000 4000
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 10

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
test NormDispIncr 1.e-3 100 0
# Creating the solution algorithm
algorithm KrylovNewton
# Creating the integration scheme, the DisplacementControl scheme
#integrator DisplacementControl $node $dof $incr <$numIter $¦¤Umin $¦¤Umax $maxLambda>
integrator DisplacementControl 9 1 0.1
# Creating the analysis object
analysis Static
# ----------- set up analysis parameters

#source LibAnalysisStaticParameters.tcl;
# initialize in case we need to do an initial stiffness iteration
initialize

# ----------- set up analysis parameters
# characteristics of pushover analysis
set Dincr [expr 0.1];	# displacement increment. you want this to be small, but not too small to slow analysis
set IDctrlNode 9;
set IDctrlDOF 1;
set analysisTypeStatic Static
variable TolStatic 1.0e-3;                        # Convergence Test: tolerance
variable Tol [expr $TolStatic];
variable maxNumIterStatic 2000;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
variable printFlagStatic 0;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
variable testTypeStatic NormDispIncr ;	# Convergence-test type: NormDispIncr, EnergyIncr
variable maxNumIterConvergeStatic 2000;	
variable printFlagConvergeStatic 0;
variable algorithmTypeStatic KrylovNewton
# perform the analysis  0  1  2   3     4    5    6   7    8   9   10    11   12   13    14   15   16   17
#                1.0   -1.0 +2.0  -2.0 +3.0 -3.0 4.0 -4.0 5.0 -5.0 -6.0  6.0 -7.0  8.0 -10.0 15.0 -25.0 40.0  
##set numSteps   {  10     40   30    40  50    60  70   80  90  100  110  120  130  150   180  250   400  650 };
##set numIters   { 200    300  400   400  450  450 500  550 550  550  550  550  550  550   550  550   550  550 };
##set increments { 0.1  -0.05  0.1  -0.1  0.1 -0.1 0.1 -0.1 0.1 -0.1 -0.1  0.1 -0.1  0.1  -0.1  0.1  -0.1  0.1 };
set iDmax { 1 3 5 10 4 14 17 7 25 30 33 40 45 } ; # vector of displacement-cycle peaks, in terms of storey drift ratio
set CycleType Full;    # you can do Full / Push / Half cycles with the proc
set Ncycles 2;				# specify the number of cycles at each peak
#  ---------------------------------    perform Static Cyclic Displacements Analysis
source GeneratePeaks.tcl

# -- STATIC PUSHOVER/CYCLIC ANALYSIS

foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep            
			puts "D1=$D1"
			set Dincr [expr $D1 - $D0]
			#puts "Dincr = $Dincr"
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			puts "Displacement Node: $IDctrlNode. Dof: $IDctrlDOF. Incr: $Dincr"
			#integrator ArcLength 1.0 0.1
			# ------------------first analyze command------------------------
			set ok [analyze 1]
			# -----------------if convergence failure-------------------------
			# max no. of iterations performed before "failure to converge" is ret'd
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent .."
				test NormDispIncr $Tol 2000 0;
				algorithm Newton -initial;
				set ok [analyze 1];
				test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;
				algorithm $algorithmTypeStatic;
				puts "Trying Newton with Initial Tangent failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Broyden ..";
				algorithm Broyden 40
				set ok [analyze 1 ]
				algorithm $algorithmTypeStatic
				puts "Trying Broyden failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Newton With LineSearch ..";
				algorithm NewtonLineSearch 0.8 
				set ok [analyze 1]
				algorithm $algorithmTypeStatic
				puts "Trying Newton With LineSearch failed to converge...";
			}
			#if {$ok != 0} {
			#	puts "Trying ArcLength ..";
			#	algorithm $algorithmTypeStatic 
			#	integrator ArcLength 1.0 0.1
			#	set ok [analyze 1]
			#	puts "Trying ArcLength With KrylovNewton failed to converge...";
			#}
			if {$ok != 0} {
				set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";
				set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] m]
				puts $putout
				return -1
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
			puts "$D0";
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
remove recorders;