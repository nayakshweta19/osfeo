######################## 
# Analysis-Sequence  1 #
######################## 
wipe;
# Start of model generation 
# ========================= 
source units&constants.tcl
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  3  -ndf  6

# NodeCoord.tcl 

# Node    tag    xCrd    yCrd 
node       1  +0.000000E+000  +0.000000E+000 0.0
node       2  +1.000000E+003  +0.000000E+000 0.0
node       3  +2.000000E+003  +0.000000E+000 0.0
node       4  +3.000000E+003  +0.000000E+000 0.0
node       5  +4.000000E+003  +0.000000E+000 0.0
node       6  +5.000000E+003  +0.000000E+000 0.0
node       7  +6.000000E+003  +0.000000E+000 0.0
node       8  +7.000000E+003  +0.000000E+000 0.0
node       9  +8.000000E+003  +0.000000E+000 0.0
node      10  +0.000000E+000  +1.000000E+003 0.0
node      11  +1.000000E+003  +1.000000E+003 0.0
node      12  +2.000000E+003  +1.000000E+003 0.0
node      13  +3.000000E+003  +1.000000E+003 0.0
node      14  +4.000000E+003  +1.000000E+003 0.0
node      15  +5.000000E+003  +1.000000E+003 0.0
node      16  +6.000000E+003  +1.000000E+003 0.0
node      17  +7.000000E+003  +1.000000E+003 0.0
node      18  +8.000000E+003  +1.000000E+003 0.0
node      19  +0.000000E+000  +2.000000E+003 0.0
node      20  +1.000000E+003  +2.000000E+003 0.0
node      21  +2.000000E+003  +2.000000E+003 0.0
node      22  +3.000000E+003  +2.000000E+003 0.0
node      23  +4.000000E+003  +2.000000E+003 0.0
node      24  +5.000000E+003  +2.000000E+003 0.0
node      25  +6.000000E+003  +2.000000E+003 0.0
node      26  +7.000000E+003  +2.000000E+003 0.0
node      27  +8.000000E+003  +2.000000E+003 0.0
node      28  +0.000000E+000  +3.000000E+003 0.0
node      29  +1.000000E+003  +3.000000E+003 0.0
node      30  +2.000000E+003  +3.000000E+003 0.0
node      31  +3.000000E+003  +3.000000E+003 0.0
node      32  +4.000000E+003  +3.000000E+003 0.0
node      33  +5.000000E+003  +3.000000E+003 0.0
node      34  +6.000000E+003  +3.000000E+003 0.0
node      35  +7.000000E+003  +3.000000E+003 0.0
node      36  +8.000000E+003  +3.000000E+003 0.0
node      37  +0.000000E+000  +4.000000E+003 0.0
node      38  +1.000000E+003  +4.000000E+003 0.0
node      39  +2.000000E+003  +4.000000E+003 0.0
node      40  +3.000000E+003  +4.000000E+003 0.0
node      41  +4.000000E+003  +4.000000E+003 0.0
node      42  +5.000000E+003  +4.000000E+003 0.0
node      43  +6.000000E+003  +4.000000E+003 0.0
node      44  +7.000000E+003  +4.000000E+003 0.0
node      45  +8.000000E+003  +4.000000E+003 0.0
node      46  +0.000000E+000  +5.000000E+003 0.0
node      47  +1.000000E+003  +5.000000E+003 0.0
node      48  +2.000000E+003  +5.000000E+003 0.0
node      49  +3.000000E+003  +5.000000E+003 0.0
node      50  +4.000000E+003  +5.000000E+003 0.0
node      51  +5.000000E+003  +5.000000E+003 0.0
node      52  +6.000000E+003  +5.000000E+003 0.0
node      53  +7.000000E+003  +5.000000E+003 0.0
node      54  +8.000000E+003  +5.000000E+003 0.0
node      55  +0.000000E+000  +6.000000E+003 0.0
node      56  +1.000000E+003  +6.000000E+003 0.0
node      57  +2.000000E+003  +6.000000E+003 0.0
node      58  +3.000000E+003  +6.000000E+003 0.0
node      59  +4.000000E+003  +6.000000E+003 0.0
node      60  +5.000000E+003  +6.000000E+003 0.0
node      61  +6.000000E+003  +6.000000E+003 0.0
node      62  +7.000000E+003  +6.000000E+003 0.0
node      63  +8.000000E+003  +6.000000E+003 0.0


# SPConstraint.tcl 

# SPC    tag    Dx    Dy 
fix       1     1     1     1    1 1 1 
fix       2     1     1     1    1 1 1 
fix       3     1     1     1    1 1 1 
fix       4     1     1     1    1 1 1 
fix       5     1     1     1    1 1 1 
fix       6     1     1     1    1 1 1 
fix       7     1     1     1    1 1 1 
fix       8     1     1     1    1 1 1 
fix       9     1     1     1    1 1 1 
fix      10     1     1     1    1 1 1 
fix      18     1     1     1    1 1 1 
fix      19     1     1     1    1 1 1 
fix      27     1     1     1    1 1 1 
fix      28     1     1     1    1 1 1 
fix      36     1     1     1    1 1 1 
fix      37     1     1     1    1 1 1 
fix      45     1     1     1    1 1 1 
fix      46     1     1     1    1 1 1 
fix      54     1     1     1    1 1 1 
fix      55     1     1     1    1 1 1 
fix      56     1     1     1    1 1 1 
fix      57     1     1     1    1 1 1 
fix      58     1     1     1    1 1 1 
fix      59     1     1     1    1 1 1 
fix      60     1     1     1    1 1 1 
fix      61     1     1     1    1 1 1 
fix      62     1     1     1    1 1 1 
fix      63     1     1     1    1 1 1 

# Define nodal masses 
# ------------------- 

# Define Multi Point Constraints 
# ------------------------------ 

# Define material(s) 
# ------------------ 
# create the material
nDMaterial ElasticIsotropic   1   2e5   0.25  6.75 

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 34.5;  #[expr (6880.+6600.+7700.+6850.)/4.0*$psi/$MPa];
set wfy 414.0; #[expr (76.+63.+60)/3.0*$ksi/$MPa];
set wE  2.0e5; #[expr (29000+28000+28700)/3.0*$ksi/$MPa];
set rou1 0.0359;
set rou2 0.0226;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
#uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
#uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2

# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial Concrete09 13 -2.000000E+001  -4.000000E-003  -2.800000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +2.400000E+004 
#uniaxialMaterial Concrete09 14 -2.000000E+001  -4.000000E-003  -2.800000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +2.400000E+004 

set pi 3.141592654;
# NDMaterial: ReinforceConcretePlateFiber
#                                        tag    rho   s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
nDMaterial FAReinforcedConcretePlateFiber 2  2.7e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0035
# Material "Material01":    matTag    K    G    sig0    sigInf    delta    H    eta 
nDMaterial  J2Plasticity       4  +1.596000E+004  +1.336200E+004  +2.000000E+001  30. 0.05    +2.000000E-002  +0.000000E+000 
# ---------------------------------------
# Define Section tag
# ---------------------------------------
#section PlateFiber secTag  ndMatTag  h
set eleSecTag 2;
set thick [expr 150.0*$mm];
section  PlateFiber  $eleSecTag   2  $thick


# Define element(s) 
# ----------------- 
# Elements.tcl 

# Element "shell":    eleTag    NodeI    NodeJ    NodeK    NodeL    secTag 
element  shell02       1       1       2      11      10       2 
element  shell02       2       2       3      12      11       2 
element  shell02       3       3       4      13      12       2 
element  shell02       4       4       5      14      13       2 
element  shell02       5       5       6      15      14       2 
element  shell02       6       6       7      16      15       2 
element  shell02       7       7       8      17      16       2 
element  shell02       8       8       9      18      17       2 
element  shell02       9      10      11      20      19       2 
element  shell02      10      11      12      21      20       2 
element  shell02      11      12      13      22      21       2 
element  shell02      12      13      14      23      22       2 
element  shell02      13      14      15      24      23       2 
element  shell02      14      15      16      25      24       2 
element  shell02      15      16      17      26      25       2 
element  shell02      16      17      18      27      26       2 
element  shell02      17      19      20      29      28       2 
element  shell02      18      20      21      30      29       2 
element  shell02      19      21      22      31      30       2 
element  shell02      20      22      23      32      31       2 
element  shell02      21      23      24      33      32       2 
element  shell02      22      24      25      34      33       2 
element  shell02      23      25      26      35      34       2 
element  shell02      24      26      27      36      35       2 
element  shell02      25      28      29      38      37       2 
element  shell02      26      29      30      39      38       2 
element  shell02      27      30      31      40      39       2 
element  shell02      28      31      32      41      40       2 
element  shell02      29      32      33      42      41       2 
element  shell02      30      33      34      43      42       2 
element  shell02      31      34      35      44      43       2 
element  shell02      32      35      36      45      44       2 
element  shell02      33      37      38      47      46       2 
element  shell02      34      38      39      48      47       2 
element  shell02      35      39      40      49      48       2 
element  shell02      36      40      41      50      49       2 
element  shell02      37      41      42      51      50       2 
element  shell02      38      42      43      52      51       2 
element  shell02      39      43      44      53      52       2 
element  shell02      40      44      45      54      53       2 
element  shell02      41      46      47      56      55       2 
element  shell02      42      47      48      57      56       2 
element  shell02      43      48      49      58      57       2 
element  shell02      44      49      50      59      58       2 
element  shell02      45      50      51      60      59       2 
element  shell02      46      51      52      61      60       2 
element  shell02      47      52      53      62      61       2 
element  shell02      48      53      54      63      62       2 

# Define time series 
# ------------------ 
# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# Start of anaysis generation 
# =========================== 

# Get Initial Stiffness 
# --------------------- 
initialize 

# Analysis: StaticDefaultCase 
# +++++++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    for {set i 1} { $i <= 63} {incr i } {
    	load $i   +0.000000E+000  +0.000000E+000 -1.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
    }
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    xL    <Px> 
} 

# Define recorder(s) 
# -------------------- 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  63 -dof  1  2  3  disp 

# create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -17500 -58000 60000
  #vrp 0 -500 250
  vup 0 0 1
  #vpn -1 -1 0.5
  viewWindow -4000 4000 -4000 4000
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 10

# Define analysis options 
# ----------------------- 
# Constraint Handler 
constraints   Plain 
# Convergence Test 
#test  NormUnbalance  +1.000000E-006    25     0     2 
# Integrator 
#integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
#algorithm  Newton 


# Convergence Test 
test  NormDispIncr  +1.000000E-002    25     0     2 
# Integrator 
integrator  DisplacementControl 32 3 +30.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 

# DOF Numberer 
numberer  RCM 
# System of Equations 
system  ProfileSPD 
# Analysis Type 
analysis  Static 

analyze 1

# perform the analysis  0  1  2   3     4    5    6   7    8   9   10    11   12   13    14   15   16   17
#                1.0   -1.0 +2.0  -2.0 +3.0 -3.0 4.0 -4.0 5.0 -5.0 -6.0  6.0 -7.0  8.0 -10.0 15.0 -25.0 40.0  
set numSteps   {  10     40   30    40  50    60  70   80  90  100  110  120  130  150   180  250   400  650 }; #   800 1200  1600  2050  2500 2960  3470 1000    10}  
set numIters   { 200    300  400   400  450  450 500  550 550  550  550  550  550  550   550  550   550  550 }; #   100  100   100   100   100  500   500  500   300}  
set increments { 0.1  -0.05  0.1  -0.1  0.1 -0.1 0.1 -0.1 0.1 -0.1 -0.1  0.1 -0.1  0.1  -0.1  0.1  -0.1  0.1 }; # -0.01 0.01 -0.01  0.01 -0.01 0.01 -0.01 0.01 -0.01}  

for { set i 0 } { $i<18 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set increment [lindex $increments $i]

    integrator DisplacementControl 32 3 $increment
    if {$numIter == 0} {
	set numIter 100
    } 
    test NormDispIncr 1e-3 $numIter 5
    analyze $numStep
    puts $i
}
