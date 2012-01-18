# Start of model generation 
# ========================= 
wipe
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  3  -ndf  6 

# Define geometry 
# --------------- 
# NodeCoord.tcl 

# Node    tag    xCrd    yCrd    zCrd 
node       1  +0.000000E+000  +0.000000E+000  +0.000000E+000 
node       2  +1.000000E+003  +0.000000E+000  +0.000000E+000 
node       3  +0.000000E+000  +0.000000E+000  +1.000000E+003 
node       4  +1.000000E+003  +0.000000E+000  +1.000000E+003 

# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Dz    Rx    Ry    Rz 
fix       1     1     1     1     1     1     1 
fix       2     1     1     1     1     1     1 

# Define nodal masses 
# ------------------- 
# NodeMass.tcl 

# Node    tag    mx    my    mz    mIx    mIy    mIz 
mass       1  +2.000000E+002  +2.000000E+002  +2.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       2  +2.000000E+002  +2.000000E+002  +2.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       3  +2.000000E+002  +2.000000E+002  +2.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       4  +2.000000E+002  +2.000000E+002  +2.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 

# Define Multi Point Constraints 
# ------------------------------ 
# MPConstraint.tcl 

# Define material(s) 
# ------------------ 
# Materials.tcl 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.900000E+004  +0.000000E+000 

# Material "CoreConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       2  -3.500000E+001  -2.000000E-003  -4.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "CoverConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "reinforcedBar":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
uniaxialMaterial  Steel02       4  +3.750000E+002  +2.000000E+005  +5.000000E-002  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "Concrete02C20":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       5  -2.000000E+001  -4.000000E-003  -2.800000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +2.400000E+004 

# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D       6  +2.000000E+005  +2.500000E-001  +0.000000E+000 

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 23.0;
set wfy 325.2;
set wE  207544.0;
set rou1 0.0023;
set rou2 0.0023;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
#uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
#uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2


# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteZ02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteZ02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
uniaxialMaterial Concrete09 13 -2.000000E+001  -4.000000E-003  -2.800000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +2.400000E+004 
uniaxialMaterial Concrete09 14 -2.000000E+001  -4.000000E-003  -2.800000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +2.400000E+004 

set pi [expr 2.0*asin(1.0)]
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMAterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
nDMaterial FAReinforcedConcretePlateFiber  15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
 

# Define section(s) 
# ----------------- 
# Sections.tcl 

# Section "ElasticDefault":    secTag    E    A    Iz  Iy  G  J 
section  Elastic       1  +2.900000E+004  +1.800000E+002  +4.860000E+003  +1.500000E+003  +1.115400E+004  +3.916000E+003 

# Section "Hollow01":    secTag 
section  Fiber       2  { 
    # PatchBox "Patch01":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2    30     2  +2.700000E+002  +1.800000E+002  +2.700000E+002  -1.800000E+002  +3.000000E+002  -1.800000E+002  +3.000000E+002  +1.800000E+002 
    patch  quad       2     2    18  -2.700000E+002  +1.800000E+002  -2.700000E+002  +1.500000E+002  +2.700000E+002  +1.500000E+002  +2.700000E+002  +1.800000E+002 
    patch  quad       2     2    18  -2.700000E+002  -1.500000E+002  -2.700000E+002  -1.800000E+002  +2.700000E+002  -1.800000E+002  +2.700000E+002  -1.500000E+002 
    patch  quad       2    30     2  -3.000000E+002  +1.800000E+002  -3.000000E+002  -1.800000E+002  -2.700000E+002  -1.800000E+002  -2.700000E+002  +1.800000E+002 
    # PatchBox "Patch02":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2    27     4  +2.100000E+002  +1.500000E+002  +2.100000E+002  -1.500000E+002  +2.700000E+002  -1.500000E+002  +2.700000E+002  +1.500000E+002 
    patch  quad       2     4    15  -2.100000E+002  +1.500000E+002  -2.100000E+002  +9.000000E+001  +2.100000E+002  +9.000000E+001  +2.100000E+002  +1.500000E+002 
    patch  quad       2     4    15  -2.100000E+002  -9.000000E+001  -2.100000E+002  -1.500000E+002  +2.100000E+002  -1.500000E+002  +2.100000E+002  -9.000000E+001 
    patch  quad       2    27     4  -2.700000E+002  +1.500000E+002  -2.700000E+002  -1.500000E+002  -2.100000E+002  -1.500000E+002  -2.100000E+002  +1.500000E+002 
    # PatchBox "Patch03":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3    22     2  +2.250000E+002  +1.350000E+002  +2.250000E+002  -1.350000E+002  +2.550000E+002  -1.350000E+002  +2.550000E+002  +1.350000E+002 
    patch  quad       3     2    12  -2.250000E+002  +1.350000E+002  -2.250000E+002  +1.050000E+002  +2.250000E+002  +1.050000E+002  +2.250000E+002  +1.350000E+002 
    patch  quad       3     2    12  -2.250000E+002  -1.050000E+002  -2.250000E+002  -1.350000E+002  +2.250000E+002  -1.350000E+002  +2.250000E+002  -1.050000E+002 
    patch  quad       3    22     2  -2.550000E+002  +1.350000E+002  -2.550000E+002  -1.350000E+002  -2.250000E+002  -1.350000E+002  -2.250000E+002  +1.350000E+002 
    # LayerStraight "Layer01":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +1.500000E+002  -2.700000E+002  +1.500000E+002 
    # LayerStraight "Layer02":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +9.000000E+001  +2.700000E+002  +9.000000E+001 
    # LayerStraight "Layer03":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  -9.000000E+001  -2.700000E+002  -9.000000E+001 
    # LayerStraight "Layer04":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  -1.500000E+002  -2.700000E+002  -1.500000E+002 
    # LayerStraight "Layer05":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  +2.700000E+002  +3.000000E+001  +2.700000E+002  -3.000000E+001 
    # LayerStraight "Layer06":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  +2.100000E+002  +3.000000E+001  +2.100000E+002  -3.000000E+001 
    # LayerStraight "Layer07":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  -2.100000E+002  +3.000000E+001  -2.100000E+002  -3.000000E+001 
    # LayerStraight "Layer08":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  -2.700000E+002  +3.000000E+001  -2.700000E+002  -3.000000E+001 
} 

# Section "PlateFiber":   secTag  ndMatTag  h
section  PlateFiber       3    6  20


# Section "ElaMembranePlateSec":    secTag    E    v    h    rho 
section  ElasticMembranePlateSection       4  +2.000000E+005  +2.500000E-001  +1.000000E-001  0.0
#section PlateFiber secTag  ndMatTag  h
section  PlateFiber  5   15   240
# Define geometric transformation(s) 
# ---------------------------------- 

# Define element(s) 
# ----------------- 
# Elements.tcl 

# Element "MITC4":    eleTag    NodeI    NodeJ    NodeK    NodeL    secTag 
element  shell02       1       1       2       4       3       5 

# Define time series 
# ------------------ 
# TimeSeries.tcl 

# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# TimeSeries "TimeSeriesX":    tsTag    dt    filePath    cFactor 
timeSeries  Path       2  -dt  +1.000000E-002  -filePath  TimeSeriesX.thf  -factor  +1.000000E-003 

# TimeSeries "TimeSeriesY":    tsTag    dt    filePath    cFactor 
timeSeries  Path       3  -dt  +1.000000E-002  -filePath  TimeSeriesY.thf  -factor  +1.000000E-003 


# Start of anaysis generation 
# =========================== 

# Get Initial Stiffness 
# --------------------- 
initialize 

# Analysis: StaticDefaultCase 
# +++++++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    load       3  1.000000E+000   +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load       4  1.000000E+000   +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
 
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 

# Define recorder(s) 
# -------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  4 -dof  1  2  3    disp 
recorder  Element -file PushoverCase_Ele_All.force -time -ele 1 globalforce
recorder  Element -file PushoverCase_Ele.stress -time -ele 1 stresses
recorder  Element -file PushoverCase_Mat_G1.stress -time -ele 1 -material 1 stresses
recorder  Element -file PushoverCase_Mat_G2.stress -time -ele 1 -material 2 stresses
recorder  Element -file PushoverCase_Mat_G3.stress -time -ele 1 -material 3 stresses
recorder  Element -file PushoverCase_Mat_G4.stress -time -ele 1 -material 4 stresses

# Define analysis options 
# ----------------------- 
# AnalysisOptn_1.tcl 

# AnalysisOptn "StaticDefault": Type: Static 
# ------------------------------------------ 
# Constraint Handler 
constraints  Plain 
# Convergence Test 
test  NormDispIncr  +1.000000E-002    25     0     2 
# Integrator 
integrator  DisplacementControl 4 1 +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

# perform the analysis  0  1  2   3     4    5    6   7    8   9   10    11   12   13    14   15   16   17
#                1.0   -1.0 +2.0  -2.0 +3.0 -3.0 4.0 -4.0 5.0 -5.0 -6.0  6.0 -7.0  8.0 -10.0 15.0 -25.0 40.0  
set numSteps   {  10     40   30    40  50    60  70   80  90  100  110  120  130  150   180  250   400  650 }; #   800 1200  1600  2050  2500 2960  3470 1000    10}  
set numIters   { 200    300  400   400  450  450 500  550 550  550  550  550  550  550   550  550   550  550 }; #   100  100   100   100   100  500   500  500   300}  
set increments { 0.1  -0.05  0.1  -0.1  0.1 -0.1 0.1 -0.1 0.1 -0.1 -0.1  0.1 -0.1  0.1  -0.1  0.1  -0.1  0.1 }; # -0.01 0.01 -0.01  0.01 -0.01 0.01 -0.01 0.01 -0.01}  

for { set i 0 } { $i<18 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set increment [lindex $increments $i]

    integrator DisplacementControl 4 1 $increment
    if {$numIter == 0} {
	set numIter 100
    } 
    test NormDispIncr 1e-3 $numIter 5
    analyze $numStep
    puts $i
}
remove recorders;
