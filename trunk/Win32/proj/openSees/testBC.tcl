wipe; # Based on hollow1.tcl
# Changes: 1. transfer into reinforcedSteel Materal model
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  2  -ndf  3 

# Define geometry 
# --------------- 

# Node    tag    xCrd    yCrd 
node       1  +0.000000E+000  +0.000000E+000 
node       2  +0.000000E+000  +3.000000E+002 

# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Rz 
fix       1     1     1     1 

# Define nodal masses 
# ------------------- 
# NodeMass.tcl 

# Mass    tag    mx    my    mIz 
#mass       2  1000 1000 5000 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.900000E+004  +0.000000E+000 

# Material "CoreConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
# C30
set Fc 20.1
set rFc 1.3;
set rFu 0.2;
set E [expr 4733.*sqrt($rFc*abs($Fc))]; #for psi: [expr 57000.*sqrt($rFc*abs($Fc))];
set fpc [expr -$rFc*abs($Fc)];
set epsc0 [expr 2.*$fpc/$E];
set fpcU [expr -$rFu*$rFc*abs($Fc)];
set epsU -0.03;
set lambda 0.1;
set fT  +0.1*$rFc*abs($Fc);
set Ets  $fT/0.002;

#uniaxialMaterial Concrete02 2  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)] $lambda [expr abs(1.*$fT)] [expr abs(1.*$Ets)];
#uniaxialMaterial  Concrete02         -3.900000E+001  -2.000000E-003  -4.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

#uniaxialMaterial ConfinedConcrete01 $tag  $secType $fpc   $Ec -epscu $epscu  $nu     $L1  $phis  $S   $fyh   $Es0   $haRatio  $mu   $phiLon  -stRatio  $stRatio
uniaxialMaterial ConfinedConcrete01    2     S1    $fpc  $E  -epscu -0.03    -nu 0.2 100.0 6.0 40.0 300.0 206000.0    0.00   200.0   8.0   -stRatio    0.85


# Material "CoverConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
set Fc 18.3;
set rFc 1.0;
set rFu 0.2;
set E [expr 4733.*sqrt($rFc*abs($Fc))]; #for psi: [expr 57000.*sqrt($rFc*abs($Fc))];
set fpc [expr -$rFc*abs($Fc)];
set epsc0 [expr 2.*$fpc/$E]
set fpcU [expr -$rFu*$rFc*abs($Fc)]
set epsU -0.03;
set lambda 0.1;
set fT  +0.1*$rFc*abs($Fc);
set Ets  $fT/0.002;
uniaxialMaterial Concrete02 3  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)] $lambda [expr abs(1.*$fT)] [expr abs(1.*$Ets)]
#uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 
#uniaxialMaterial  Concrete01 3  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)]

# Material "reinforcedBar":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
#uniaxialMaterial  Steel02       5  +3.100000E+002  +2.000000E+005  +2.000000E-002  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "reinforcedSteel360":    matTag    Fy    Fu    E    Eh    epsh    epsu    <-GABuck  lsr  beta  r  gamma>    <-DMBuck  lsr  <alpha>>    <-CMFatigue  Cf  alpha  Cd>    <-IsoHard  <a1  <limit>>>    <-MPCurveParams  R1  R2  R3> 
#uniaxialMaterial  ReinforcingSteel       4  +3.600000E+002  +3.350000E+002  +2.000000E+005  +5.000000E+002  +5.000000E-002  +1.200000E-001  -GABuck  +6.000000E+000  +1.000000E+000  +4.000000E-001  +5.000000E-001  -DMBuck  +6.000000E+000  +1.000000E+000  -CMFatigue  +2.600000E-001  +5.060000E-001  +3.890000E-001  -IsoHard  +4.300000E+000  +1.000000E-002  -MPCurveParams  +3.330000E-001  +1.800000E+001  +4.000000E+000 
set MaterialTag 4;
set fY 300.;
set fU [expr 1.3*$fY];
set Es 2.0e5;
set Esh [expr 0.1*$Es];
set esh 0.008;
set eult 0.016;
set lsr1 5.0;
set beta 2.0;
set r 0.0; #0.4;
set gama 0.5;
set lsr2 8.0
set alpha1 0.95
set Cf 0.26;
set alpha2 0.506;
set Cd 0.3; #0.389;
set a1 4.3;
set limit 0.2; #0.01;
set R1 0.333;
set R2 8;
set R3 4;

uniaxialMaterial ReinforcingSteel $MaterialTag [expr 1.*$fY] [expr 1.*$fU] [expr 1.*$Es] [expr 1.*$Esh] [expr 1.*$esh] [expr 1.*$eult] \
																	-GABuck [expr 1.*$lsr1] [expr 1.*$beta] [expr 1.*$r] [expr 1.*$gama] \
																	-CMFatigue [expr 1.*$Cf] [expr 1.*$alpha2] [expr 1.*$Cd] \
																	-IsoHard [expr 1.*$a1] [expr 1.*$limit] \
																	-MPCurveParams [expr 1.*$R1] [expr 1.*$R2] [expr 1.*$R3]
#																	-DMBuck [expr 1.*$lsr2] [expr 1.*$alpha1] \
#																	-GABuck [expr 1.*$lsr1] [expr 1.*$beta] [expr 1.*$r] [expr 1.*$gama] \

set wfc 28.0;
set wfyv 560;
set wfyh1 450;
set wE 190000.0;
set rou1 0;
set rou2 0;
set rouv 0.0246;
set rouh1 0.00246; # #10@70

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

# Define section(s) 
# ----------------- 
# Sections.tcl 

# Section "ElasticDefault":    secTag    E    A    Iz 
section  Elastic       1  +2.900000E+004  +1.800000E+002  +4.860000E+003 

# Section "Hollow01":    secTag 
section  Fiber       2  { 
    # PatchBox "Patch01":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3    12     1  +2.700000E+002  +1.800000E+002  +2.700000E+002  -1.800000E+002  +3.000000E+002  -1.800000E+002  +3.000000E+002  +1.800000E+002 
    patch  quad       3     1    18  -2.700000E+002  +1.800000E+002  -2.700000E+002  +1.500000E+002  +2.700000E+002  +1.500000E+002  +2.700000E+002  +1.800000E+002 
    patch  quad       3     1    18  -2.700000E+002  -1.500000E+002  -2.700000E+002  -1.800000E+002  +2.700000E+002  -1.800000E+002  +2.700000E+002  -1.500000E+002 
    patch  quad       3    12     1  -3.000000E+002  +1.800000E+002  -3.000000E+002  -1.800000E+002  -2.700000E+002  -1.800000E+002  -2.700000E+002  +1.800000E+002 
    # PatchBox "Patch02":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2    10     2  +2.100000E+002  +1.500000E+002  +2.100000E+002  -1.500000E+002  +2.700000E+002  -1.500000E+002  +2.700000E+002  +1.500000E+002 
    patch  quad       2     2    14  -2.100000E+002  +1.500000E+002  -2.100000E+002  +9.000000E+001  +2.100000E+002  +9.000000E+001  +2.100000E+002  +1.500000E+002 
    patch  quad       2     2    14  -2.100000E+002  -9.000000E+001  -2.100000E+002  -1.500000E+002  +2.100000E+002  -1.500000E+002  +2.100000E+002  -9.000000E+001 
    patch  quad       2    10     2  -2.700000E+002  +1.500000E+002  -2.700000E+002  -1.500000E+002  -2.100000E+002  -1.500000E+002  -2.100000E+002  +1.500000E+002 
    # PatchBox "Patch03":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3     6     1  +1.800000E+002  +9.000000E+001  +1.800000E+002  -9.000000E+001  +2.100000E+002  -9.000000E+001  +2.100000E+002  +9.000000E+001 
    patch  quad       3     2    12  -1.800000E+002  +9.000000E+001  -1.800000E+002  +6.000000E+001  +1.800000E+002  +6.000000E+001  +1.800000E+002  +9.000000E+001 
    patch  quad       3     2    12  -1.800000E+002  -6.000000E+001  -1.800000E+002  -9.000000E+001  +1.800000E+002  -9.000000E+001  +1.800000E+002  -6.000000E+001 
    patch  quad       3     6     1  -2.100000E+002  +9.000000E+001  -2.100000E+002  -9.000000E+001  -1.800000E+002  -9.000000E+001  -1.800000E+002  +9.000000E+001 
    # LayerStraight "Layer01":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +1.500000E+002  -2.700000E+002  +1.500000E+002 
    # LayerStraight "Layer02":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +9.000000E+001  -2.700000E+002  +9.000000E+001 
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

set t1 35.0;
set NStrip1 2; # thickness 1
set t2 80.0;
set NStrip2 4; # thickness 2
set t3 42.0;
set NStrip3 3; # thickness 3
set np 1; # int. points
set C 0.4; # center of rotation 

#section definition 
#section CSMMFiber2d 3 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { 
section FiberInt 3 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { 
	fiber -40 0 110 1;
	fiber -30 0 120 1;
	fiber -20 0 210 2;
	fiber -10 0 220 2;
	fiber   0 0 230 2;
	fiber  10 0 240 2;
	fiber  20 0 310 3;
	fiber  30 0 320 3;
	fiber  40 0 330 3;
}

# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran.tcl Corotational PDelta Linear

# GeoTran    type    tag 
geomTransf  LinearInt       1 

# Define element(s) 
# ----------------- 
# Elements.tcl 
# dispBeamColumn or nonlinearBeamColumn
# Element "Beam":            eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  nonlinearBeamColumn       1       1       2     $np     3     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
element dispBeamColumnInt 1 1 2 $np 3 1 $C

# Define time series 
# ------------------ 
# TimeSeries.tcl 

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
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1  "Linear" { 
    # Load    nodeTag    LoadValues 
    load       2  +0.000000E+000  -2.800000E+005  +0.000000E+000 
 
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    xL    <Px> 
} 

# Define recorder(s) 
# -------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  Node_DefoShape_Dsp.out  -time -nodeRange 1  5 -dof  1  2  3  disp 

# Define analysis options 
# ----------------------- 
# AnalysisOptn_3.tcl 

# AnalysisOptn "StaticDefault": Type: Static 
# ------------------------------------------ 
# Constraint Handler 
constraints  Plain 
# Convergence Test 
test  EnergyIncr  +1.000000E-004    25     0     2 
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

analyze     1 