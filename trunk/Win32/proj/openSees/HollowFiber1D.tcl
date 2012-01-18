wipe; # Based on hollow1.tcl
# Changes: 1. transfer into reinforcedSteel Materal model
# Create ModelBuilder 
# ------------------- 

model  BasicBuilder  -ndm  2  -ndf  3 

# Define geometry 
# --------------- 

# Node    tag    xCrd    yCrd 
node       1  +0.000000E+000  +0.000000E+000 
#node       2  +0.000000E+000  +4.000000E+002 
node       3  +0.000000E+000  +1.240000E+002 
#node       4  +0.000000E+000  +12.000000E+002 
node       5  +0.000000E+000  +2.680000E+003 

# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Rz 
fix       1     1     1     1 

# Define nodal masses 
# ------------------- 
# NodeMass.tcl 

# Mass    tag    mx    my    mIz 
#mass       2  0 0 0 
#mass       3  0 0 0
#mass       4  0 0 0
#mass       5  0 0 0

# Define Multi Point Constraints 
# ------------------------------ 
# MPConstraint.tcl 

# Define material(s) 
# ------------------ 
# Materials.tcl 
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

uniaxialMaterial ConcreteL02 2  $fpc -0.03 

#uniaxialMaterial ConfinedConcrete01 $tag  $secType $fpc   $Ec -epscu $epscu  $nu     $L1  $phis  $S   $fyh   $Es0   $haRatio  $mu   $phiLon  -stRatio  $stRatio
#uniaxialMaterial ConfinedConcrete01    2     S1    $fpc  $E  -epscu -0.03    -nu 0.2 100.0 6.0 40.0 300.0 206000.0    0.00   200.0   8.0   -stRatio    0.85


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
#uniaxialMaterial Concrete02 3  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)] $lambda [expr abs(1.*$fT)] [expr abs(1.*$Ets)]
#uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 
uniaxialMaterial  Concrete01 1  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)]

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
#                                               -DMBuck [expr 1.*$lsr2] [expr 1.*$alpha1] \
#                                               -GABuck [expr 1.*$lsr1] [expr 1.*$beta] [expr 1.*$r] [expr 1.*$gama] \

#uniaxialMaterial Concrete01 1 -30.800000000000004 -0.002345151796691333 -24.640000000000004 -0.04 
#uniaxialMaterial Concrete02 2 -36.4 -0.00254944878790797 -7.28 -0.04 0.1 3.64 1820.0
uniaxialMaterial ReinforcingSteel 3 430.0 494.49999999999994 199948.00558280002 19994.80055828 0.008 0.03  -GABuck 6.0 2.0 0.1 0.25  -CMFatigue 0.2 0.5 0.4  -IsoHard 4.3 0.01  -MPCurveParams 0.333 18.0 4.0
#uniaxialMaterial ReinforcingSteel 4 580.0 783.0 199948.00558280002 19994.80055828 0.008 0.016  -DMBuck 6.0 1.0  -CMFatigue 0.2 0.5 0.4  -IsoHard 4.3 0.01  -MPCurveParams 0.333 18.0 4.0
uniaxialMaterial ConcreteL01 5 -40.0 -0.0025
uniaxialMaterial ConcreteL02 6 -40.0 -0.0025
uniaxialMaterial Steel02 7 430.0 199948.00558280002 0.05 18.0 0.925 0.15 0.0 1.0 0.0 1.0 0.0

# Define section(s) 
# ----------------- 
# Sections.tcl 
section Fiber 1 -GJ 1189799672666.6667 { 
layer straight 3 9 50.3 -337.0 -213.0 337.0 -213.0
layer straight 3 0 50.3 -438.0 -162.0 -387.0 -162.0
layer straight 3 2 50.3 387.0 213.0 387.0 162.0
patch quad 1 8 1 375.0 -150.0 -375.0 -150.0 -387.0 -162.0 387.0 -162.0
patch quad 2 6 12 -387.0 162.0 -438.0 162.0 -438.0 -162.0 -387.0 -162.0
layer straight 3 2 50.3 -438.0 213.0 -438.0 162.0
layer straight 3 2 50.3 -387.0 213.0 -387.0 162.0
patch quad 2 6 12 438.0 162.0 387.0 162.0 387.0 -162.0 438.0 -162.0
patch quad 2 16 8 387.0 213.0 -387.0 213.0 -387.0 162.0 387.0 162.0
patch quad 1 16 1 450.0 225.0 -450.0 225.0 -438.0 213.0 438.0 213.0
patch quad 1 24 1 450.0 225.0 438.0 213.0 438.0 -213.0 450.0 -225.0
layer straight 3 3 50.3 -387.0 -112.0 -387.0 112.0
layer straight 3 0 50.3 -387.0 213.0 -438.0 213.0
patch quad 1 24 1 -438.0 213.0 -450.0 225.0 -450.0 -225.0 -438.0 -213.0
layer straight 3 0 50.3 438.0 213.0 387.0 213.0
patch quad 1 1 12 387.0 162.0 375.0 150.0 375.0 -150.0 387.0 -162.0
layer straight 3 2 50.3 387.0 -162.0 387.0 -213.0
layer straight 3 2 50.3 -438.0 -162.0 -438.0 -213.0
layer straight 3 2 50.3 438.0 213.0 438.0 162.0
patch quad 1 1 12 -375.0 150.0 -387.0 162.0 -387.0 -162.0 -375.0 -150.0
patch quad 2 16 12 -387.0 213.0 -438.0 213.0 -438.0 162.0 -387.0 162.0
layer straight 3 2 50.3 -387.0 -162.0 -387.0 -213.0
layer straight 3 3 50.3 438.0 -112.0 438.0 112.0
patch quad 2 16 8 387.0 -162.0 -387.0 -162.0 -387.0 -213.0 387.0 -213.0
layer straight 3 9 50.3 -337.0 -162.0 337.0 -162.0
layer straight 3 3 50.3 387.0 -112.0 387.0 112.0
patch quad 2 6 12 438.0 213.0 387.0 213.0 387.0 162.0 438.0 162.0
layer straight 3 3 50.3 -438.0 -112.0 -438.0 112.0
layer straight 3 0 50.3 -387.0 162.0 -438.0 162.0
layer straight 3 0 50.3 387.0 -213.0 438.0 -213.0
layer straight 3 2 50.3 438.0 -162.0 438.0 -213.0
patch quad 1 8 1 387.0 162.0 -387.0 162.0 -375.0 150.0 375.0 150.0
patch quad 2 16 12 -387.0 -213.0 -387.0 -162.0 -438.0 -162.0 -438.0 -213.0
patch quad 1 16 1 438.0 -213.0 -438.0 -213.0 -450.0 -225.0 450.0 -225.0
layer straight 3 0 50.3 387.0 -162.0 438.0 -162.0
layer straight 3 9 50.3 337.0 213.0 -337.0 213.0
layer straight 3 0 50.3 438.0 162.0 387.0 162.0
layer straight 3 0 50.3 -438.0 -213.0 -387.0 -213.0
layer straight 3 9 50.3 337.0 162.0 -337.0 162.0
patch quad 2 16 12 438.0 -162.0 387.0 -162.0 387.0 -213.0 438.0 -213.0
};
section Fiber 2 -GJ 1085357445780.4637 { 
layer straight 4 2 50.3 -59.6 -225.0 59.6 -225.0
layer straight 4 0 50.3 -155.0 -155.0 -85.0 -155.0
layer straight 4 2 50.3 85.0 225.0 85.0 155.0
patch quad 1 8 1 60.0 -130.0 -60.0 -130.0 -85.0 -155.0 85.0 -155.0
patch quad 2 6 12 -85.0 155.0 -155.0 155.0 -155.0 -155.0 -85.0 -155.0
layer straight 4 2 50.3 -155.0 225.0 -155.0 155.0
layer straight 4 2 50.3 -85.0 225.0 -85.0 155.0
patch quad 2 6 12 155.0 155.0 85.0 155.0 85.0 -155.0 155.0 -155.0
patch quad 2 16 8 85.0 225.0 -85.0 225.0 -85.0 155.0 85.0 155.0
patch quad 1 16 1 180.0 250.0 -180.0 250.0 -155.0 225.0 155.0 225.0
patch quad 1 24 1 180.0 250.0 155.0 225.0 155.0 -225.0 180.0 -250.0
layer straight 4 4 50.3 -85.0 -129.6 -85.0 129.6
layer straight 4 0 50.3 -85.0 225.0 -155.0 225.0
patch quad 1 24 1 -155.0 225.0 -180.0 250.0 -180.0 -250.0 -155.0 -225.0
layer straight 4 0 50.3 155.0 225.0 85.0 225.0
patch quad 1 1 12 85.0 155.0 60.0 130.0 60.0 -130.0 85.0 -155.0
layer straight 4 2 50.3 85.0 -155.0 85.0 -225.0
layer straight 4 2 50.3 -155.0 -155.0 -155.0 -225.0
layer straight 4 2 50.3 155.0 225.0 155.0 155.0
patch quad 1 1 12 -60.0 130.0 -85.0 155.0 -85.0 -155.0 -60.0 -130.0
patch quad 2 16 12 -85.0 225.0 -155.0 225.0 -155.0 155.0 -85.0 155.0
layer straight 4 2 50.3 -85.0 -155.0 -85.0 -225.0
layer straight 4 4 50.3 155.0 -129.6 155.0 129.6
patch quad 2 16 8 85.0 -155.0 -85.0 -155.0 -85.0 -225.0 85.0 -225.0
layer straight 4 2 50.3 -59.6 -155.0 59.6 -155.0
layer straight 4 4 50.3 85.0 -129.6 85.0 129.6
patch quad 2 6 12 155.0 225.0 85.0 225.0 85.0 155.0 155.0 155.0
layer straight 4 4 50.3 -155.0 -129.6 -155.0 129.6
layer straight 4 0 50.3 -85.0 155.0 -155.0 155.0
layer straight 4 0 50.3 85.0 -225.0 155.0 -225.0
layer straight 4 2 50.3 155.0 -155.0 155.0 -225.0
patch quad 1 8 1 85.0 155.0 -85.0 155.0 -60.0 130.0 60.0 130.0
patch quad 2 16 12 -85.0 -225.0 -85.0 -155.0 -155.0 -155.0 -155.0 -225.0
patch quad 1 16 1 155.0 -225.0 -155.0 -225.0 -180.0 -250.0 180.0 -250.0
layer straight 4 0 50.3 85.0 -155.0 155.0 -155.0
layer straight 4 2 50.3 59.6 225.0 -59.6 225.0
layer straight 4 0 50.3 155.0 155.0 85.0 155.0
layer straight 4 0 50.3 -155.0 -225.0 -85.0 -225.0
layer straight 4 2 50.3 59.6 155.0 -59.6 155.0
patch quad 2 16 12 155.0 -155.0 85.0 -155.0 85.0 -225.0 155.0 -225.0
};
# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran.tcl Corotational PDelta Linear

# GeoTran    type    tag 
set transfTag 1;
geomTransf  Corotational       $transfTag 

# Define element(s) 
# ----------------- 
#Priestley96
set ksi 6.8972285944571905; set fY [expr 60*$ksi]; set Db 1.0; set ElementLength  8.000000E+002; #1.600000E+003;
set E 3e4; set A 630000; set Iz 5146875000; set Iy 16790625000; set G 1e8; set J 1e6;

set Lp [expr 0.08*$ElementLength+0.15*($fY/$ksi)*$Db]; set Lp2 [expr 0.3*($fY/$ksi)*$Db]; if {$Lp < $Lp2} {; set Lp $Lp2;};
set Lpi [expr $Lp]; 

# dispBeamColumn or nonlinearBeamColumn
# Element "Beam":            eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       1       1       3     2     1     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element beamWithHinges $eleTag $iNode $jNode $secTagI $Lpi $secTagJ $Lpj $E $A $Iz $Iy $G $J $transfTag <-mass $massDens> <-iter $maxIters $tol>
#element beamWithHinges        1       1       3     4    $Lpi   4 $Lpi $E $A $Iz $Iy $G $J 1 -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element beamWithHinges        1       1       3     4     $Lpi  4 1 $E $A $Iz $transfTag -mass +4.665600E+000  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       2       3       5     2     1     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element beamWithHinges        2       3       5     4    $Lpi   4 $Lpi $E $A $Iz $Iy $G $J 1 -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element beamWithHinges        2       3       5     4    1   4 1 $E $A $Iz $transfTag -mass +4.665600E+000  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  dispBeamColumn       3       3       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  dispBeamColumn       4       4       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 


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
    load       5  +0.000000E+000  -2.500000E+005  +0.000000E+000 
 
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

# Reset for next analysis case 
# ---------------------------- 
loadConst -time 0.0 
remove recorders 
wipeAnalysis 

# Define load pattern 
# ------------------- 
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       2  "Linear" { 
    # Load    nodeTag    LoadValues 
    load       5  1.000000E+005  +0.000000E+000   +0.000000E+000 
 
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    xL    <Px> 
} 

# Creating a recorder to monitor nodal displacements
#recorder Element -file hollow1_Col_1_se.dat -time -ele 1 section fiber 0. 150.0 stressStrain
#recorder Element -file hollow1_ele_p1_S.dat -time -ele 1 2 3 4 integrPoint 1 stresses
#recorder Element -file hollow1_ele_p1_E.dat -time -ele 1 2 3 4 integrPoint 1 strains

recorder Node -file hollow2_nodeReaction.dat -time -node 1 -dof 1 2 3 reaction
recorder Node -file hollow2_node_disp_Output.dat -time -node 5 -dof 1 2 disp;

# ----------- set up analysis parameters
# characteristics of pushover analysis
set Dincr [expr 1.0];	# displacement increment. you want this to be small, but not too small to slow analysis
set IDctrlNode 5;
set IDctrlDOF 1;

# --------------------------------------------------------------------------------------------------
# static analysis parameters

# CONSTRAINTS handler
#          Plain Constraints (only for homogeneous equations)
#          Lagrange Multipliers
#          Penalty Method (rigidDiaphragm)
#          Transformation Method
set RigidDiaphragm "OFF";
variable constraintsTypeStatic Plain;  # default;
if {  [info exists RigidDiaphragm] == 1} {
	if {$RigidDiaphragm=="ON"} {
		variable constraintsTypeStatic Lagrange;	#     for large model, try Transformation
	};	# if rigid diaphragm is on
};	# if rigid diaphragm exists
constraints $constraintsTypeStatic

# DOF NUMBERER
#          Plain
#          RCM
set numbererTypeStatic Plain
numberer $numbererTypeStatic 

# SYSTEM
#          ProfileSPD -- Direct profile solver for symmetric positive definite matrices 
#          BandGeneral -- Direct solver for banded unsymmetric matrices 
#          BandSPD -- Direct solver for banded symmetric positive definite matrices 
#          SparseGeneral -- Direct solver for unsymmetric sparse matrices 
#          SparseSPD -- Direct solver for symmetric sparse matrices 
#          UmfPack -- Direct UmfPack solver for unsymmetric matrices 
#          FullGeneral
set systemTypeStatic BandGeneral;		# try UmfPack for large model  
system $systemTypeStatic 

# TEST: # convergence test to 
#          NormUnbalance -- Specifies a tolerance on the norm of the unbalanced load at the current iteration 
#          NormDispIncr -- Specifies a tolerance on the norm of the displacement increments at the current iteration 
#          EnergyIncr-- Specifies a tolerance on the inner product of the unbalanced load and displacement increments at the current iteration 
#          RelativeNormUnbalance --
#          RelativeNormDispIncr --
#          RelativeEnergyIncr --
variable TolStatic 1.0e-3;                        # Convergence Test: tolerance
variable Tol [expr $TolStatic];
variable maxNumIterStatic 1000;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
variable printFlagStatic 5;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
variable testTypeStatic NormDispIncr ;	# Convergence-test type: NormDispIncr, EnergyIncr
test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;
# for improved-convergence procedure:
variable maxNumIterConvergeStatic 2000;	
variable printFlagConvergeStatic 0;


# Solution ALGORITHM: 
#          Linear -- Uses the solution at the first iteration and continues 
#          Newton -- Uses the tangent at the current iteration to iterate to convergence 
#          ModifiedNewton -- Uses the tangent at the first iteration to iterate to convergence 
#          NewtonLineSearch -- 
#          KrylovNewton -- 
#          BFGS -- 
#          Broyden -- 
variable algorithmTypeStatic KrylovNewton
algorithm $algorithmTypeStatic;        

# Static INTEGRATOR:
#          LoadControl -- Specifies the incremental load factor to be applied to the loads in the domain 
#          DisplacementControl -- Specifies the incremental displacement at a specified DOF in the domain 
#          Minimum Unbalanced Displacement Norm -- Specifies the incremental load factor such that the residual displacement norm in minimized 
#          Arc Length -- Specifies the incremental arc-length of the load-displacement path 
# Transient INTEGRATOR: -- determine the next time step for an analysis including inertial effects 
#          Newmark -- The two parameter time-stepping method developed by Newmark 
#          HHT -- The three parameter Hilbert-Hughes-Taylor time-stepping method 
#          Central Difference -- Approximates velocity and acceleration by centered finite differences of displacement 
integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr
#integrator ArcLength $Dincr [expr $Dincr/10.0]

# ANALYSIS
#          Static Analysis -- solves the KU=R problem, without the mass or damping matrices. 
#          Transient Analysis -- solves the time-dependent analysis. The time step in this type of analysis is constant. The time step in the output is also constant. 
#          variableTransient Analysis -- performs the same analysis type as the Transient Analysis object. The time step, however, is variable. This method is used when 
#                 there are convergence problems with the Transient Analysis object at a peak or when the time step is too small. The time step in the output is also variable.
set analysisTypeStatic Static
analysis $analysisTypeStatic 

# initialize in case we need to do an initial stiffness iteration
initialize

# Perform the Cycle pushover analysis
# characteristics of cyclic analysis	
set iDmax {5 10 15 20 25 30 35 40 45 50 60 70 80 100 120 140 160} ; #1 3 5 10 4 14 17 7 25 30 33 40#20 25 30 35 40 45 50 60 70 80 100 120 140 160 240 };  # vector of displacement-cycle peaks, in terms of storey drift ratio
set CycleType Full;    # you can do Full / Push / Half cycles with the proc
set Ncycles 1;				# specify the number of cycles at each peak
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