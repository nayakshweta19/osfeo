# quad3D.tcl 
wipe;
source Units&Constants.tcl;
######################## 
# Analysis-Sequence  1 #
######################## 

# Start of model generation 
# ========================= 

# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  3  -ndf  6 

# Define geometry 
# --------------- 
# NodeCoord.tcl 

# Node    tag    xCrd    yCrd    zCrd 
node       1  +0.000000E+000  +0.000000E+000  +0.000000E+000 
node       2  +0.000000E+000  +0.000000E+000  +1.200000E+003 
node       3  +8.000000E+002  +0.000000E+000  +0.000000E+000 
node       4  +8.000000E+002  +0.000000E+000  +1.200000E+003 
node       5  +0.000000E+000  +6.000000E+002  +0.000000E+000 
node       6  +0.000000E+000  +6.000000E+002  +1.200000E+003 
node       7  +8.000000E+002  +6.000000E+002  +0.000000E+000 
node       8  +8.000000E+002  +6.000000E+002  +1.200000E+003 
node       9  +0.000000E+000  +0.000000E+000  +0.000000E+000 
node      10  +0.000000E+000  +0.000000E+000  +1.200000E+003 
node      11  +8.000000E+002  +0.000000E+000  +0.000000E+000 
node      12  +8.000000E+002  +0.000000E+000  +1.200000E+003 
node      13  +0.000000E+000  +6.000000E+002  +0.000000E+000 
node      14  +0.000000E+000  +6.000000E+002  +1.200000E+003 
node      15  +8.000000E+002  +6.000000E+002  +0.000000E+000 
node      16  +8.000000E+002  +6.000000E+002  +1.200000E+003 


# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Dz    Rx    Ry    Rz 
fix       1     1     1     1     1     1     1 
fix       3     1     1     1     1     1     1 
fix       5     1     1     1     1     1     1 
fix       7     1     1     1     1     1     1 
fix       9     1     1     1     1     1     1 
fix      11     1     1     1     1     1     1 
fix      13     1     1     1     1     1     1 
fix      15     1     1     1     1     1     1 

# Define nodal masses 
# ------------------- 
# NodeMass.tcl 

# Node    tag    mx    my    mz    mIx    mIy    mIz 
mass       2  +3.000000E+002  +3.000000E+002  +3.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       4  +3.000000E+002  +3.000000E+002  +3.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       6  +3.000000E+002  +3.000000E+002  +3.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 
mass       8  +3.000000E+002  +3.000000E+002  +3.000000E+002  +1.000000E+000  +1.000000E+000  +1.000000E+000 

# Define Multi Point Constraints 
# ------------------------------ 
# MPConstraint.tcl 

# Equal DOF: MP01:    mNodeTag    sNodeTag    dof 
equalDOF       2      10  1  2  3  4  5  6 

# Equal DOF: MP03:    mNodeTag    sNodeTag    dof 
equalDOF       6      14  1  2  3  4  5  6 

# Equal DOF: MP04:    mNodeTag    sNodeTag    dof 
equalDOF       8      16  1  2  3  4  5  6 

# Equal DOF: MP2:    mNodeTag    sNodeTag    dof 
equalDOF       4      12  1  2  3  4  5  6 

# Define material(s) 
# ------------------ 
# Materials.tcl 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.900000E+004  +0.000000E+000 

# Material "CoverConcreteC35":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       2  -2.340000E+007  -2.000000E-003  -4.680000E+006  -3.343100E-003  +1.000000E-001  +2.200000E+006  +1.000000E+020 

# Material "CoreConcreteC35":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       3  -3.280000E+007  -2.803000E-003  -6.560000E+006  -1.205000E-001  +1.000000E-001  +2.200000E+006  +5.418947E+009 

# Material "CoverConcreteC30":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       4  -2.010000E+007  -2.000000E-003  -4.020000E+006  -4.178600E-003  +1.000000E-001  +2.010000E+006  +1.000000E+020 

# Material "CoreConcreteC30":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       5  -2.950000E+007  -2.935300E-003  -5.900000E+006  -1.107800E-001  +1.000000E-001  +2.010000E+006  +5.025000E+009 

# Material "SteelHRB335":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
uniaxialMaterial  Steel02       6  +3.350000E+008  +2.000000E+005  +5.000000E-002  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "ElasticTorsion":    matTag    E    eta  
uniaxialMaterial  Elastic       7  +1.000000E+012  +0.000000E+000 


# Define section(s) 
# ----------------- 
source RCsection.tcl
# Section "ElasticDefault":    secTag    E    A    Iz  Iy  G  J 
section  Elastic       1  +2.900000E+004  +1.800000E+002  +4.860000E+003  +1.500000E+003  +1.115400E+004  +3.916000E+003 

# ===============
# # Section "Beam300X500":    secTag Fiber section 1
# ===============
# Locaton:1-5th floor
# 0.5m*0.5m
# 3 bar with diameter 25 each said

set  numBars   3
set  tatolArea [expr 1472*pow($mm,2)]  
set  barArea   [expr  $tatolArea/$numBars]
#               id h   b    numBars  barArea coreID coverID steelID  fc  E
RcColumnSection 2 0.5 0.5 $numBars $barArea 3 2 6 23400000.0 23400000000 

# Section "GJ":    secTag    E    A    Iz  Iy  G  J 
section  Elastic       3  +2.000000E+005  +1.000000E+000  +1.000000E+000  +1.000000E+000  +6.000000E+004  +1.000000E+000 


 # Section "Beam300X500t":    secTag    matTag    string ...    -section secTag 
section  Aggregator       4  7 T  -section 2


# =================
# Section "Column500X500":    secTag # Fiber section 5
# =================
# Location:all
# 0.5m*0.3m
# Up   said with 3 bar diameter 20
# Down said with 3 bar diameter 20
set  numBarsUp1   3
set  tatolArea    [expr  942*pow($mm,2)]  
set  barAreaUp1   [expr  $tatolArea/$numBarsUp1]
set  numBarsDown  3  
set  tatolArea    [expr  942*pow($mm,2)]  
set  barAreaDown  [expr  $tatolArea/$numBarsDown]
set  Lo           0.3;#beam span
set  BfIn         [expr $Lo];#the wide of beam flange   
#             id h b numBarsUp1 barAreaUp2 numBarsDown barAreaDown coreID coverID steelID {shape T} Bf fc E 
RcBeamSection 5 0.5 0.3  $numBarsUp1 $barAreaUp1 $numBarsDown  $barAreaDown 3 2 6 T $BfIn 23400000.0 23400000000  


 # Section "Column500X500t":    secTag    matTag    string ...    -section secTag 
section  Aggregator       6  7 T  -section 5


# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran    type    tag    vec_xz 
geomTransf  Linear       1  1  0  0 
# GeoTran    type    tag    vec_xz 
geomTransf  Linear       2  0  1  0 
# GeoTran    type    tag    vec_xz 
geomTransf  Linear       3  +0.000000E+000  -1.000000E+000  +0.000000E+000 
# GeoTran    type    tag    vec_xz 
geomTransf  Linear       4  +0.000000E+000  -1.000000E+000  +0.000000E+000 

# Define element(s) 
# ----------------- 
# Elements.tcl 

# Element "ColumnEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       1       1       2     5     4     2  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "ColumnEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       2       3       4     5     4     2  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "BeamEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       3       2       4     5     6     3  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "ColumnEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       4       5       6     5     4     2  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "ColumnEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       5       7       8     5     4     2  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "BeamEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       6       6       8     5     6     4  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "BeamEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       7       2       6     5     6     1  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Element "BeamEle":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  nonlinearBeamColumn       8       4       8     5     6     1  -mass +0.000000E+000  -iter   10  +1.000000E-008 

# Material "ElasticIso":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D       8  +2.000000E+50  +2.000000E-001  +1.000000E+003 

# Element "quad":    eleTag    NodeI    NodeJ    NodeK    NodeL    h    type    matTag    pressure    rho    b1    b2 
element  quad3D       9       9      11      12      10  +1.500000E+001  "PlaneStress"       8  +0.000000E+000  +0.000000E+000  +0.000000E+000 +0.000000E+000 

# Element "quad":    eleTag    NodeI    NodeJ    NodeK    NodeL    h    type    matTag    pressure    rho    b1    b2 
element  quad3D      10       7       5       6       8  +1.500000E+001  "PlaneStress"       8  +0.000000E+000  +0.000000E+000  +0.000000E+000 +0.000000E+000 

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
# LoadPattern_2.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    load       2  +1.500000E+003  +2.500000E+003  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load       4  +1.500000E+003  +2.500000E+003  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load       6  +1.500000E+003  +2.500000E+003  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load       8  +1.500000E+003  +2.500000E+003  +0.000000E+000  +0.000000E+000  +0.000000E+000  +0.000000E+000 
 
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 
# Define analysis options 
# ----------------------- 
# Create the system of equation, a sparse solver with partial pivoting
system BandGeneral

# Create the constraint handler
constraints Plain

# Create the DOF numberer
numberer Plain

# Create the convergence test
test NormDispIncr 1.0e-2 10 1

# Create the solution algorithm
algorithm KrylovNewton

# Create the integration scheme, the DisplacementControl scheme
integrator LoadControl 0.1

# Create the analysis object
analysis Static

# initialize in case we need to do an initial stiffness iteration
initialize

#Eigen Vector outputs
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_1.out  -time  -nodeRange 1  16 -dof  1  2  3  4  5  6  eigen1 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_2.out  -time  -nodeRange 1  16 -dof  1  2  3  4  5  6  eigen2 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_3.out  -time  -nodeRange 1  16 -dof  1  2  3  4  5  6  eigen3 

set eigFID [open EigenDefaultCase_Node_EigenVector_EigenVal.out w] 
puts $eigFID [eigen  generalized  genBandArpack     3] 
close $eigFID 

analyze  1  0.0001 

remove recorders;

# Define recorder(s) 
# -------------------- 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  16 -dof  1  2  3  4  5  6  disp 

analyze     1 

remove recorders;

# Clean up 
# -------- 