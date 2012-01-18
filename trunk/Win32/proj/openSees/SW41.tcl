# Create ModelBuilder with 2 dimensions and 2 DOF/node
# 修改自 SW4 for nonlinearBS material
# 改为一个单元水平荷载
# 英制
# 为了对比弹性多轴本构对quad单元的刚度，增加了弹性多轴本构

wipe;
model basic -ndm 2 -ndf 2

# define UNITS
set in 25.4; # define basic units — output units  mm
set kip 4453.74; # define basic units — output units  N
set sec 1.; # define basic units - output units
set ft  [expr 12.*$in]; # define engineering units
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set lbf [expr $psi*$in*$in]; # pounds force
set pcf [expr $lbf/pow($ft,3)]; # pounds per cubic foot
set in2 [expr $in*$in]; # inchA2
set in4 [expr $in*$in*$in*$in]; # inchA4
set cm [expr $in/2.54]; # centimeter, needed for displacement input in
set PI [expr 2*asin(1.0)]; # define constants
set g [expr 32.2*$ft/pow($sec,2)]; # gravitational acceleration

# define input displacement history
set DHfile1 "TimeSeriesY.thf";

## create the material
## MATERIAL parameters
# nominal concrete compressive strength

set CompressiveStress [expr -5.35*$ksi];
set fc [expr $CompressiveStress];
set Ec [expr 57000.*$psi*sqrt(abs($fc)/$psi)];

# unconfmed concrete
set fclU $fc; # UNCONFINED concrete (todeschini parabolic model), maximum stress
set epslU -0.0022; # strain at maximum strength of unconfmed concrete
set fc2U [expr 0.2*$fclU]; # ultimate stress
set eps2U -0.01; # strain at ultimate stress
set lambda 0.1; # ratio between unloading slope at $eps2 and initial slope $Ec

# tensile-strength properties
set ftU [expr -0.14*$fclU]; # tensile strength +tension
set Ets [expr $ftU/0.002]; # tension softening stiffness

# DEFINE MATERIAL tags
set IDMat 1;
set IDMatFlg 2;
set IDMatSlab 3;
# back-bone steel stress-strain curve parameters for all materials
set Fy [expr 60*$ksi]; # STEEL yield stress
set Es [expr 29000*$ksi]; # modulus of steel
set epsY [expr $Fy/$Es]; # steel yield strain
set epsSH [expr 8*$epsY];
set Esh [expr 0.02*$Es]; # tangent stiffness at onset of StrainHardening
set Bs 0.01; # strain-hardening ratio
set Fy1 [expr 1.01*$Fy]; # steel stress post-yield
set epsY1 [expr $epsY+($Fy1-$Fy)/($Bs*$Es)]; # steel strain post-yield
set Fu [expr 1.16*$Fy];
set epsU 1.0; # ultimate strain of steel#

# hysteretic steel material parameters -- baseline
set pinchX 0.6; # pinching parameter for hysteretic model
set pinchY 0.7; # pinching parameter for hysteretic model

set damage1 0.1; # damage parameter for hysteretic model
set damage2 0.2; # damage parameter for hysteretic model
set betaMUsteel 0.4; # degraded unloading stiffness for hysteretic material based on MUA(-beta)

# steel02 and steel03 parameters — baseline
set R0 18;  # control the transition from elastic to plastic branches,
set cR1 0.925; # control the transition from elastic to plastic branches,
set cR2 0.15; # control the transition from elastic to plastic branches
set a2 0.1; # isotropic hardening parameter, associated with al
set a1 [expr $a2*($Fy/$Es)]; # isotropic hardening parameter, increase of comp.
                            # yield envelope as proportion of yield strength after a plastic strain
set a4 0.1; #ft isotropic hardening parameter, associated with a3
set a3 [expr $a4*($Fy/$Es)]; # isotropic hardening parameter, increase of tension yield envelope
                             # as proportion of yield strength after a plastic strain

# Reinforcing direction, D, and ratio, R, and diameter, d.
set rD1 [expr 0.0*$PI];
set rR1 0.0339;
set rdmm1 [expr 0.2362*$in]; #mm

set rD2 [expr 0.5*$PI];
set rR2 0.0250;
set rdmm2 [expr 0.2362*$in]; #mm

set FrR1 0.0123;
set FrR2 0.0686;
set Frd1 0.2362;
set Frd2 0.4724;

set SR1 0.003;
set SR2 0.06;

#MATERIAL assignments

nDMaterial NonlinearBS $IDMat $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $rR1 $rdmm1 $rD2 $rR2 $rdmm2 \
						$Fy $epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

nDMaterial NonlinearBS $IDMatFlg $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $FrR1 $Frd1 $rD2 $FrR2 $Frd2 \
						$Fy $epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

nDMaterial NonlinearBS $IDMatSlab $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $SR1 $rdmm1 $rD2 $SR2 $rdmm2 \
						$Fy $epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D 6 30000 0.2  0 
# Section "ElaMembranePlateSec":    secTag    E    v    h    rho 
#section  ElasticMembranePlateSection       4  +2.000000E+005  +2.500000E-001  +1.000000E-001  0.0

# Define parameters
set Quad "quad"

if { $Quad == "quad" } {
  set eleArgs "PlaneStress2D $IDMat";
}

# Node    tag    xCrd    yCrd 
node       1       0       0
node       2 [expr 1000.*$in] 0
node       3       0 [expr 1000.*$in]
node       4 [expr 1000.*$in] [expr 1000.*$in] 

# SPC    tag    Dx    Dy    Dz    Rx    Ry    Rz 
fix       1     1     1       
fix       2     1     1       
fix       3     0     0       
fix       4     0     0      

# Node    tag    mx    my    mz    mIx    mIy    mIz 
mass       1 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       2 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       3 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       4 [expr 1000.*$lbf] [expr 1000.*$lbf]

# ELEMENT DEFINITION
element quad     1       1       2       4       3     [expr 12.*$in]  "PlaneStress2D"   $IDMat
# fail element bbarQuad 1       1       2       4       3  $IDMat
# Element "quad":    eleTag    NodeI    NodeJ    NodeK    NodeL    h    type    matTag    pressure    rho    b1    b2 
#element  quad    1        1       2       4       3    [expr 6.*$in]  "PlaneStress"    6  +0.000000E+000  +0.000000E+000  +0.000000E+000 +0.000000E+000 
# Start of anaysis generation 
# =========================== 

# Get Initial Stiffness 
# --------------------- 
initialize 
puts "Model Built"

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -1250 -1580 1600
  vrp 0 -1500 2500
  vup 0 1 0
  vpn 0.4 0.6 0.8
  viewWindow -3000 6000 -3000 6000
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 1

#
# Apply the displacement history and Perform a transient analysis
#
# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    load       3  1.000000E+000   0.000000E+000 
    load       4  1.000000E+000   0.000000E+000 
 
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 

# Define recorder(s) 
# -------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 3  4 -dof  1  2  3  disp 

# Define analysis options 
# ----------------------- 
# AnalysisOptn_1.tcl 

# AnalysisOptn "StaticDefault": Type: Static 
# ------------------------------------------ 
# Constraint Handler 
constraints  Plain 
# Convergence Test 
test  NormDispIncr  +1.000000E-003   50     0     2 
# Integrator 
integrator  DisplacementControl 4 1 +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  SparseGEN; #UmfPack 
# Analysis Type 
analysis  Static 

# perform the analysis
# analyze 17380
set numSteps   {  1000     4000   3000    4000   6000    8000 12000  16000  20500  25000 29600  34700 10000    10 };    #}
set numIters   {  200     200   200    200   200    100  100   100   100   100  500   500  500   300 };    #}
set increments {  0.1 -0.05  .1   -0.1  0.1  -0.1 0.1 -0.1  0.1 -0.1 0.1 -0.1 0.1 -0.1 };     #}

for { set i 0 } { $i<14 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set increment [lindex $increments $i]
 
    integrator DisplacementControl 4 1 $increment
    if {$numIter == 0} {
	set numIter 1
    } 
    test NormDispIncr 1e-3 $numIter 5
    analyze $numStep
    puts $i
    print ele
}
remove recorders;
