# Create ModelBuilder with 2 dimensions and 2 DOF/node
# 修改自 SW4 for nonlinearBS material
# 改为一个单元水平荷载
# 英制
# 为了对比弹性多轴本构对quad单元的刚度，增加了弹性多轴本构

wipe;
model basic -ndm 2 -ndf 2

source units&constants.tcl;  # units mm, N, s

set PI [expr 2*asin(1.0)]; # define constants
#set g [expr 32.2*$ft/pow($sec,2)/$mm]; # gravitational acceleration

# define input displacement history
set DHfile1 "TimeSeriesY.thf";

## create the material
## MATERIAL parameters
# nominal concrete compressive strength

set CompressiveStress [expr -5.35*$ksi*$mm2];
set fc [expr $CompressiveStress];
set Ec [expr 2.0*$fc/-0.003]; #57*sqrt(-$fc)];

# unconfmed concrete
set fclU $fc; # UNCONFINED concrete (todeschini parabolic model), maximum stress
set epslU -0.003; # strain at maximum strength of unconfmed concrete
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
set Fy [expr 60*$ksi/$MPa]; # STEEL yield stress
set Es [expr 29000*$ksi/$MPa]; # modulus of steel
set epsY [expr $Fy/$Es]; # steel yield strain
set epsSH [expr 8*$epsY];
set Esh [expr 0.02*$Es]; # tangent stiffness at onset of StrainHardening
set Bs 0.01; # strain-hardening ratio
set Fy1 [expr 1.01*$Fy]; # steel stress post-yield
set epsY1 [expr $epsY+($Fy1-$Fy)/($Bs*$Es)]; # steel strain post-yield
set Fu [expr 1.16*$Fy];
set epsU 1.0; # ultimate strain of steel#

#MATERIAL assignments
set rou1 0;
set rou2 0;
set rouv 0.0246;
set rouh1 0.00246; # #10@70
# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $Fy $Es $fc $rouv
uniaxialMaterial SteelZ01 12 $Fy $Es $fc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f’c ec0
uniaxialMaterial ConcreteL02 14 [expr -$fc] -0.003
uniaxialMaterial ConcreteL02 15 [expr -$fc] -0.003

#nDMaterial CSMMRCPlaneStress $IDMat 0.0 11 12 14 15 [expr 0.5*$PI] [expr 0.0*$PI] $rouv $rouh1 $fc $Fy $Es 0.002
                                                                                                              #aggr dx  dy
nDMaterial MCFTRCPlaneStress $IDMat 0.0 11 12 14 15 [expr 0.5*$PI] [expr 0.0*$PI] $rouv $rouh1 $fc $Fy $Es 0.002 10 140 140.


# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D 6 2e4 0.2  0 
# Section "ElaMembranePlateSec":    secTag    E    v    h    rho 
#section  ElasticMembranePlateSection       4  +2.000000E+005  +2.500000E-001  +1.000000E-001  0.0

# Define parameters
set Quad "quad"

if { $Quad == "quad" } {
  set eleArgs "PlaneStress2D $IDMat";
}

# Node    tag    xCrd    yCrd 
node       1       0       0
node       2 [expr 200] 0
node       3       0 [expr 200]
node       4 [expr 200] [expr 200] 

# SPC    tag    Dx    Dy    Dz    Rx    Ry    Rz 
fix       1     1     1       
fix       2     0     1       
fix       3     0     0       
fix       4     0     0      

# Node    tag    mx    my    mz    mIx    mIy    mIz 
mass       1 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       2 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       3 [expr 1000.*$lbf] [expr 1000.*$lbf]
mass       4 [expr 1000.*$lbf] [expr 1000.*$lbf]

# ELEMENT DEFINITION
element quad     1       1       2       4       3     [expr 6.*$in/$mm]  "PlaneStress2D"  $IDMat
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
  prp -1750 -5800 6000
  #vrp 0 -500 250
  vup 0 1 0
  #vpn -1 -1 0.5
  viewWindow -400 400 -400 400
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
pattern Plain 2 1 {
	load 3 0.0 1.e3;
	load 4 0.0 1.e3;
}
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
system  UmfPack 
# Analysis Type 
analysis  Static 

analyze 1

loadConst;

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
system  UmfPack 
# Analysis Type 
analysis  Static 

# perform the analysis
# analyze 17380
set numSteps   {  100    400   300 }; #   400   600    800 1200  1600  2050  2500 2960  3470 1000    10 };    #}
set numIters   {  200    200   200 }; #   200   200    100  100   100   100   100  500   500  500   300 };    #}
set increments {  0.1   0.05   0.1 }; #  -0.1   0.1   -0.1  0.1  -0.1   0.1  -0.1  0.1  -0.1  0.1  -0.1 };     #}

for { set i 0 } { $i<3 } { incr i 1 } {
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
    #print ele
}
#remove recorders;
