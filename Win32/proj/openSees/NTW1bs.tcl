# NTW1bs.tcl all shell element with MITC element, with block command
# no colomn elements, so there are 6DOF/node
# the appended dofs there is fix in x-, y- axials
# use FAReinforcedConcrtePlateFiber material
wipe;
source Units&Constants.tcl;
######################## 
# Analysis-Sequence  1 #
######################## 

# Start of model generation 
# ========================= 
set LunitTXT "mm"
set FunitTXT "N"
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  3  -ndf  6 

# Define material(s) 
# ------------------ 
# Materials.tcl 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.000000E+005  +0.000000E+000 

# Material "CoreConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       2  -3.500000E+001  -2.000000E-003  -4.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "CoverConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "reinforcedBar":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
uniaxialMaterial  Steel02       4  +3.750000E+002  +2.000000E+005  +5.000000E-002  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D       5  +2.000000E+005  +2.500000E-001  +7.600000E+003 

# Material "TortionMat":    matTag    E    eta  
uniaxialMaterial  Elastic       6  +2.000000E+009  +0.000000E+000 

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E   ; #Units: L=mm, F=N, 
set wfc [expr     32.0];
set wfy [expr    410.0];
set wE  [expr    2.0e5];
set rou1 0.0059;
set rou2 0.0026;
set t [expr $in*6.0/$mm];

# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2

# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
uniaxialMaterial ConcreteZ02  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteZ02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 

set pi [expr 2.0*asin(1.0)];
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE $eps0
nDMaterial FAReinforcedConcretePlateFiber 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
#nDMaterial RAReinforcedConcretePlateFiber 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

#section PlateFiber secTag  ndMatTag  h
section  PlateFiber  2   15  $t


# Define section(s) 
# ----------------- 
# Sections.tcl 

# Section "ElasticDefault":    secTag    E    A    Iz  Iy  G  J 
section  Elastic       1  +2.000000E+005  +1.800000E+002  +4.860000E+003  +1.500000E+003  +1.115400E+004  +3.916000E+003 


# Section "PlateFiber":    secTag    matTag    h 
section  PlateFiber       3       5  +1.000000E-002 

# Section "WebTip":    secTag 
section  Fiber       4  { 
    # PatchBox "CoverCon":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3    10     1  [expr -3.000*$in/$mm]  [expr -10.875*$in/$mm]  [expr +3.000*$in/$mm]  [expr -10.875*$in/$mm]  [expr +3.000*$in/$mm]  [expr +10.500*$in/$mm]  [expr -3.000*$in/$mm]  [expr +10.500*$in/$mm]
    patch  quad       3     1    40  [expr -3.000*$in/$mm]  [expr -10.500*$in/$mm]  [expr -2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr -2.625*$in/$mm]  [expr +10.500*$in/$mm]  [expr -3.000*$in/$mm]  [expr +10.500*$in/$mm]
    patch  quad       3     1    40  [expr +2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr +3.000*$in/$mm]  [expr -10.500*$in/$mm]  [expr +3.000*$in/$mm]  [expr +10.500*$in/$mm]  [expr +2.625*$in/$mm]  [expr +10.500*$in/$mm]
    #patch  quad      3    22     8  [expr +3.000*$in/$mm]  [expr -10.875*$in/$mm]  [expr -3.000*$in/$mm]  [expr -10.125*$in/$mm]  [expr -3.000*$in/$mm]  [expr -10.125*$in/$mm]  [expr +3.000*$in/$mm]  [expr -10.875*$in/$mm]
    # PatchQuad "CoreCon":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2     4    14  [expr -2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr +2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr +2.625*$in/$mm]  [expr +10.500*$in/$mm]  [expr -2.625*$in/$mm]  [expr +10.500*$in/$mm] 

    fiber  [expr -2.625*$in/$mm]  [expr +10.500*$in/$mm]  [expr 1.1*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr +10.500*$in/$mm]  [expr 1.1*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr +7.000*$in/$mm ]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr +7.000*$in/$mm ]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr +3.500*$in/$mm ]  [expr 3.1*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr +3.500*$in/$mm ]  [expr 3.1*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr +0.000*$in/$mm ]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr +0.000*$in/$mm ]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr -3.500*$in/$mm ]  [expr 4.4*$in2/$mm2]       4
    fiber  [expr +2.625*$in/$mm]  [expr -3.500*$in/$mm ]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr -7.000*$in/$mm ]  [expr 3.1*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr -7.000*$in/$mm ]  [expr 3.1*$in2/$mm2]       4 
    fiber  [expr -2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr 4.4*$in2/$mm2]       4 
    fiber  [expr +2.625*$in/$mm]  [expr -10.500*$in/$mm]  [expr 4.4*$in2/$mm2]       4 
} 

# Section "FlageTip":    secTag 
section  Fiber       5  { 
    # PatchBox "CoverCon":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3    30     1  [expr -9.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr +9.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr +9.000*$in/$mm]  [expr +3.000*$in/$mm]  [expr -9.000*$in/$mm]  [expr +3.000*$in/$mm]
    patch  quad       3     1     8  [expr -9.000*$in/$mm]  [expr -2.625*$in/$mm]  [expr -8.625*$in/$mm]  [expr -2.625*$in/$mm]  [expr -8.625*$in/$mm]  [expr +2.625*$in/$mm]  [expr -9.000*$in/$mm]  [expr +2.625*$in/$mm]
    patch  quad       3    30     1  [expr -9.000*$in/$mm]  [expr -3.000*$in/$mm]  [expr +9.000*$in/$mm]  [expr -3.000*$in/$mm]  [expr +9.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr -9.000*$in/$mm]  [expr +2.625*$in/$mm]
    # PatchQuad "CoreCon":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2    28     8  [expr -8.625*$in/$mm]  [expr -2.625*$in/$mm]  [expr +9.000*$in/$mm]  [expr -2.625*$in/$mm]  [expr +9.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr -8.625*$in/$mm]  [expr +2.625*$in/$mm] 
        # PatchQuad "CoreCon":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    #Fiber "FlageTip":    y    z    A    matTag 
    fiber  [expr -8.625*$in/$mm]  [expr +2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr +9.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr +5.500*$in/$mm]  [expr +2.625*$in/$mm]  [expr +3.1*$in2/$mm2]       4 
    fiber  [expr +2.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr -1.500*$in/$mm]  [expr +2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr -5.000*$in/$mm]  [expr +2.625*$in/$mm]  [expr +3.1*$in2/$mm2]       4 
    fiber  [expr -8.625*$in/$mm]  [expr -2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr +9.000*$in/$mm]  [expr -2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr +5.500*$in/$mm]  [expr -2.625*$in/$mm]  [expr +3.1*$in2/$mm2]       4 
    fiber  [expr +2.000*$in/$mm]  [expr -2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr -1.500*$in/$mm]  [expr -2.625*$in/$mm]  [expr +4.4*$in2/$mm2]       4 
    fiber  [expr -5.000*$in/$mm]  [expr -2.625*$in/$mm]  [expr +3.1*$in2/$mm2]       4 
} 

# Section "TSec":    secTag    matTag    string 
section  Uniaxial       6       1  T 


 # Section "WebTipT":    secTag    matTag    string ...    -section secTag 
section  Aggregator       7  6 T  -section 4


 # Section "FlageTipT":    secTag    matTag    string ...    -section secTag 
section  Aggregator       8  6 T  -section 5

# Define geometry 
# --------------- 
# generate the nodes and elements
set nxFlage 2;  # divied by 2, for middle wall
set nyFlage 4; # divided by 4, for 4 story
set flageNodeStart 1;
set flageEleStart 1;
# element "Shell2d":    eleTag    NodeI    NodeJ    NodeK    NodeL    secTag 
#element  shell02       1       1       11       15       5       $secTag;
# block2d $nx $ny $e1 $n1 element (element arguments) {
set cmd "block2D $nxFlage $nyFlage $flageNodeStart $flageEleStart shell02  \" 2 \"  {
    1   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    2   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    3   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
    4   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
}"
eval $cmd;

set nxWeb 2;  
set nyWeb 4;  # divided by 4, for 4 story
set webNodeStart [expr ($nxFlage+1)*($nyFlage+1)+1];
set webEleStart  [expr $nxFlage*$nyFlage+1];
#element  quad3d       9      16       6       10      20     $secTag;
set cmd "block2D $nxWeb $nyWeb $webNodeStart $webEleStart shell02  \" 2 \"  {
    1   [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm]
    2   [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    3   [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
    4   [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.88E+2*$in/$mm]
}"
eval $cmd;

set nodeIDs [getNodeTags]
set endNodeTag [lindex $nodeIDs end]
set eleIDs [getEleTags]
set endEleTag [lindex $eleIDs end]

set fileModel [open "smodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  6 \n"
foreach nid  $nodeIDs {
  puts $fileModel "node $nid [nodeCoord $nid]"
}

# Define Single Point Constraints 
# ------------------------------- 
fixZ       0.0     1     1     1     1     1     1 
puts $fileModel "fixZ       0.0     1     1     1     1     1     1 ;"

if { $nyWeb != $nyFlage } {
  puts "ERROR, nyWeb != $nyFlage ! "
}
# Constraints; Define nodal masses
for { set i 1 } { $i <= $nyFlage} {incr i 1} {
  equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3 4 5 6
  puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3 4 5 6"
  #set cNodes "";
  for { set j 0 } { $j <= $nxFlage} {incr j 1} {
    if { $j != $nxFlage/2 } {
      #fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0 0 0 0
      #puts $fileModel "fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0 0 0 0"
    }
    mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1.0 1000. 1.0 1.0 1.0
    puts $fileModel "mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1.0 1000. 1.0 1.0 1.0"
    #lappend cNodes [expr $flageNodeStart+($nxFlage+1)*$i+$j];
  }
  for { set j 0 } { $j <= $nxWeb} {incr j 1} {
    if { $j != $nxWeb } {
      #fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0 0 0 0
      #puts $fileModel "fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0 0 0 0"
    }
    mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1.0 1000. 1000. 1.0 1.0 1.0
    puts $fileModel "mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1.0 1000. 1000. 1.0 1.0 1.0"
    #lappend cNodes [expr $webNodeStart+($nxWeb+1)*$i+$j];
  }
  #set cmd "rigidDiaphragm 3 [expr $webNodeStart+($nxWeb+1)*($i+1)-1] $cNodes"
  #eval $cmd;
  #puts $fileModel "rigidDiaphragm 3 [expr $webNodeStart+($nxWeb+1)*($i+1)-1] $cNodes"
}

set eleIDs [getEleTags]
set endEleTag [lindex $eleIDs end]
foreach nid  $eleIDs {
  puts $fileModel "element shell02 $nid [eleNodes $nid]"
}
close $fileModel

# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran    type    tag    vec_xz 
set webTransf 1;
set flageTransf1 2;
set flageTransf2 3;

geomTransf  Linear   $webTransf  0  1  0 
geomTransf  Linear   $flageTransf1 0 1 0
geomTransf  Linear   $flageTransf2 0 -1 0

# nodes for bottom column elements
#node      21  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      22  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
#node      23  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      24  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
#node      25  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      26  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 7.20E+1*$in/$mm];

#fix nodes for columen elements
#fix      21     1     1     1    1     1     1 
#fix      23     1     1     1    1     1     1 
#fix      25     1     1     1    1     1     1 
#just for eigen analysis robustic
#mass      21  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000
#mass      22  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000
#mass      23  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000
#mass      24  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000
#mass      25  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000
#mass      26  +1.000000E+000  +1.000000E+000  +1.000000E+000  +1.00000E+000  +1.00000E+000  +1.00000E+000

# Define Multi Point Constraints 
# ------------------------------ 
# MPConstraint.tcl constrain the column element dofs and the shells'

# Equal DOF: MPFlage1:    mNodeTag    sNodeTag    dof 
#equalDOF       2      22  1  2  3 

# Equal DOF: MPFlage2:    mNodeTag    sNodeTag    dof 
#equalDOF      12      24  1  2  3  

# Equal DOF: MPWeb1:    mNodeTag    sNodeTag    dof 
#equalDOF      17      26  1  2  3  

# Equal DOF: MPWeb1:    mNodeTag    sNodeTag    dof 
# asumming the plane plate
#rigidDiaphragm 3    7  2   12   17;  #22 24 26
#rigidDiaphragm 3    8  3   13   18;
#rigidDiaphragm 3    9  4   14   19;
#rigidDiaphragm 3   10  5   15   20;

# ----------------- 
# Elements.tcl 
set ColumnType nonlinearBeamColumn
# Element "FlageCol":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  $ColumnType      13      21      22     2     8     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006 

# Element "FlageCol":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  $ColumnType      14      23      24     2     8     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006 

# Element "WebCol":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  $ColumnType      15      25      26     2     7     $webTransf  -mass +0.000000E+000  -iter   25  +1.000000E-006 

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

# Analysis: StaticDefaultCase   /////////////////////////////////////////////////////////////
# +++++++++++++++++++++++++++   /////////////////////////////////////////////////////////////

# Define load pattern 
# ------------------- 
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $nxFlage} {incr j 1} {
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000"
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000"
    }
    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
}

# Define recorder(s) 
# -------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1 $endNodeTag -dof  1  2  3  4  5  6 disp 

# Define analysis options   /////////////////////////////////////////////////////////////
# -----------------------   /////////////////////////////////////////////////////////////
# AnalysisOptn_1.tcl 

# AnalysisOptn "StaticDefault": Type: Static 
# ------------------------------------------ 
# Constraint Handler 
constraints  Transformation 
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
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

# Analysis: EigenDefaultCase 
# ++++++++++++++++++++++++++ 

# Define recorder(s) 
# -------------------- 
# Recorder_2.tcl 

# Node Recorder "EigenVector":    fileName    <nodeTag>    dof    respType 
#recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_1.out  -time  -nodeRange 1 $endNodeTag -dof  1  2  3  eigen1 
#recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_2.out  -time  -nodeRange 1 $endNodeTag -dof  1  2  3  eigen2 
#recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_3.out  -time  -nodeRange 1 $endNodeTag -dof  1  2  3  eigen3 
#recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_4.out  -time  -nodeRange 1 $endNodeTag -dof  1  2  3  eigen4 

# Define analysis options     /////////////////////////////////////////////////////////////
# -----------------------     /////////////////////////////////////////////////////////////
# AnalysisOptn_2.tcl          /////////////////////////////////////////////////////////////

# AnalysisOptn "EigenDefault": Type: Eigen 
# ---------------------------------------- 
# Constraint Handler 
constraints  Transformation 
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  Newton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

#set eigFID [open EigenDefaultCase_Node_EigenVector_EigenVal.out w] 
set lambda [eigen 15] 
#puts $eigFID $lambda
#close $eigFID 
record;
foreach lam $lambda {
  if { $lam > 0.0 } {
    lappend omega [expr sqrt($lam)]
    lappend f [expr sqrt($lam)/(2*$pi)]
    lappend T [expr (2*$pi)/sqrt($lam)]
    puts "Tn = [expr (2*$pi)/sqrt($lam)] s";
  } else {
    lappend omega -1.0
    lappend f "complex"
    lappend T "complex"
    puts "Tn = complex";
  }
}

# get values of eigenvectors for translational DOFs
#---------------------------------------------------
set fileEig [open "eigFile.eig" "w"]
for { set k 1 } { $k <= 15 } { incr k } {
  set tempStr "[lindex $T [expr $k-1]]"
  foreach tg $nodeIDs {
    lappend tempStr [lindex [nodeEigenvector $tg $k 1] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 2] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 3] 0]
  }
  puts $fileEig $tempStr;
}
close $fileEig 

# Clean up 
# -------- 
remove recorders;
#wipe 

# Reset for next analysis case 
# ---------------------------- 
loadConst -time 0.0 
wipeAnalysis 

# Analysis: PushoverCase 
# ++++++++++++++++++++++ 

# Define time series 
# ------------------ 
# TimeSeries.tcl 

# TimeSeries "TimeSeriesX":    tsTag    dt    filePath    cFactor 
timeSeries  Path       2  -dt  +1.000000E-002  -filePath  TimeSeriesX.thf  -factor  +1.000000E-003 

# TimeSeries "TimeSeriesY":    tsTag    dt    filePath    cFactor 
timeSeries  Path       3  -dt  +1.000000E-002  -filePath  TimeSeriesY.thf  -factor  +1.000000E-003 

# LoadPattern_3.tcl 

# LoadPattern "PushoverPattern":    patternTag    tsTag 
##pattern  Plain       3    "Linear"  { 
##    # Load    nodeTag    LoadValues 
## 
##    # SP    nodeTag    dofTag    DispValue 
##    sp 5  1 1.0
##    sp 10 1 1.0
##    sp 15 1 1.0
##    sp 25 1 1.0
##    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
## 
##    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
##} 

pattern  Plain       2    "Linear"  { 
    # Load    nodeTag    LoadValues 
     for { set j 0 } { $j <= $nxFlage} {incr j 1} {
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  +0.000000E+000  0.000000E+000  +0.000000E+000
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  0.000000E+000  0.000000E+000  0.000000E+000
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    # SP    nodeTag    dofTag    DispValue 
    #sp 5  1 1.
    #sp 10 1 1.
    #sp 15 1 1.
    #sp 25 1 1.
    #sp 5  2 1.
    #sp 10 2 1.
    #sp 15 2 1.
    #sp 25 2 1.
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 

recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.dis  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Element -file PushoverCase_Ele_All.stress -time -eleRange 1 $endEleTag stress
recorder  Element -file PushoverCase_Ele_All.strain -time -eleRange 1 $endEleTag strain

# perform the analysis  /////////////////////////////////////////////////////////////
set IDctrlNode [expr $webNodeStart+($nxWeb+1)*$nyFlage];
set IDctrlDOF 2;
set Dincr 0.01;
# Constraint Handler Transformation
constraints  Plain
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

# set up analysis parameters   /////////////////////////////////////////////////////////////
#source LibAnalysisStaticParameters.tcl;  
#source GeneratePeaks.tcl;
#  perform Static Cyclic Displacements Analysis   /////////////////////////////////////////////////////////////

set iDmax1 [expr 0.08*$in/$mm]
set iDmax2 [expr 0.12*$in/$mm]
set iDmax3 [expr 0.3*$in/$mm]
set iDmax4 [expr 0.4*$in/$mm]
set iDmax5 [expr 0.6*$in/$mm]
set iDmax6 [expr 1.1*$in/$mm]

# perform the analysis    /////////////////////////////////////////////////////////////
#                0.08   -0.08  0.12   -0.12   0.3   -0.3   0.4   -0.4   0.6   -0.6   1.1   -1.1
set numSteps     13970; #{   508   1016  1270    1524  2667   3810  4445   5080  6350   7620 10795  13970 }
set numIters       500; #{   100    200   200     200   200    100   200    200   200    200   500    500 }
set increments   0.004; #{ 0.004 -0.004 0.004  -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 }
 
for {set i 0} {$i<1} {incr i 1} {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set increment [lindex $increments $i]
 
    integrator DisplacementControl 13 1 $increment
    if {$numIter == 0} {
	set numIter 1
    } 
    test NormDispIncr 1e-3 $numIter 5
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
        set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
        puts $putout
        test NormDispIncr 1e-3 $numIter 5
        set ok [analyze 1]
      }; # end if
      incr j 1;
      record;
      puts $end
    }
}
