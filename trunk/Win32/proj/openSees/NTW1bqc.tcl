# NTW1q1.tcl all shell element with quad3d element, no colomn elements, so there 3DOF/node
# the appended dofs there is fix in x-, y- axials
# use FAReinforcedConcrtePlaneStress material
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
# ---------------------------------- 
# Define material(s) 
# ---------------------------------- 
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
nDMaterial  ElasticIsotropic3D       5  +2.000000E+005  +2.500000E-001  +0.000000E+000 

# Material "TortionMat":    matTag    E    eta  
uniaxialMaterial  Elastic       6  +2.000000E+009  +0.000000E+000 

# ---------------------------------- 
# Define section(s) 
# ---------------------------------- 

# Section "WebTip":    secTag 
section  Fiber     4  { 
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
section  Fiber     5  { 
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
set webSec 7;
section  Aggregator    $webSec  6 T  -section 4

# Section "FlageTipT":    secTag    matTag    string ...    -section secTag 
set flageSec 8
section  Aggregator    $flageSec  6 T  -section 5

# ---------------------------------- 
# Define geometry 
# ---------------------------------- 

# generate the nodes and elements
set nxFlage 2;  # divied by 2, for middle wall
set nyFlage 4; # divided by 4, for 4 story

set nxWeb 1;  
set nyWeb 4;  # divided by 4, for 4 story

set nFrame 4; # 4 story for beam and column

# define NODAL COORDINATES

node 1  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
node 2  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
node 3  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
node 4  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm]
node 5  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm]
node 6  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm]
node 7  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm]
node 8  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 7.20E+1*$in/$mm]
node 9  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm]
node 10 [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm]
node 11 [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm]
node 12 [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 1.44E+2*$in/$mm]
node 13 [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm]
node 14 [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm]
node 15 [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm]
node 16 [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.16E+2*$in/$mm]
node 17 [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
node 18 [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
node 19 [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
node 20 [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.88E+2*$in/$mm]
fix  1     1     1     1    1     1     1 
fix  2     1     1     1    1     1     1 
fix  3     1     1     1    1     1     1 
fix  4     1     1     1    1     1     1 

# ---------------------------------- 
# define ELEMENTS
# ---------------------------------- 
# Define Beam-Column Elements
set np 2;	# number of Gauss integration points for nonlinear curvature distribution
set ColumnType nonlinearBeamColumn;
set BeamType dispBeamColumn;
# ---------------------------------- 
# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran    type    tag    vec_xz 
set webTransf 1;
set flageTransf1 2;
set flageTransf2 3;

geomTransf  Linear   $webTransf  0  1  0 
geomTransf  Linear   $flageTransf1 0 1 0
geomTransf  Linear   $flageTransf2 0 -1 0

# ---------------------------------- 
# Define element(s) 
# ---------------------------------- 
element  $ColumnType      1      1      5     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      2      2      6     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      3      3      7     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      4      4      8     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      5      5      9     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      6      6     10     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      7      7     11     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      8      8     12     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType      9      9     13     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     10     10     14     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     11     11     15     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     12     12     16     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     13     13     17     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     14     14     18     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     15     15     19     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006 
element  $ColumnType     16     16     20     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006 	# columns


# ---------------------------------- 
# Define model basic for 3-Dim 3-DOF issues
# ---------------------------------- 
model  BasicBuilder  -ndm  3  -ndf  3

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------

# set fc fy E   ; #Units: L=mm, F=N, 
set wfc [expr     32.0];
set wfy [expr    430.0];
set wE  [expr    2.0e5];
set rou1 0.0039;
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
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0 $eps0
nDMaterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025
#nDMaterial FAReinforcedConcretePlateFiber 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025

#section PlateFiber secTag  ndMatTag  h
section  PlateFiber  2   15  $t

# ---------------------------------- 
# Define geometry and node and element in 3D 3dof
# ---------------------------------- 

set flageNodeStart 21;
set flageEleStart 21;

set webNodeStart [expr ($nxFlage+1)*($nyFlage+1)+21];
set webEleStart  [expr $nxFlage*$nyFlage+21];

#element  quad3d       1       1       11       15       5     $t  "PlaneStress"   15;
# block2d $nx $ny $e1 $n1 element (element arguments) {
set cmd "block2D $nxFlage $nyFlage $flageNodeStart $flageEleStart quad3d  \" $t  PlaneStress  15 \"  {
    1   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    2   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    3   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
    4   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
}"
eval $cmd;

#element  quad3d       9      16       6       10      20     $t  "PlaneStress"   15;
set cmd "block2D $nxWeb $nyWeb $webNodeStart $webEleStart quad3d  \" $t  PlaneStress  15 \"  {
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

set fileModel [open "qmodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  3 \n"
foreach nid  $nodeIDs {
  puts $fileModel "node $nid [nodeCoord $nid]"
}

puts $fileModel "element  $ColumnType      1      1      5     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      2      2      6     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      3      3      7     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      4      4      8     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      5      5      9     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      6      6     10     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      7      7     11     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      8      8     12     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType      9      9     13     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     10     10     14     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     11     11     15     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     12     12     16     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     13     13     17     $np     $flageSec     $flageTransf2  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     14     14     18     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     15     15     19     $np     $flageSec     $flageTransf1  -mass +0.000000E+000  -iter   25  +1.000000E-006" 
puts $fileModel "element  $ColumnType     16     16     20     $np     $webSec       $webTransf     -mass +0.000000E+000  -iter   25  +1.000000E-006"

# ---------------------------------- 
# Define Single Point Constraints 
# ---------------------------------- 
for { set i 0 } { $i <= $nxWeb } { incr i 1 } {
  fix [expr $webNodeStart+$i] 1 1 1;
  puts $fileModel "fix [expr $webNodeStart+$i] 1 1 1;"
}
for { set i 0 } { $i <= $nxFlage } { incr i 1 } {
  fix [expr $flageNodeStart+$i] 1 1 1;
  puts $fileModel "fix [expr $flageNodeStart+$i] 1 1 1;"
}
if { $nyWeb != $nyFlage } {
  puts "ERROR, nyWeb != $nyFlage ! "
}
# Constraints & Define nodal masses
set nePerStory [expr $nyFlage/$nFrame]
for { set i 1 } { $i <= $nyFlage} {incr i 1} {
  if { $i%$nePerStory == 0 } {
    # web and flage public node
    equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3
    puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1]	[expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3"
    # web/flage and column public nodes
    equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1]	[expr int($i/$nePerStory)*4+2]  1 2 3
    equalDOF  [expr $flageNodeStart+($nxFlage+1)*$i]	[expr int($i/$nePerStory)*4+1]  1 2 3
    equalDOF  [expr $flageNodeStart+($nxFlage+1)*($i+1)-1]	[expr int($i/$nePerStory)*4+3]  1 2 3
    equalDOF  [expr $webNodeStart+($nxWeb+1)*$i]	[expr int($i/$nePerStory)*4+4]  1 2 3
    puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1]	[expr int($i/$nePerStory)*4+2]  1 2 3"
    puts $fileModel "equalDOF  [expr $flageNodeStart+($nxFlage+1)*$i]	[expr int($i/$nePerStory)*4+1]  1 2 3"
    puts $fileModel "equalDOF  [expr $flageNodeStart+($nxFlage+1)*($i+1)-1]	[expr int($i/$nePerStory)*4+3]  1 2 3"
    puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*$i]	[expr int($i/$nePerStory)*4+4]  1 2 3"
  } else {
    equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3
    puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1]	[expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3"
  }

  for { set j 0 } { $j <= $nxFlage} {incr j 1} {
    if { $j != $nxFlage/2 } {
      fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0
      puts $fileModel "fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0"
    }
    mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1.0 1000.
    puts $fileModel "mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1.0 1000."
  }
  for { set j 0 } { $j <= $nxWeb} {incr j 1} {
    if { $j != $nxWeb } {
      fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0
      puts $fileModel "fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0"
    }
    mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1.0 1000. 1000.
    puts $fileModel "mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1.0 1000. 1000."
  }
}

set eleIDs [getEleTags]
set endEleTag [lindex $eleIDs end]
foreach nid  $eleIDs {
  puts $fileModel "element quad $nid [eleNodes $nid]"
}
close $fileModel

# ---------------------------------- 
# Define time series 
# ---------------------------------- 
# TimeSeries.tcl 

# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# ---------------------------------- 
# Start of anaysis generation 
# ---------------------------------- 
# ---------------------------------- 
# Get Initial Stiffness 
# ---------------------------------- 
initialize 
# ---------------------------------- 
# Analysis: StaticDefaultCase 
# ---------------------------------- 
# ---------------------------------- 
# Define load pattern 
# ---------------------------------- 
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $nxFlage} {incr j 1} {
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004 
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004 
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }

    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
}
# ---------------------------------- 
# Define recorder(s) 
# ---------------------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1 $endNodeTag  -dof  1  2  3  disp 
# ---------------------------------- 
# Define analysis options 
# ---------------------------------- 
# AnalysisOptn_1.tcl 

# AnalysisOptn "StaticDefault": Type: Static 

# Constraint Handler Lagrange
constraints  Plain;  
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
# ---------------------------------- 
# Reset for next analysis case 
# ---------------------------------- 
loadConst -time 0.0 
remove recorders 
wipeAnalysis 
# ---------------------------------- 
# Analysis: EigenDefaultCase 
# ---------------------------------- 
# ---------------------------------- 
# Define recorder(s) 
# ---------------------------------- 
# Recorder_2.tcl 

# Node Recorder "EigenVector":    fileName    <nodeTag>    dof    respType 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_1.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen1 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_2.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen2 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_3.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen3 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_4.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen4 
# ---------------------------------- 
# Define analysis options 
# ---------------------------------- 
# AnalysisOptn_2.tcl 
# ---------------------------------- 
# AnalysisOptn "EigenDefault": Type: Eigen 
# ---------------------------------- 
# Constraint Handler Lagrange
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
# ---------------------------------- 
# Clean up 
# ---------------------------------- 
remove recorders;
#wipe 
# ---------------------------------- 
# Reset for next analysis case 
# ---------------------------------- 
loadConst -time 0.0 
wipeAnalysis 
# ---------------------------------- 
# Analysis: PushoverCase 
# ---------------------------------- 
# ---------------------------------- 
# Define time series 
# ---------------------------------- 
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
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +1.000000E+000  0.000000E+000 
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +1.000000E+000  0.000000E+000 
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 

recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.dis  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Element -file PushoverCase_Ele_All.stress -time -eleRange 1 $endEleTag stress
recorder  Element -file PushoverCase_Ele_All.strain -time -eleRange 1 $endEleTag strain

# ---------------------------------- 
# perform the analysis
# ---------------------------------- 
set IDctrlNode [expr $webNodeStart+($nxWeb+1)*$nyFlage];
set IDctrlDOF 2;
set Dincr 0.01;
# Constraint Handler Lagrange Transformation
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

# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;  # constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator
source GeneratePeaks.tcl;
#  ---------------------------------    perform Static Cyclic Displacements Analysis
set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";  # format for screen/file output of DONE/PROBLEM analysis
set iDmax "  1 3 5 10";
set CycleType Full;      # you can do Full / Push / Half cycles with the proc
set Ncycles 1;      # specify the number of cycles at each peak
set Fact 1.0;

foreach Dmax $iDmax {
  set iDstep [GeneratePeaks $Dmax $Dincr $CycleType $Fact];  # this proc is defined above
  for {set i 1} {$i <= $Ncycles} {incr i 1} {
    set zeroD 0
    set D0 0.0
    foreach Dstep $iDstep {
      set D1 $Dstep
      set Dincr [expr $D1 - $D0]
      integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
      analysis Static
      # ----------------------------------------------first analyze command------------------------
      set ok [analyze 1]
      # ----------------------------------------------if convergence failure-------------------------
      if {$ok != 0} {
        # if analysis fails, we try some other stuff
        # performance is slower inside this loop  global maxNumIterStatic;      # max no. of iterations performed before "failure to converge" is ret'd
        if {$ok != 0} {
          puts "Trying Newton with Initial Tangent .."
          test NormDispIncr   $Tol 2000 0
          algorithm Newton -initial
          set ok [analyze 1]
          test $testTypeStatic $TolStatic      $maxNumIterStatic    0
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          puts "Trying Broyden .."
          algorithm Broyden 8
          set ok [analyze 1 ]
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          puts "Trying NewtonWithLineSearch .."
          algorithm NewtonLineSearch 0.8 
          set ok [analyze 1]
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
          puts $putout
          return -1
        }; # end if
      }; # end if
      # -----------------------------------------------------------------------------------------------------
      set D0 $D1;      # move to next step
    }; # end Dstep
  };    # end i
};  # end of iDmaxCycl
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
  puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
  puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}

