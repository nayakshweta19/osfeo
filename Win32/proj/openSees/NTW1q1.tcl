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
model  BasicBuilder  -ndm  3  -ndf  3

# Define geometry 
# --------------- 
# NodeCoord.tcl 

# Node    tag    xCrd    yCrd    zCrd 
node       1  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
node       2  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
node       3  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm];
node       4  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm];
node       5  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm];
node       6  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
node       7  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
node       8  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm];
node       9  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm];
node      10  [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm];
node      11  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
node      12  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
node      13  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 1.44E+2*$in/$mm];
node      14  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.16E+2*$in/$mm];
node      15  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm];
node      16  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm];  #21
node      17  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 7.20E+1*$in/$mm];  #22
node      18  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 1.44E+2*$in/$mm];  #23
node      19  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.16E+2*$in/$mm];  #24
node      20  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.88E+2*$in/$mm];  #25

# nodes for bottom column elements
#node      21  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      22  [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
#node      23  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      24  [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 7.20E+1*$in/$mm];
#node      25  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm];
#node      26  [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 7.20E+1*$in/$mm];

# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Dz 
fix       1     1     1     1  ; #  1     1     1 
fix       6     1     1     1  ; #  1     1     1 
fix      11     1     1     1  ; #  1     1     1 
fix      16     1     1     1  ; #  1     1     1 

fix       2     0     1     0  ;
fix       3     0     1     0  ;
fix       4     0     1     0  ;
fix       5     0     1     0  ;
fix      12     0     1     0  ;
fix      13     0     1     0  ;
fix      14     0     1     0  ;
fix      15     0     1     0  ;
fix      17     1     0     0  ;
fix      18     1     0     0  ;
fix      19     1     0     0  ;
fix      20     1     0     0  ;
#fix nodes for columen elements
#fix      21     1     1     1    1     1     1 
#fix      23     1     1     1    1     1     1 
#fix      25     1     1     1    1     1     1 

# Define nodal masses
# ------------------- 
# NodeMass.tcl 

# Node    tag    mx    my    mz 
mass       1  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000  
mass       2  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000  
mass       3  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000  
mass       4  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass       5  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass       6  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass       7  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass       8  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass       9  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      10  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      11  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      12  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      13  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      14  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      15  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      16  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      17  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      18  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      19  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000
mass      20  +5.800000E+003  +5.800000E+003  +5.800000E+003 ; # +1.00000E+000  +1.00000E+000  +1.00000E+000

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
set wfy [expr    430.0];
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
uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 

set pi [expr 2.0*asin(1.0)];
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0 $eps0
#nDMaterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025
nDMaterial FAReinforcedConcretePlateFiber 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025

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


# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran    type    tag    vec_xz 
set webTransf 1;
set flageTransf1 2;
set flageTransf2 3;

#geomTransf  Linear   $webTransf  0  1  0 
#geomTransf  Linear   $flageTransf1 0 1 0
#geomTransf  Linear   $flageTransf2 0 -1 0

# Define element(s) 
# ----------------- 
# Elements.tcl 
set ColumnType nonlinearBeamColumn
# element quad eleID node1 node2 node3 node4 thick  type         matID     pressure    rho    b1    b2 
element  quad3d       1       1       6       7       2     $t  "PlaneStress"   15;
element  quad3d       2       2       7       8       3     $t  "PlaneStress"   15;
element  quad3d       3       3       8       9       4     $t  "PlaneStress"   15;
element  quad3d       4       4       9      10       5     $t  "PlaneStress"   15;
element  quad3d       5       6      11      12       7     $t  "PlaneStress"   15;
element  quad3d       6       7      12      13       8     $t  "PlaneStress"   15;
element  quad3d       7       8      13      14       9     $t  "PlaneStress"   15;
element  quad3d       8       9      14      15      10     $t  "PlaneStress"   15;
element  quad3d       9      16       6       7      17     $t  "PlaneStress"   15;
element  quad3d      10      17       7       8      18     $t  "PlaneStress"   15;
element  quad3d      11      18       8       9      19     $t  "PlaneStress"   15;
element  quad3d      12      19       9      10      20     $t  "PlaneStress"   15;
# element "Shell2d":    eleTag    NodeI    NodeJ    NodeK    NodeL    secTag 
#element  shell02       1       1       6       7       2       2
#element  shell02       2       2       7       8       3       2
#element  shell02       3       3       8       9       4       2
#element  shell02       4       4       9      10       5       2
#element  shell02       5       6      11      12       7       2
#element  shell02       6       7      12      13       8       2
#element  shell02       7       8      13      14       9       2
#element  shell02       8       9      14      15      10       2
#element  shell02       9      16       6       7      17       2
#element  shell02      10      17       7       8      18       2
#element  shell02      11      18       8       9      19       2
#element  shell02      12      19       9      10      20       2

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

# Analysis: StaticDefaultCase 
# +++++++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
# LoadPattern_1.tcl 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    load       5  +0.000000E+000  +0.000000E+000  -4.000000E+004   +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load      10  +0.000000E+000  +0.000000E+000  -6.000000E+004   +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load      15  +0.000000E+000  +0.000000E+000  -4.000000E+004   +0.000000E+000  +0.000000E+000  +0.000000E+000 
    load      20  +0.000000E+000  +0.000000E+000  -4.000000E+004   +0.000000E+000  +0.000000E+000  +0.000000E+000 

    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
}

# Define recorder(s) 
# -------------------- 
# Recorder_1.tcl 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  20 -dof  1  2  3  4  5  6 disp 

# Define analysis options 
# ----------------------- 
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
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_1.out  -time  -nodeRange 1  20 -dof  1  2  3  eigen1 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_2.out  -time  -nodeRange 1  20 -dof  1  2  3  eigen2 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_3.out  -time  -nodeRange 1  20 -dof  1  2  3  eigen3 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_4.out  -time  -nodeRange 1  20 -dof  1  2  3  eigen4 


# Define analysis options 
# ----------------------- 
# AnalysisOptn_2.tcl 

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

set eigFID [open EigenDefaultCase_Node_EigenVector_EigenVal.out w] 
puts $eigFID [eigen generalized 15] 
close $eigFID 
record;

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
    load  5  0.0  1.0  0   0 0 0       
    load 10  0.0  1.0  0   0 0 0 
    load 15  0.0  1.0  0   0 0 0 
    load 20  0.0  1.0  0   0 0 0 
    #load  5 0   1.0  0   0 0 0
    #load 10 0   1.0  0   0 0 0
    #load 15 0   1.0  0   0 0 0 
    #load 25 0   1.0  0   0 0 0 
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

recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  20 -dof  1  2  3  disp 

# ---------------------------------------- 
# Constraint Handler 
constraints  Penalty 10e12 10e12;  
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl 10 1 0.01
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

# perform the analysis
set IDctrlNode 10;
set IDctrlDOF 2;
set Dincr 0.01;
# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;  # constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator
source GeneratePeaks.tcl;
#  ---------------------------------    perform Static Cyclic Displacements Analysis
set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";  # format for screen/file output of DONE/PROBLEM analysis
set iDmax "0.5  1 3 5 10";
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

