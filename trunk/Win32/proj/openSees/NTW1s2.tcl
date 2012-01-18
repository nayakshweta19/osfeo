# ----------------------------
# Start of model generation
# ----------------------------
wipe;
source Units&Constants.tcl;
# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 3 -ndf 6

# create the material
nDMaterial ElasticIsotropic   1   4.0e4   0.2  6.75e-6;

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 34.5;  #[expr (6880.+6600.+7700.+6850.)/4.0*$psi/$MPa];
set wfy 414.0; #[expr (76.+63.+60)/3.0*$ksi/$MPa];
set wE  1.95e5; #[expr (29000+28000+28700)/3.0*$ksi/$MPa];
set rou1 0.0259;
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
#uniaxialMaterial Concrete09 13 -20.1  -0.0021  -5.6  -0.03  0.14  +2.6  +1300 
#uniaxialMaterial Concrete09 14 -20.1  -0.0021  -5.6  -0.03  0.14  +2.6  +1300 

set pi 3.141592654;
# NDMaterial: ReinforceConcretePlateFiber
#                                        tag    rho   s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

#nDMaterial RAReinforcedConcretePlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
nDMaterial CSMMRCPlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

# ---------------------------------------
# Define Section tag
# ---------------------------------------
#section PlateFiber secTag  ndMatTag  h
set elastEleSecTag 2;
set thick [expr $in*6.0/$mm];
section  PlateFiber  $elastEleSecTag   1  $thick

set eleSecTag 3;
section  PlateFiber  $eleSecTag   2  $thick  

# ---------------------------------------
# Define geometry
# ---------------------------------------

# define some  parameters

set mass 0.0;

set xWidth [expr 7.8E+1*$in/$mm];
set yWidth [expr 5.4E+1*$in/$mm];
set height [expr 2.88E+2*$in/$mm];

set numX 2; # num elements in x direction
set numY 2; # num elements in y direction ,divied by 2, for middle wall
set numZ 4; # num elements in z direction ,divied by 4, for 4 story\
#     ^
#     | Y
#       
#     |
#     |
#     |
#     |--------------  -> X  
#     |
#     |
#     |
#
set nodeNum 1

set incrX [expr $xWidth/(1.0*$numX)]
set incrY [expr $yWidth/(1.0*$numY)]
set incrZ [expr $height/(1.0*$numZ)]

set fileModel [open "smodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  6 \n"

#web wall
set yLoc [expr 2.7E+1*$in/$mm];
set zLoc 0.0;
for {set i 0} {$i <= $numZ} {incr i 1} {
    set xLoc 0.0;
    for {set j 0} {$j <= $numX} {incr j 1} {
	    node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass $mass $mass $mass
	    puts $fileModel "node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass $mass $mass $mass"
	    incr nodeNum
	    set xLoc [expr $xLoc + $incrX]
    }
    set zLoc [expr $zLoc + $incrZ]
}

set eleTag 1
for {set j 0} {$j <$numZ} {incr j 1} {
    set iNode [expr 1 + $j*($numX+1)]
    set jNode [expr $iNode + 1]
    set kNode [expr $iNode + $numX + 2]
    set lNode [expr $kNode - 1]
    if {  $j == [expr $numZ-1] } { ;  # $j == 0 ||
        for {set i 0} {$i <$numX} {incr i 1} {
	        element shell02 $eleTag $iNode $jNode $kNode $lNode $elastEleSecTag 
	        puts $fileModel "element shell $eleTag $iNode $jNode $kNode $lNode $elastEleSecTag"
	        incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
        }
    } else {
    	for {set i 0} {$i <$numX} {incr i 1} {
	        element shell02 $eleTag $iNode $jNode $kNode $lNode $eleSecTag 
	        puts $fileModel "element shell $eleTag $iNode $jNode $kNode $lNode $eleSecTag"
	        incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
        }
	}
}
for {set i 1} {$i <= [expr $numX+1] } {incr i 1} {
    fix $i 1 1 1 1 1 1
    puts $fileModel "fix $i 1 1 1 1 1 1"
}

#for {set i [expr $numX+2]} {$i <= [expr $nodeNum-1] } {incr i 1} {
#    fix $i 0 1 0 0 0 0 
#    puts $fileModel "fix $i 0 1 0 0 0 0"
#}
# end web wall

# flage wall
set xLoc 0.0;
set zLoc 0.0;

set nodeYdirnStart $nodeNum

for {set i 0} {$i <= $numZ} {incr i 1} {
    set yLoc 0.0;
    for {set j 0} {$j <= $numY} {incr j 1} {
	    node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass $mass $mass $mass
	    puts $fileModel "node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass $mass $mass $mass"
	    incr nodeNum
	    set yLoc [expr $yLoc + $incrY]
    }
    set zLoc [expr $zLoc + $incrZ]
}

for {set j 0} {$j <$numZ} {incr j 1} {
    set iNode [expr $nodeYdirnStart + $j*($numY+1)]
    set jNode [expr $iNode + 1]
    set kNode [expr $iNode + $numY + 2]
    set lNode [expr $kNode - 1]
    if { $j == [expr $numZ-1] } { ; # $j == 0 ||
        for {set i 0} {$i <$numY} {incr i 1} {
	        element shell02 $eleTag $iNode $jNode $kNode $lNode $elastEleSecTag
	        puts $fileModel "element shell $eleTag $iNode $jNode $kNode $lNode $elastEleSecTag"
	        incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
        }
    } else {
    	for {set i 0} {$i <$numY} {incr i 1} {
	        element shell02 $eleTag $iNode $jNode $kNode $lNode $eleSecTag
	        puts $fileModel "element shell $eleTag $iNode $jNode $kNode $lNode $eleSecTag"
	        incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
        }
    }
}

for {set i 1; set j $nodeYdirnStart} {$i <= [expr $numY+1]} {incr i 1; incr j 1} {
    fix $j 1 1 1 1 1 1 
    puts $fileModel "fix $j 1 1 1 1 1 1"
}

#for {set i [expr $nodeYdirnStart+$numY+1]} {$i <= [expr $nodeNum-1] } {incr i 1} {
#    fix $i 1 0 0 0 0 0 
#    puts $fileModel "fix $i 1 0 0 0 0 0"
#}
# end flage wall

# constraint public nodes
for {set i [expr 1+$numX+1]; set j [expr $nodeYdirnStart+1+$numY+$numY/2]} {$i < $nodeYdirnStart} {incr i [expr $numX+1]; incr j [expr $numY+1]} {
    #remove sp $i 2
    #remove sp $j 1
    equalDOF $i $j 1 2 3
   # puts $fileModel "remove sp $i 2"
    #puts $fileModel "remove sp $j 1"
    puts $fileModel "equalDOF $i $j 1 2 3"
}

#---------------------------------------------------
# Perform the eigenVector output
#---------------------------------------------------
set lambda [eigen 15];
record; puts "no graivity analysis period and frequency"
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

timeSeries Linear 1
pattern  Plain       1    "Linear"  { 
puts $fileModel "pattern  Plain       1 {"
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $numY} {incr j 1} {
      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 
      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 "
    }
    for { set j 1 } { $j <= $numX} {incr j 1} {
      load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 
      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 "
    }
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
puts $fileModel "}"
}

# Constraint Handler Lagrange
constraints  Plain;  
# Convergence Test 
test  NormDispIncr  +1.000000E-003    250     0   
#test EnergyIncr  1.0e-12    10         0
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations BandGeneral SparseGeneral UmfPackProfileSPD
system  BandGeneral 
# Analysis Type 
analysis  Static 
#---------------------------------------------------
# Perform the analysis
#---------------------------------------------------
set ok [analyze   10]
puts "Gravity Static Analysis OK= $ok (0 for ok, other for error)";

loadConst -time 0.0
#---------------------------------------------------
# Perform the eigenVector output
#---------------------------------------------------
set lambda [eigen 15];
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
set nodeIDs [getNodeTags];
set endNodeTag [lindex $nodeIDs end];
set eleIDs [getEleTags];
set endEleTag [lindex $eleIDs end];
# get values of eigenvectors for translational DOFs
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

# create the display
recorder display nodeNum 10 10 800 600 -wipe
# next three commmands define viewing system, all values in global coords
prp -80000 -10000 40000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 0 1;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -4000 4000 -4000 4000; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 -1 10; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# create the display
recorder display Deform 810 10 800 600 -wipe
# next three commmands define viewing system, all values in global coords
prp -10000 -100000 70000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 0 1;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -4000 4000 -4000 4000; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 5 1; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
loadConst -time 0.0

# ---------------------------------- 
# Analysis: PushoverCase 
# ---------------------------------- 
# ---------------------------------- 
# Define time series 
# ---------------------------------- 

# TimeSeries "TimeSeriesX":    tsTag    dt    filePath    cFactor 
timeSeries  Path       2  -dt  +1.000000E-002  -filePath  TimeSeriesX.thf  -factor  +1.000000E-003 

# TimeSeries "TimeSeriesY":    tsTag    dt    filePath    cFactor 
timeSeries  Path       3  -dt  +1.000000E-002  -filePath  TimeSeriesY.thf  -factor  +1.000000E-003 

# LoadPattern "PushoverPattern":    patternTag    tsTag 

pattern  Plain       3    "Linear"  { 
puts $fileModel "pattern  Plain       3  Linear {"
	# Load    nodeTag    LoadValues 

    # SP    nodeTag    dofTag    DispValue 
    for { set j 0 } { $j <= $numY} {incr j 1} {
      #load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      #puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #load on 1 direction
      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #sp [expr $nodeYdirnStart+($numY+1)*$numZ+$j] 1 1.0
      #puts $fileModel "sp [expr $nodeYdirnStart+($numY+1)*$numZ+$j] 1 1.0"
    }
    for { set j 1 } { $j <= $numX} {incr j 1} {
      #load [expr 1+($numX+1)*$numZ+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      #puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000 +1.000000E+000   0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #load on 1 direction
      load [expr 1+($numX+1)*$numZ+$j]  +1.000000E+000 +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +1.000000E+000  +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #sp [expr 1+($numX+1)*$numZ+$j] 1 1.0 
      #puts $fileModel "sp [expr 1+($numX+1)*$numZ+$j] 1 1.0"
    }

    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
    puts $fileModel "}"
} 

#pattern  Plain       2    "Linear"  { 
#puts $fileModel "pattern  Plain       2  Linear {"
#    # Load    nodeTag    LoadValues 
#    for { set j 0 } { $j <= $numY} {incr j 1} {
#      # load on 2 direction
#      #load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      #puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#      #load on 1 direction
#      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#    }
#    for { set j 1 } { $j <= $numX} {incr j 1} {
#      #load on 2 direction
#      #load [expr 1+($numX+1)*$numZ+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      #puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000 +1.000000E+000   0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#      #load on 1 direction
#      load [expr 1+($numX+1)*$numZ+$j]  +1.000000E+000 +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +1.000000E+000  +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#    }
#    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
# 
#    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
#puts $fileModel "}"
#} 
close $fileModel 
set  r_1 [recorder  Node  -file  output_Node_DefoShape_Dsp.dis  -time -nodeRange 1  $endNodeTag -dof  1  2  3 disp]
#recorder  Element -file output_Ele_All.stress -time -eleRange 1 $endEleTag stress
#recorder  Element -file output_Ele_All.strain -time -eleRange 1 $endEleTag strain
set  r_2 [recorder Node -xml output_disp_11.xml -nodeRange 1  $endNodeTag -time -dof 1 2 3 4 5 6 disp]
#set  r_2 [recorder Node -file output_accel_11.txt -nodeRange 1  $endNodeTag  -time -dof 1 2 accel]
set  r_4 [recorder Element -xml output_stress_1_11.xml -time  -eleRange 1 $endEleTag stresses]
#set  r_5 [recorder Element -file output_stress_2_11.txt -time  -eleRange 1 $endEleTag material 2 forces]
#set  r_6 [recorder Element -file output_stress_3_11.txt -time  -eleRange 1 $endEleTag material 3 forces]
#set  r_7 [recorder Element -file output_stress_4_11.txt -time  -eleRange 1 $endEleTag material 4 forces]
set  r_8 [recorder Element -xml output_strain_1_11.xml -time  -eleRange 1 $endEleTag strains]
#set  r_9 [recorder Element -file output_strain_2_11.txt -time  -eleRange 1 $endEleTag material 2 deformations]
#set r_10 [recorder Element -file output_strain_3_11.txt -time  -eleRange 1 $endEleTag material 3 deformations]
#set r_11 [recorder Element -file output_strain_4_11.txt -time  -eleRange 1 $endEleTag material 4 deformations]
# ---------------------------------- 
# perform the analysis
# ---------------------------------- 
set IDctrlNode [expr 1+($numX+1)*$numZ];
set IDctrlDOF 1;
set Dincr 0.1;
# Constraint Handler Lagrange Transformation
constraints  Plain
# Convergence Test 
test  NormDispIncr  +1.000000E-002    250     0
#test EnergyIncr  1.0e-12    10         0
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
puts "integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr"
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  SparseGeneral 
# Analysis Type 
analysis  Static 

set rcdcmd "recorder  Node  -file  PushoverCase_KeyNode_Dsp.out  -time -node [expr 1+($numX+1)*$numZ] -dof $IDctrlDOF  disp "
eval $rcdcmd;

#  perform Static Cyclic Displacements Analysis   /////////////////////////////////////////////////////////////

set iDmax1 [expr 0.08*$in/$mm];  # 2.032mm
set iDmax2 [expr 0.12*$in/$mm];  # 3.048mm
set iDmax3 [expr 0.3*$in/$mm];  #  7.62mm
set iDmax4 [expr 0.4*$in/$mm];  # 10.16mm
set iDmax5 [expr 0.6*$in/$mm];  # 15.24mm
set iDmax6 [expr 1.1*$in/$mm];  # 27.94mm

# perform the analysis    /////////////////////////////////////////////////////////////
#                0.08   -0.08  0.12   -0.12   0.3   -0.3   0.4   -0.4   0.6   -0.6   1.1   -1.1
#set numSteps     {   508   1016  1270    1524  2667   3810  4445   5080  6350   7620 10795  13970 }
set numIters     {     0    200   400     400   400    400   400    400   400    400   400    400 }
#set increments   { 0.004 -0.004 0.004  -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 }
set numSteps     { 20   40   50   60  106  152 178  204 254  304 432  560 }
set increments   { 0.1 -0.1 0.1  -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 }
	
for { set i 0 } { $i<12 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set Dincr [lindex $increments $i]
    integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
    if {$numIter == 0} {
	set numIter 100
    }
    test NormDispIncr 1.0e-3 $numIter 5
    #analyze $numStep

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
        set putout "PROBLEM: $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF]";
        puts $putout
        test NormDispIncr 1e-3 $numIter 5
        set ok [analyze 1]
      }; # end if
      incr j 1;
      #record;
      #puts $end
    }

    puts $i

}
#remove recorders;

