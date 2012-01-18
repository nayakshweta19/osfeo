# ----------------------------
# Start of model generation
# ----------------------------
wipe;
source Units&Constants.tcl;
# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 3 -ndf 3

# create the material
nDMaterial ElasticIsotropic   1   1000   0.25  6.75 

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
uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
 
set pi 3.141592654
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag    rho   s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
nDMaterial FAReinforcedConcretePlaneStress 2  2.3e-3 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

# Define geometry
# ---------------

# define some  parameters

set mass 1.0;

set Quad  quad3d;

set eleArgs "PlaneStress2D";
set thick [expr $in*6.0/$mm];

set xWidth [expr 7.8E+1*$in/$mm];
set yWidth [expr 5.4E+1*$in/$mm];
set height [expr 2.88E+2*$in/$mm];

set numX 3; # num elements in x direction
set numY 4; # num elements in y direction ,divied by 2, for middle wall
set numZ 8; # num elements in z direction ,divied by 4, for 4 story

set nodeNum 1

set incrX [expr $xWidth/(1.0*$numX)]
set incrY [expr $yWidth/(1.0*$numY)]
set incrZ [expr $height/(1.0*$numZ)]

set yLoc [expr 2.7E+1*$in/$mm];
set zLoc 0.0;

set fileModel [open "qmodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  3 \n"
#web wall
for {set i 0} {$i <= $numZ} {incr i 1} {
    set xLoc 0.0;
    for {set j 0} {$j <= $numX} {incr j 1} {
	node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass
	puts $fileModel "node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass"
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
    for {set i 0} {$i <$numX} {incr i 1} {
	#element quad3d $eleTag $iNode $jNode $kNode $lNode $thick "PlaneStress2D" 1
	element SSPquad $eleTag $iNode $jNode $kNode $lNode 1 "PlaneStress" $thick
	puts $fileModel "element quad $eleTag $iNode $jNode $kNode $lNode $thick \"PlaneStress2D\" 1"
	incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
	
    }
}

incr nodeNum -1

timeSeries Linear 1
pattern Plain 1 1 {
    load $nodeNum 10.0 0.0 0.0
}

for {set i 1} {$i <= $numX} {incr i 1} {
    fix $i 1 1 1
    puts $fileModel "fix $i 1 1 1"
}

for {set i [expr $numX+1]} {$i <= $nodeNum} {incr i 1} {
    fix $i 0 1 0
    puts $fileModel "fix $i 0 1 0"
}

set xLoc 0.0;
set zLoc 0.0;

incr nodeNum
set nodeYdirnStart $nodeNum

for {set i 0} {$i <= $numZ} {incr i 1} {
    set yLoc 0.0;
    for {set j 0} {$j <= $numY} {incr j 1} {
	node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass
	puts $fileModel "node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass"
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
    for {set i 0} {$i <$numY} {incr i 1} {
	#element quad3d $eleTag $iNode $jNode $kNode $lNode $thick "PlaneStress2D" 1
	element SSPquad $eleTag $iNode $jNode $kNode $lNode 1 "PlaneStress" $thick
	puts $fileModel "element quad $eleTag $iNode $jNode $kNode $lNode $thick \"PlaneStress2D\" 1"
	incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
    }
}

incr nodeNum -1

for {set i 1; set j $nodeYdirnStart} {$i <= $numY} {incr i 1; incr j 1} {
    fix $j 1 1 1
    puts $fileModel "fix $i 1 1 1"
}

for {set i [expr $nodeYdirnStart+$numY]} {$i <= $nodeNum} {incr i 1} {
    fix $i 1 0 0
    puts $fileModel "fix $i 1 0 0"
}

for {set i [expr 1 + $numX]; set j [expr $nodeYdirnStart+$numY]} {$i < $nodeYdirnStart} {incr i $numX; incr j $numY} {
    remove sp $i 2
    remove sp $j 1
    equalDOF $i $j 1 2 3
    puts $fileModel "remove sp $i 2"
    puts $fileModel "remove sp $j 1"
    puts $fileModel "equalDOF $i $j 1 2 3"
}

integrator LoadControl  1.0  1   1.0   10.0
test EnergyIncr  1.0e-12    10         0
algorithm Newton
numberer RCM
constraints Plain 
system ProfileSPD
analysis Static

# Perform the analysis
analyze   10     

# create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -17500 -58000 60000
  #vrp 0 -500 250
  vup 0 0 1
  #vpn -1 -1 0.5
  viewWindow -4000 4000 -4000 4000
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 10

# --------------------------
# End of recorder generation
# --------------------------

# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
setTime 0.0

# Now remove the loads and let the beam vibrate
remove loadPattern 1

# Create the transient analysis
test EnergyIncr  1.0e-12    10         0
algorithm Newton
numberer RCM
constraints Plain 
integrator Newmark 0.5 0.25
system ProfileSPD
#integrator GeneralizedMidpoint 0.50
analysis Transient

# Perform the transient analysis (20 sec)
#       numSteps  dt
analyze  500     0.5

