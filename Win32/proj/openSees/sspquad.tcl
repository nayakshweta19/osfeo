#########################################################
#					                #
# Coarse-mesh cantilever beam analysis.  The beam is    #
# modeled with only 4 elements and uses anti-symmetry.  #
#							#
# ---> Basic units used are kN and meters		#
#							#
#########################################################
 
wipe
 
model BasicBuilder -ndm 2 -ndf 2
 
# beam dimensions
set L 600.0
set D 600.0
 
# define number and size of elements 
set nElemX 2
set nElemY 2
set nElemT [expr $nElemX*$nElemY]
set sElemX [expr $L/$nElemX]
set sElemY [expr $D/$nElemY]
 
set nNodeX [expr $nElemX + 1]
set nNodeY [expr $nElemY + 1]
set nNodeT [expr $nNodeX*$nNodeY]
 
# create the nodes
set nid   1
set count 0.0
for {set j 1} {$j <= $nNodeY} {incr j 1} {
	for {set i 1} {$i <= $nNodeX} {incr i 1} {
		node $nid [expr 0.0 + $count*$sElemX] [expr ($j-1)*$sElemY]
		set nid   [expr $nid + 1]
		set count [expr $count + 1]
	}
	set count 0.0
}
 
# boundary conditions
fixX 0.0   1 1
#fix [expr $nElemY*$nNodeX+1]   1 0
#for {set k 2} {$k <= $nNodeX} {incr k 1} {
#	fix $k 1 0
#}
 
# define material
set matID 1
set E     2e5
set nu    0.2
nDMaterial ElasticIsotropic $matID $E $nu

# set fc fy E
set wfc 34.5;  #[expr (6880.+6600.+7700.+6850.)/4.0*$psi/$MPa];
set wfy 414.0; #[expr (76.+63.+60)/3.0*$ksi/$MPa];
set wE 2.0e5;  #[expr (29000+28000+28700)/3.0*$ksi/$MPa];
set rou1 0.0059;
set rou2 0.0226;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
#uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
#uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2


# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 

set pi 3.141592654;
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag    rho   s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 2  2.7e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
nDMaterial CSMMRCPlaneStress 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025
set matID 2
# create elements
set thick 150.
set b1 0.0
set b2 0.0
set count 1
for {set j 1} {$j <= $nNodeY} {incr j 1} {
	for {set i 1} {$i <= $nNodeX} {incr i 1} {
		if {($i < $nNodeX) && ($j < $nNodeY)} {
			set nI [expr $i+($j-1)*$nNodeX]
			set nJ [expr $i+($j-1)*$nNodeX+1]
			set nK [expr $i+$j*$nNodeX+1]
			set nL [expr $i+$j*$nNodeX]
			element SSPquad $count  $nI $nJ $nK $nL $matID "PlaneStress" $thick $b1 $b2
# 			element quad    $count  $nI $nJ $nK $nL $thick "PlaneStress" $matID 0.0 2e-7 $b1 $b2
			set count [expr $count+1]
		}
	}
}
 
# create recorders
set step 0.1
 
recorder Node -time -file results/d1p1m1.out -dT $step -nodeRange 1 $nNodeT -dof 1 2 disp
recorder Element -eleRange 1 $nElemT -time -file results/s1p1m1.out  -dT $step  stress
recorder Element -eleRange 1 $nElemT -time -file results/e1p1m1.out  -dT $step  strain

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp 0 0 1000
  #vrp 0 -500 250
  vup 1 0 0
  #vpn -1 -1 0.5
  viewWindow -400 800 -400 800
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 1

# create loading
set P -300.0;
 
pattern Plain 3 {Series -time {0 10 15} -values {0 1 1} -factor 1} { 
	load [expr $nNodeX*1]   0.0 [expr 0.1875*$P]
	load [expr $nNodeX*2]   0.0 [expr 0.3125*$P]
	load [expr $nNodeX*3]   0.0 [expr 0.1875*$P]
	#load $nNodeT   0.0 [expr 0.1875*$P]
	#load $nNodeX   0.0 [expr 0.3125*$P]
	#load [expr $nNodeX+1]   0.0 [expr -0.1875*$P]
}
 
# create analysis
 
integrator LoadControl 0.1
numberer RCM
system SparseGeneral
constraints Transformation
test NormDispIncr 1e-3 200 0
algorithm Newton
analysis Static
 
analyze 250
 
#wipe