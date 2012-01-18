# -------------------------------------------------------------------------
# Create ModelBuilder for 2D element (with two-dimensions and 2 DOF/node)
# -------------------------------------------------------------------------
 wipe
model basic -ndm 2 -ndf 2
set L 1000;
set H 1000;
set deltaL 1000;
set deltaH 1000;
set nL 1; #[expr $L/$deltaL];
set nH 1; #[expr $H/$deltaH];
set t 350;
# Create nodes & add to Domain - command: node nodeId xCrd yCrd
 
set j 0;
while {$j < [expr $nH+1]} {
    set i 0;
    while {$i < [expr $nL+1]} {
	node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]
	puts "node [expr $j*($nL+1)+$i+1] [expr $i*$deltaL] [expr $j*$deltaH]"
	set i [expr $i+1]
    }
    set j [expr $j+1]
}
 
# fix one end as a pin, the other end as a roller
fix 1 1 1
fix 2 1 1

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 30.0;
set wfy 235.0;
set wE 2.06e5;
set rou1 0.003;
set rou2 0.003;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
#uniaxialMaterial Steel01 $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
#uniaxialMaterial Steel01 11 $wfy  $wE  0.02
#uniaxialMaterial Steel01 12 $wfy  $wE  0.02

# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial Concrete01 $matTag $fpc $epsc0 $fpcu $epsU
uniaxialMaterial Concrete01 13 [expr -$wfc] -0.0025 [expr -1.3*$wfc] -0.003
uniaxialMaterial Concrete01 14 [expr -$wfc] -0.0025 [expr -1.3*$wfc] -0.003

set pi 3.141592654
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
nDMaterial FAReinforcedConcretePlaneStress 15  0.005 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
puts "nDMaterial FAReinforcedConcretePlaneStress 15  0.05 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002"

# -----------------------------------------------------------------
#  Define 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------

set j 0;
while {$j < $nH} {
    set i 0;
    while {$i < $nL} {
     # Create quad elements - command:
     # element quad eleID node1 node2 node3 node4 thick  type         matID
 
	element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15
	puts "element quad [expr $j*$nL+$i+1] [expr $j*($nL+1)+$i+1] [expr $j*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+2] [expr ($j+1)*($nL+1)+$i+1] $t PlaneStress 15"
	set i [expr $i+1]
    }
    set j [expr $j+1]
}
 
#4.3
# ---------------------------
# Define horizontal loads
# ---------------------------
 
set N 5000.0;
 
# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
 
    # Create the nodal load - command: load nodeID xForce yForce
    load 3 [expr  $N] 0 
    load 4 [expr  $N] 0 
}
 
#4.4
# ------------------------------
# End of model generation
# ------------------------------
 
# ------------------------------
# Start of analysis generation
# ------------------------------
 
# Create the system of equation, a sparse solver with partial pivoting
  system BandGeneral
 
 
# Create the constraint handler
  constraints Plain
 
# Create the DOF numberer
  numberer Plain
 
 
# Create the convergence test
  test NormDispIncr 1.0e-4 20 1
 
 
# Create the solution algorithm
  algorithm KrylovNewton
 
# Create the integration scheme, the DisplacementControl scheme
  integrator LoadControl 0.1
 
# Create the analysis object
  analysis Static
 
# initialize in case we need to do an initial stiffness iteration
  initialize
 
# ------------------------------
# End of analysis generation
# ------------------------------
 
# Create a recorder to monitor nodal displacements
recorder Node -file e1.out -time -node 3 4 -dof 1 2 3 disp
 
# perform the analysis
set ok [analyze 10]

set IDctrlNode 3
set IDctrlDOF 1
set controlDisp 0.001
set Dmax 0.1
set Tol 1e-5

if {$ok != 0} {
	# if analysis fails, we try some other stuff, performance is slower inside this loop
	set Dstep 0.0;
	set ok 0
	while {$Dstep <= 1.0 && $ok == 0} {
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr $controlDisp/$Dmax]
		set ok [analyze 1 ]
		# if analysis fails, we try some other stuff
		# performance is slower inside this loop	global maxNumIterStatic;	
		# max no. of iterations performed before "failure to converge" is ret'd
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 2000 0
			algorithm Newton -initial
			set ok [analyze 1]
			test $testType $Tol      $maxNumIter    0
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 ]
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch 0.8
			set ok [analyze 1]
			algorithm $algorithmType
		}
	};	# end while loop
};      # end if ok !0
