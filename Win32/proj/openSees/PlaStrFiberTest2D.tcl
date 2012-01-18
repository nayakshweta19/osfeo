wipe;
model basic -ndm 2 -ndf 3;

set pi 3.141592654;
set nDMatTag 6;
set PlaneStressMatTag 20;
set testMatTag 101;
#nDMaterial ElasticIsotropic $matTag $E $v
nDMaterial ElasticIsotropic 1 10e3 0.2;
nDMaterial PlaneStress 2 1;
nDMaterial PlaneStrain 3 1;

#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sig0    30.0
#---final saturation yield stress
        set sigInf   60.5
#---exponential hardening parameter eta>= 0 for rate independent case
        set delta    0.1
#---linear hardening parameter
        set H        0.05
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sigY    5.0
#---failure surface and associativity
	set rho     0.398
	set rhoBar  0.398
#---isotropic hardening
	set Kinf    0.0
	set Ko 	    0.0
	set delta1  0.0
#---kinematic hardening
	set H 	    0.0
	set theta   1.0
#---tension softening
	set delta2  0.0
 
#--material models
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
#nDMaterial DruckerPrager  $nDMatTag 7.8e-6 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

##nDMaterial PlaneStress $PlaneStressMatTag $nDMatTag

# set fc fy E
set wfc 76.0;
set wfyv 500.0;
set wfyh1 450.0;
set wE 1.96e5;
set rou1 0;
set rou2 0;
set rouv 0.0352;
set rouh1 0.0367;

# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteL02 14 [expr -$wfc] -0.003
uniaxialMaterial ConcreteL02 15 [expr -$wfc] -0.003

# NDMaterial: FAFourSteelPCPlaneStress
#                                                  tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 19:52 2011-10-30rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelPCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteetPCPlaneStress           matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteetPCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAFourSteelRCPlaneStress  (100)
#                   $PlaneStressMatTag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelRCPlaneStress $PlaneStressMatTag 0.0 11 12 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 $rouv $rouh1 $wfc 0.0 $wfyv $wE 0.003

#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
nDMaterial CSMMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002

#nDMaterial MCFTRCPlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi] [expr 0.0*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002 10 140 140.

nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

set np 4; # int. points
#set C 0.4; # center of rotation 

#section definition 
#section CSMMFiber2d 3 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 21;
#	fiber2d  10 0 1 21;
#};
section Timoshenko 3 { ;
	layer2d straight $testMatTag 5 5376.33 -116 0 116 0
	#fiber2d -100 0 5376 $testMatTag;
	#fiber2d  -50 0 5376 $testMatTag;
	#fiber2d    0 0 5376 $testMatTag;
	#fiber2d   50 0 5376 $testMatTag;
	#fiber2d  100 0 5376 $testMatTag;
};
#section FiberInt 4 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 2;
#	fiber2d  10 0 1 2;
#};

#nodes
set n 2;
node 1 0 0.0;
node 2 0 1016;
fix 1 1 1 1;
fix 2 0 0 1;
mass 1 1e6 1e6 0;
mass 2 1e6 1e6 0;

# GeoTran    type    tag 
#geomTransf  Timoshenko  1;
geomTransf PDelta 1
#geomTransf LinearInt 1
#element definition
#element Timoshenko 1 1 2 $np 3 1;
#element dispBeamColumnInt 1 1 2 1 4 1 $C
#element dispBeamColumn 1 1 2 $np 3 1
element nonlinearBeamColumn 1 1 2 $np 3 1
set Lpi 0.2;
set Lpj 0.2;
set E 3.03e4;
set A 64516.0;
set Iz 346859521;
set Iy 346859521;
set G 1.01e4;
set J 6.88e9;
#element beamWithHinges 1 1 2 3 $Lpi 3 $Lpj $E $A $Iz 1

# Set axial load 
pattern Plain 1 Constant {;
	load 2   0.0 -489e3 0.0;
};

initialize;
integrator LoadControl 0 1 0 0;
system SparseGEN; #ItpackSparseGEN
# Convergence Test 
test  NormDispIncr  +1.000000E-004    250     0
#test EnergyIncr  1.0e-12    10         0
numberer Plain;
constraints Plain;
algorithm ModifiedNewton -initial;
analysis Static;

# perform the gravity load analysis, 
analyze 1

## Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0;
remove recorders;

# set lateral load 
pattern Plain 2 Linear { ;
	load 2 1.0 0.0 0.0;
} 

# Create a recorder to monitor nodal displacement and element forces
recorder Node -file nodeTop.out -node 2 -dof 1 disp
recorder Element -file elesX.out -time -ele 1 globalForce

set IDctrlNode 2;
set IDctrlDOF 1;
set Dincr 0.05;
# Constraint Handler Lagrange Transformation
constraints  Plain
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  SparseGEN
# Analysis Type 
analysis  Static 

# perform the analysis    /////////////////////////////////////////////////////////////
#                0.08   -0.08  0.12   -0.12   0.3   -0.3   0.4   -0.4   0.6   -0.6   1.1   -1.1
set numIters     {   200    200   400     400   400    400   400    400   400    400   400    400 }
#set numSteps     {   508   1016  1270    1524  2667   3810  4445   5080  6350   7620 10795  13970 }
#set numSteps     { 20   40   50   60  106  152 178  204 254  304 432  560 }
set numSteps    { 20 40 60 80 100 160 240 320 400 500 600 750 1000}
#set increments   { 0.004 -0.004 0.004  -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 }
set increments   { 0.1 -0.1 0.1  -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 }
#set increments   { 1 -1 1  -1 1 -1 1 -1 1 -1 1 -1 1 }

for { set i 0 } { $i<12 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set Dincr [lindex $increments $i]
    integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
    puts "integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr, numSteps $i"
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
      puts "nodeDisp 2 1 [nodeDisp 2 1]"
    }

    puts $i
}
#remove recorders;

