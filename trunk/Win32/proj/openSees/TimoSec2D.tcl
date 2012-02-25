wipe;					# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

# MATERIAL parameters -------------------------------------------------------------------

set nDMatTag 1;
nDMaterial  ElasticIsotropic3D $nDMatTag 5.e5 0.3 0.0

set k       27777.78; #---bulk modulus
set G       9259.26;  #---shear modulus
set sig0    3.0;      #---yield stress
set sigInf  6.5;      #---final saturation yield stress
set delta   0.1;      #---exponential hardening parameter eta>= 0 for rate independent case
set H       0.05;     #---linear hardening parameter
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

set nDBeamTag 11;
set k       27777.78; #---bulk modulus
set G       9259.26;  #---shear modulus
set sigY    5.0;      #---yield stress
set rho     0.398;    #---failure surface and associativity
set rhoBar  0.398;    #
set Kinf    0.0;      #---isotropic hardening
set Ko 	    0.0
set delta1  0.0
set H 	    0.0;      #---kinematic hardening
set theta   1.0
set delta2  0.0;      #---tension softening
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
#nDMaterial DruckerPrager $nDMatTag 7.8 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

#nDMaterial BeamFiber $nDBeamTag $nDMatTag;

# back-bone steel stress-strain curve parameters for all materials
set IDMat 101;
set Fy [expr 60*$ksi];         # STEEL yield stress
set Es [expr 29000*$ksi];      # modulus of steel
set epsY [expr $Fy/$Es];       # steel yield strain
set epsSH [expr 8*$epsY];
set Esh [expr 0.02*$Es];       # tangent stiffness at onset of StrainHardening
set Bs 0.01;                   # strain-hardening ratio
set Fy1 [expr 1.01*$Fy];       # steel stress post-yield
set epsY1 [expr $epsY+($Fy1-$Fy)/($Bs*$Es)];    # steel strain post-yield
set Fu [expr 1.16*$Fy];
set epsU 1.0;                   # ultimate strain of steel#

set CompressiveStress [expr -5.35*$ksi];        # nominal concrete compressive strength
set fc [expr $CompressiveStress];
set Ec [expr 57000.*$psi*sqrt(abs($fc)/$psi)];

set fclU $fc;              # UNCONFINED concrete (todeschini parabolic model), maximum stress
set epslU -0.0022;         # strain at maximum strength of unconfmed concrete
set fc2U [expr 0.2*$fclU]; # ultimate stress
set eps2U -0.01;           # strain at ultimate stress
set lambda 0.1;            # ratio between unloading slope at $eps2 and initial slope $Ec
set ftU [expr -0.14*$fclU]; # tensile strength +tension
set Ets [expr $ftU/0.002];  # tension softening stiffness

# hysteretic steel material parameters -- baseline
set pinchX 0.6;              # pinching parameter for hysteretic model
set pinchY 0.7;              # pinching parameter for hysteretic model
set damage1 0.1;             # damage parameter for hysteretic model
set damage2 0.2;             # damage parameter for hysteretic model
set betaMUsteel 0.4;         # degraded unloading stiffness for hysteretic material based on MUA(-beta)
# steel02 and steel03 parameters ¡ª baseline
set R0 18;                   # control the transition from elastic to plastic branches,
set cR1 0.925;               # control the transition from elastic to plastic branches,
set cR2 0.15;                # control the transition from elastic to plastic branches
set a2 0.1;                  # isotropic hardening parameter, associated with al
set a1 [expr $a2*($Fy/$Es)]; # isotropic hardening parameter, increase of comp.
                             # yield envelope as proportion of yield strength after a plastic strain
set a4 0.1;                  # ft isotropic hardening parameter, associated with a3
set a3 [expr $a4*($Fy/$Es)]; # isotropic hardening parameter, increase of tension yield envelope
                             # as proportion of yield strength after a plastic strain

set rD1 [expr 0.0*$PI];      # Reinforcing direction, D, and ratio, R, and diameter, d.
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

nDMaterial NonlinearBS $IDMat $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $rR1 $rdmm1 $rD2 $rR2 $rdmm2 $Fy $epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

set pi 3.141592654;
set PlaneStressMatTag 20;
set wfc  30.0; #19.0;
set wfyv 458.0;
set wfyh1 300.0;
set wE 2.e5;
set epsc0 0.00215;
set fpcu  46.0;
set epsU  0.025;
set lambda 0.1;
set ft  1.72; #[expr 0.14*$wfc];
set Ets [expr $ft/0.002];
set rou1 0;
set rou2 0;
set db1 16;
set db2 18;
set b 0.02
set rouv 0.01785;
set rouh1 0.00713;

uniaxialMaterial MCFTSteel03 11 $wfyv $wE $b
uniaxialMaterial MCFTSteel03 12 $wfyh1 $wE $b
uniaxialMaterial MCFTConcrete03 14 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets
uniaxialMaterial MCFTConcrete03 15 [expr -$wfc] -$epsc0 -$fpcu -$epsU $lambda $ft $Ets

nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial DSFMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

uniaxialMaterial SteelZ01 101 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 102 $wfyh1 $wE $wfc $rouh1

#uniaxialMaterial ConcreteL02 104 [expr -$wfc] -0.003
#uniaxialMaterial ConcreteL02 105 [expr -$wfc] -0.003
uniaxialMaterial ManderConcreteB01 104 [expr -$wfc] -0.0025 30000. 0.40 0.60 [expr 0.1*$wfc] 8.0e-5 0.5 2.9 -spall 2.1
uniaxialMaterial ManderConcreteB01 105 [expr -$wfc] -0.0025 30000. 0.40 0.60 [expr 0.1*$wfc] 8.0e-5 0.5 2.9 -spall 2.1

#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
nDMaterial CSMMRCPlaneStress 22 0.0 101 102 104 105 [expr 0.5*$pi] [expr 0.0*$pi] $rouv $rouh1 $wfc $wfyv $wE -0.002

set beamSecTag 1;
#section definition 
#section Fiber $beamSecTag { ;
#	fiber -10 0 1 $nDBeamTag;
#	fiber  10 0 1 $nDBeamTag;
#};

set uDBeamTag 2;
nDMaterial PlaneStressRCFiber $uDBeamTag 22; #$IDMat; #$PlaneStressMatTag; #$nDMatTag;

#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
#section WSection2d $beamSecTag $nDBeamTag  300 10 120 20 16 12
section Timoshenko $beamSecTag { ;
	fiber2d -10 0 100 $uDBeamTag;
	fiber2d  -5 0 100 $uDBeamTag;
	fiber2d   0 0 100 $uDBeamTag;
	fiber2d   5 0 100 $uDBeamTag;
	fiber2d  10 0 100 $uDBeamTag;
};

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0
node 1002 0.0 0.0

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 1
fix 1002 0 0 0

# Define element
#                         tag ndI ndJ  secTag
element zeroLengthSection  2001   1001   1002  $beamSecTag

# Create recorder
recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)

# Define constant axial load
set axialLoad -5000;
pattern Plain 3001 "Constant" {
	load 1002 $axialLoad 0.0 0.0
}

# Define analysis parameters
integrator LoadControl 0.01
system SparseGeneral -piv;	# Overkill, but may need the pivoting!
test EnergyIncr  1.0e-9 10
numberer Plain
constraints Plain
algorithm Newton
analysis Static

# Do one analysis for constant axial load
analyze 100
loadConst 0.0

# Define reference moment
pattern Plain 3002 "Linear" {
	load 1002 0.0 1.0 0.0 
}
# Compute curvature increment
set maxK 10.0;
set numIncr 1000000;
set dK [expr $maxK/$numIncr];

# Use displacement control at node 1002 for section analysis, dof 3
integrator DisplacementControl 1002 2 $dK 1 $dK $dK

# Do the section analysis
set ok [analyze $numIncr]

# ----------------------------------------------if convergence failure-------------------------
set IDctrlNode 1002
set IDctrlDOF 3
set Dmax $maxK
set Dincr $dK
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr  
set maxNumIterStatic 6
set algorithmTypeStatic Newton

global LunitTXT;					# load time-unit text
if {  [info exists LunitTXT] != 1} {set LunitTXT "Length"};		# set blank if it has not been defined previously
	
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Curv=%.4f /%s";	# format for screen/file output of DONE/PROBLEM analysis

if {$ok != 0} {  
	# if analysis fails, we try some other stuff, performance is slower inside this loop
	set Dstep 0.0;
	set ok 0
	while {$Dstep <= 1.0 && $ok == 0} {	
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr $controlDisp/$Dmax]
		set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
		if {$ok != 0} {;				# reduce step size if still fails to converge
			set Nk 4;			# reduce step size
			set DincrReduced [expr $Dincr/$Nk];
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
				set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
				if {$ok != 0} {  
					# if analysis fails, we try some other stuff
					# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
					puts "Trying Newton with Initial Tangent .."
					test NormDispIncr   $TolStatic      2000 0
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
				if {$ok != 0} {;				# stop if still fails to converge
					puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					return -1
				}; # end if
			}; # end for
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;	# bring back to original increment
		}; # end if
	};	# end while loop
};      # end if ok !0
# -----------------------------------------------------------------------------------------------------
    
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}