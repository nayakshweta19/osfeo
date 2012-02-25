# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source libUnits.tcl;

# MATERIAL parameters -------------------------------------------------------------------

set nDMatTag 1;
# Material "Steel":    matTag    E    v    rho 
#nDMaterial  ElasticIsotropic3D $nDMatTag 5.e5 0.2 78.5
#---bulk modulus
	set k       27777.78
#---shear modulus
	set G       9259.26
#---yield stress
	set sig0    3.0
#---final saturation yield stress
        set sigInf   6.5
#---exponential hardening parameter eta>= 0 for rate independent case
        set delta    0.1
#---linear hardening parameter
        set H        0.05
nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

set nDBeamTag 11;

#--define material parameters for the model
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
#nDMaterial DruckerPrager $nDMatTag 7.8 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

#nDmaterial for beam fiber 3D
#nDMaterial BeamFiber $nDBeamTag $nDMatTag;
#                Dodd_Restrepo $tag $Fy $Fsu $ESH $ESU $Youngs $ESHI $FSHI <$OmegaFac $Conv>
#uniaxialMaterial Dodd_Restrepo $nDBeamTag 200. 275. 0.002 0.015 2.e5 0.005 245.
#                         $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
#uniaxialMaterial Steel01 $nDBeamTag 500 0.002 0.03;
	
set beamSecTag 1;
#section definition 
#section Fiber $beamSecTag { ;
#	fiber -10 0 1 $nDBeamTag;
#	fiber  10 0 1 $nDBeamTag;
#};
set uDBeamTag 2;
nDMaterial BeamFiber2d $uDBeamTag $nDMatTag;
#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
#section WSection2d $beamSecTag $nDBeamTag  300 10 120 20 16 12
section Timoshenko $beamSecTag { ;
	fiber2d -10 0 10 $uDBeamTag;
	fiber2d  -5 0 10 $uDBeamTag;
	fiber2d   0 0 10 $uDBeamTag;
	fiber2d   5 0 10 $uDBeamTag;
	fiber2d  10 0 10 $uDBeamTag;
};

# --------------------------------------------------------------------------------------------------
set iNode 1001
set jNode 1002
# Define two nodes
node 1001 0.0 0.0;
node 1002 0.0 200.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 1
fix 1002 0 0 0

set transfTag 1;
#geomTransf Timoshenko $transfTag
geomTransf Linear $transfTag;
set np 5; # int. points
# Define element
set eleTag 2001
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $secTag
#forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag
#element nonlinearBeamColumn 2001 1001 1002 $np $beamSecTag $transfTag;
#element dispBeamColumn 2001 1001 1002 $np $beamSecTag $transfTag;
element Timoshenko2d01 2001 1001 1002 $np $beamSecTag $transfTag;
#element beamWithHinges $eleTag $iNode $jNode $secTagI $Lpi $secTagJ $Lpj $E $A $Iz $Iy $G $J $transfTag
set Lpi 0.16;
set Lpj 0.16;
set E 2.100000E+011;
set A 2.480000E-003;
set Iz 4.510000E-005;
#element beamWithHinges $eleTag $iNode $jNode $beamSecTag $Lpi $beamSecTag $Lpj $E $A $Iz $transfTag

# Create recorder
recorder Node -file data/Mphi.out -time -node 1002 -dof 1 2 3 disp;	# output moment (col 1) & curvature (col 2)
recorder Node -file data/Reaction.out -time -node 1001 -dof 1 2 3 reaction;

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -750 -800 800
  #vrp 0 -500 250
  vup 0 1 0
  #vpn -1 -1 0.5
  viewWindow -400 400 -400 400
}
port -1 1 -1 1
projection 1
fill 0
display 1 -1 1

## Define constant axial load
#set axialLoad 1000;
#pattern Plain 3001 "Constant" {;
#	load 1002 0.0 0.0 $axialLoad 0.0 0.0 0.0;
#};
#
## Define analysis parameters
integrator LoadControl 0 1 0 0;
system SparseGEN;	# Overkill, but may need the pivoting!
test EnergyIncr  1.0e-9 10;
numberer Plain;
constraints Plain;
algorithm Newton;
analysis Static;
#analyze 1
## Do one analysis for constant axial load
# Compute curvature increment
set maxK 50.0;
set numIncr 100;
set dK [expr $maxK/$numIncr];

set IDctrlNode 1002;
set IDctrlDOF 1;
set Dmax $maxK;
set Dincr $dK;
set TolStatic 1.e-9;
set testTypeStatic EnergyIncr;
set maxNumIterStatic 6;
set algorithmTypeStatic Newton;

# Define reference moment
pattern Plain 3002 "Linear" {;
	load $IDctrlNode 1.0 0.0 0.0;
};

# Use displacement control at node 1002 for section analysis, dof 3
integrator DisplacementControl $IDctrlNode $IDctrlDOF $dK 1 $dK $dK

# Do the section analysis
set ok [analyze $numIncr]
