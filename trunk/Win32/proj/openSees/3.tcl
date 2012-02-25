wipe;
model BasicBuilder -ndm 2 -ndf 2;
set dataDir Data;
file mkdir $dataDir;
source libUnits.tcl;

set pi 3.141592654;
set nDMatTag 6;
set PlaneStressMatTag 20;
set testMatTag 101;
set testSecTag 201;
#nDMaterial ElasticIsotropic $matTag $E $v
nDMaterial ElasticIsotropic 1 1.0e4 0.2;
nDMaterial PlaneStress 2 1;
nDMaterial PlaneStrain 3 1;

set k       27777.78; #---bulk modulus
set G       9259.26; #---shear modulus
set sig0    30.0; #---yield stress
set sigInf   60.5; #---final saturation yield stress
set delta    0.1; #---exponential hardening parameter eta>= 0 for rate independent case
set H        0.05; #---linear hardening parameter
#nDMaterial J2Plasticity $nDMatTag $k $G $sig0 $sigInf $delta $H

set k       27777.78; #---bulk modulus
set G       9259.26; #---shear modulus
set sigY    5.0; #---yield stress
set rho     0.398; #---failure surface and associativity
set rhoBar  0.398; 
set Kinf    0.0; #---isotropic hardening
set Ko 	    0.0
set delta1  0.0
set H 	    0.0; #---kinematic hardening
set theta   1.0
set delta2  0.0; #---tension softening
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
#nDMaterial DruckerPrager  $nDMatTag 7.8e-6 $k $G $sigY $rho $rhoBar $Kinf $Ko $delta1 $delta2 $H $theta

##nDMaterial PlaneStress $PlaneStressMatTag $nDMatTag

set wfc 30.0;
set wfyv 300.0;
set wfyh1 450.0;
set wE 1.96e5;
set epsc0 0.002;
set fpcu  46.0;
set epsU  0.025;
set lambda 0.1;
set ft  [expr 0.14*$wfc];
set Ets [expr $ft/0.002];
set rou1 0;
set rou2 0;
set db1 16;
set db2 18;
set b 0.02
set rouv 0.0152;
set rouh1 0.0167;

uniaxialMaterial MCFTSteel03 11 $wfyv $wE $b
uniaxialMaterial MCFTSteel03 12 $wfyh1 $wE $b

#uniaxialMaterial MCFTConcrete01 14 [expr -$wfc] -0.0022
#uniaxialMaterial MCFTConcrete01 15 [expr -$wfc] -0.0022
uniaxialMaterial MCFTConcrete03 14 [expr -$wfc] -$epsc0 $fpcu $epsU $lambda $ft $Ets
uniaxialMaterial MCFTConcrete03 15 [expr -$wfc] -$epsc0 $fpcu $epsU $lambda $ft $Ets

nDMaterial MCFTRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.
#nDMaterial DSFMRCPlaneStress $PlaneStressMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $db1 $db2 $wfc $wfyv $wfyh1 $wE 0.0025 10 75 75.

nDMaterial PlaneStressFiber $testMatTag $PlaneStressMatTag;

# --------------------------------------------------------------------------------------------------
node 1001 0.0 0.0;
node 1002 1000.0 0.0;
node 1003 1000.0 1000.0
node 1004 0.0 1000.0;

# Fix all degrees of freedom except axial and bending
fix 1001 1 1 ; #1
fix 1002 0 1 ; #1

# Define element
#                         tag ndI ndJ  secTag
#element zeroLengthSection  2001   1001   1002  $testSecTag
element quad 2001 1001 1002 1003 1004 100.0  "PlaneStress" $PlaneStressMatTag
# Create recorder
recorder Node -file disp.out -time -node 1003 1004 -dof 1 2 disp;
set numElem 1;
recorder Element -eleRange 1 $numElem -time -file stress.out  stress
recorder Element -eleRange 1 $numElem -time -file strain.out  strain

# Define constant axial load
set axialLoad -100;
pattern Plain 3001 "Constant" {
    #load 1002 $axialLoad 0.0 0.0
    load 1003   0.0 $axialLoad 
    load 1004   0.0 $axialLoad 
}

# Define analysis parameters
integrator LoadControl 0 1 0 0
system UmfPack;	# Overkill, but may need the pivoting!SparseGeneral -piv
#test EnergyIncr  1.0e-4 10
test NormDispIncr 1.0e-6 250
numberer Plain
constraints Plain
algorithm Newton
analysis Static