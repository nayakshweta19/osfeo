wipe;

model basic -ndm 2 -ndf 2;

node 1 0 0 -mass 0.125 0.125;
node 2 1 0 -mass 0.125 0.125;
node 3 1 1 -mass 0.125 0.125;
node 4 0 1 -mass 0.125 0.125;

#--define material parameters for the model
#---bulk modulus
	set k       27777.78;
#---shear modulus
	set G       9259.26;
#---yield stress
	set sigY    5.0;
#---failure surface and associativity
	set rho     0.398;
	set rhoBar  0.398;
#---isotropic hardening
	set Kinf    0.0;
	set Ko 	    0.0;
	set delta1  0.0;
#---kinematic hardening
	set H 	    0.0;
	set theta   1.0;
#---tension softening
	set delta2  0.0;
 
#--material models
#	   type	         tag  k   G   sigY   rho   rhoBar   Kinf   Ko   delta1   delta2   H   theta   
nDMaterial DruckerPrager 1    $k  $G  $sigY  $rho  $rhoBar  $Kinf  $Ko  $delta1  $delta2  $H  $theta;

#nDMaterial MultiAxialCyclicPlasticity 1 1 200 100 100 0.1 0.2 0.3 0.5 0.2 0.3 0.4 0.1 0.2;

nDMaterial PlaneStrain 2 1;

#nDMaterial PlaneStress 3 1;

element quad 1 1 2 3 4 0.2 PlaneStrain 2;
element nineNodeMixedQuad
#element quad 2 1 2 3 4 0.1 PlaneStress 3;

#fix 1 1 1;
#fix 2 1 1;


