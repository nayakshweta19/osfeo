wipe;
model basic -ndm 2 -ndf 3 
set E 2.0e5; 
uniaxialMaterial Elastic 1 $E;
set A 1;
set Iz 2;
set secTag 1;
section Elastic $secTag $E $A $Iz;
node 1 0 0.0;
node 2 0 1000.0;

mass 1 1e6 1e6 0;
mass 2 1e6 1e6 0;

fix 1 1 1 1;
fix 2 0 0 0;
	
# GeoTran    type    tag 
set transfTag 1;
#geomTransf  PDelta $transfTag;
geomTransf  Timoshenko $transfTag;
set np 3;
element Timoshenko 1 1 2 $np $secTag $transfTag;
#element dispBeamColumnInt 1 1 2 $np 2 1 $C
#element dispBeamColumnInt 3 4 5 $np 2 1 $C 
#element dispBeamColumn 1 1 2 $np $secTag $transfTag;


# Set axial load 
pattern Plain 1 Constant {
	load 2 -85 0.0  0.0
}

initialize;
integrator LoadControl 0 1 0 0;
system BandGeneral;
test  EnergyIncr  +1.000000E-004    25     0     2 
numberer Plain;
constraints Plain;
algorithm ModifiedNewton -initial;
analysis Static;

# perform the gravity load analysis, 
analyze [expr 1]
