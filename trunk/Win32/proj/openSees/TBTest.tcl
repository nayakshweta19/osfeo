wipe;
model basic -ndm 2 -ndf 3 
set tolAx 1.0e-3;
set iterAx 100;
set tolLatNew 1.0e-6;
set iterLatNew 100;
set tolLatIni 1.0e-5;
set iterLatIni 1000; 
set dUi 0.1; # Displacement increment
set maxU 100; # Max. Displacement

#concrete parameters
set fc -25.0;
set eo -0.003;
set r 1.88; 
set k 1.0;
set alphaC 0.32;
set fcr 2.14;
set ecr 0.00008;
set b 0.4;
set alphaT 0.08;
set fcc -18.9;
set ecc -0.0038;
set rcc 1.88; 

#steel parameters
set fyRed 345.0;
set by 0.02;
set Ry 1.564;
set roy 0.00327;
set fyBRed 325.0;
set byB 0.02; 
set RyB 25.0;
set royB 0.0293;
set fxRed 325.0;
set bx 0.02;
set Rx 1.564;
set rox 0.00327;

#concrete (confined and unconfined)
uniaxialMaterial Concrete06 1 $fcc $ecc $rcc $k $alphaC $fcr $ecr $b $alphaT;
uniaxialMaterial Concrete06 2 $fc $eo $r $k $alphaC $fcr $ecr $b $alphaT;
#uniaxialMaterial Concrete01 1 $fcc  $ecc [expr 0.2*$fcc] -0.009
#uniaxialMaterial Concrete01 2 $fc $eo [expr 0.2*$fcc] -0.009

# steel (bound., web and horiz. reinforcement)
set E 2.0e5; 
uniaxialMaterial Steel02 1003 $fyBRed $E $byB $RyB 0.9 0.1 0 0.1 0 0.1 ; 
uniaxialMaterial Steel02 1004 $fyRed $E $by $Ry 0.9 0.1 0 0.1 0 0.1 ; 
uniaxialMaterial Steel02 1005 $fxRed $E $bx $Rx 0.9 0.1 0 0.1 0 0.1 ;

# Define cross-section 
set t1 360;
set NStrip1 1; # thickness 1
set t2 120;
set NStrip2 4; # thickness 2
set t3 360;
set NStrip3 1; # thickness 3
geomTransf Timoshenko 1
set np 4; # int. points

#section definition 
section CSMMFiber2d 2 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { 
#vertical fibers: fiber yLoc zLoc area matTag
fiber -240    0 43200 1;   fiber -240    0 503   1003;
fiber  -90  135 21600 1;   fiber  -90  135 251.5 1003;
fiber   90  135 21600 1;   fiber   90  135 251.5 1003;
fiber  -90 -135 21600 1;   fiber  -90 -135 251.5 1003;
fiber   90 -135 21600 1;   fiber   90 -135 251.5 1003;
fiber  240    0 43200 1;   fiber  240    0 503   1003;
#horiz. reinf.: Hfiber yLoc zLoc area matTag
Hfiber 0 0 1680 1005;
}

#nodes
node 1 0 0.0;
node 2 0 1000.0;

mass 1 1e6 1e6 0;
mass 2 1e6 1e6 0;

#element definition
element Timoshenko 1 1 2 $np 2 1;

fix 1 1 1 1;
fix 2 0 0 0;

# Set axial load 
pattern Plain 1 Constant {
	load 2 0.0 -2000.0 0.0
}

initialize;
integrator LoadControl 0 1 0 0;
system BandGeneral;
test NormUnbalance $tolAx $iterAx 0;
numberer Plain;
constraints Plain;
algorithm ModifiedNewton -initial;
analysis Static;


# perform the gravity load analysis, 
analyze [expr 1]

## Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0
remove recorders

# set lateral load 
pattern Plain 2 Linear { 
	load 2 1.0 0.0 1.0;
} 
system BandGeneral;
constraints Plain;
numberer Plain;

# Create a recorder to monitor nodal displacement and element forces
recorder Node -file TBnodeTop.out -time -node 2 -dof 1 2 disp

# recorder for element1 section1 steel stress/strain and section force-def.
#recorder Element -file Sect_FandD.out -ele 1 section 1 forceAndDeformation
#recorder Element -file Sect_eX.out -ele 1 section 1 eX
#recorder Element -file Sect_eY.out -ele 1 section 1 eY
#recorder Element -file Sect_sX.out -ele 1 section 1 sX
#recorder Element -file Sect_sY.out -ele 1 section 1 sY
set dU $dUi;
set numSteps [expr int($maxU/$dU)];
integrator DisplacementControl 2 1 $dU 1 $dU $dU
test NormDispIncr $tolLatNew $iterLatNew; 
algorithm KrylovNewton
analysis Static


#source cycle2.tcl
#### end

set ok [analyze $numSteps];
set jump 1;
if {$ok != 0} {
	set currentDisp [nodeDisp 2 1]
	set ok 0
	while {abs($currentDisp) < abs($maxU)} {
		set ok [analyze 1]
		puts "\n Trying.. $currentDisp\n"
		# if the analysis fails try initial tangent iteration
		if {$ok != 0} {
			puts "\n regular newton failed .. try an initial stiffness"
			test NormDispIncr $tolLatIni $iterLatIni 0; 
			algorithm ModifiedNewton -initial ; 
			set ok [analyze 1]
			puts "\n Trying.. $currentDisp\n"
			if {$ok == 0} {
				puts " that worked .. back to regular newton \n"
				set jump 1
				integrator DisplacementControl 2 1 $dU 1 
			} else {
				puts "\n that didn't worked .. Try next point\n"
				set jump [expr $jump+1];
				integrator DisplacementControl 2 1 [expr $dU*$jump] 1 
			}
			test NormDispIncr $tolLatNew $iterLatNew 
			algorithm Newton
			puts "[nodeDisp 2 1]\t[nodeDisp 2 2]\t[nodeDisp 2 3]"
		} else {
			set jump 1
			integrator DisplacementControl 2 1 $dU 1 
		}
		set currentDisp [expr $currentDisp+$dU]
	}
}
#### end
