wipe;
model basic -ndm 2 -ndf 3 
set tolAx 1.0e-3;
set iterAx 100;
set tolLatNew 1.0e-6;
set iterLatNew 100;
set tolLatIni 1.0e-5;
set iterLatIni 1000; 
set dUi 0.02; # Displacement increment
set maxU 1.6; # Max. Displacement

#concrete parameters
set fc -6.21;
set eo -0.0021;
set r 2.91; 
set k 1.36;
set alphaC 0.32;
set fcr 0.294;
set ecr 0.00008;
set b 0.4;
set alphaT 0.08;
set fcc -6.91;
set ecc -0.0038;
set rcc 1.88; 

#steel parameters
set fyRed 59.15;
set by 0.02;
set Ry 1.564;
set roy 0.00327;
set fyBRed 57.33;
set byB 0.02; 
set RyB 25.0;
set royB 0.0293;
set fxRed 59.15;
set bx 0.02;
set Rx 1.564;
set rox 0.00327;

#concrete (confined and unconfined)
#uniaxialMaterial Concrete06 1 $fcc $ecc $rcc 1 $alphaC $fcr $ecr $b $alphaT;
#uniaxialMaterial Concrete06 2 $fc $eo $r $k $alphaC $fcr $ecr $b $alphaT;
#uniaxialMaterial Concrete01 1 -72.4 -0.0024 -14.0 -0.008
#uniaxialMaterial Concrete01 2 -72.4 -0.0024 -14.0 -0.008

# steel (bound., web and horiz. reinforcement)
set E 29000.0; 
#uniaxialMaterial Steel02 1003 $fyBRed $E $byB $RyB 0.9 0.1 0 0.1 0 0.1 ; 
#uniaxialMaterial Steel02 1004 $fyRed $E $by $Ry 0.9 0.1 0 0.1 0 0.1 ; 
#uniaxialMaterial Steel02 1005 $fxRed $E $bx $Rx 0.9 0.1 0 0.1 0 0.1 ;

# set fc fy E
set wfc 28.0;
set wfyv 560;
set wfyh1 450;
set wE 190000.0;
set rou1 0;
set rou2 0;
set rouv 0.0246;
set rouh1 0.00246; # #10@70
# UniaxialMaterial: steelZ01
# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteL02 14 [expr -$wfc] -0.003
uniaxialMaterial ConcreteL02 15 [expr -$wfc] -0.003
set pi 3.141592654
# NDMaterial: FAFourSteelPCPlaneStress
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelPCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteetPCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteetPCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAFourSteelRCPlaneStress  (100)
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 $rouv $rouh1 $wfc 0.0 $wfyv $wE 0.003

#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteelRCPlaneStress 21 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
nDMaterial CSMMRCPlaneStress 21 0.0 11 12 14 15 [expr 0.5*$pi] [expr 0.0*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002


# Define cross-section 
set t1 6.0;
set NStrip1 2; # thickness 1
set t2 4.0;
set NStrip2 4; # thickness 2
set t3 6.0;
set NStrip3 2; # thickness 3
geomTransf LinearInt 1
set np 1; # int. points
set C 0.4; # center of rotation 

#section definition 
section CSMMFiber2d 2 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { 
fiber2d -23 0 10 21; fiber -23 0 5 21;
fiber2d -19 0 10 21; fiber -19 0 5 21; 
fiber2d -14 0 33 21;
fiber2d -6 0 33 21;
fiber2d 6 0 33 21;
fiber2d 14 0 33 21;
fiber2d 19 0 10 21; fiber 19 0 5 21;
fiber2d 23 0 10 21; fiber 23 0 5 21;
}

#nodes
set n 2
node 1 0 0.0;
node 2 0 20;

mass 1 1e6 1e6 0
mass 2 1e6 1e6 0

#element definition
element dispBeamColumnInt 1 1 2 $np 2 1 $C 
#element dispBeamColumnInt 2 3 4 $np 2 1 $C
#element dispBeamColumnInt 3 4 5 $np 2 1 $C 
#element dispBeamColumnInt 4 5 6 $np 2 1 $C
#element dispBeamColumnInt 5 6 7 $np 2 1 $C
#element dispBeamColumnInt 6 7 8 $np 2 1 $C
#element dispBeamColumnInt 7 8 9 $np 2 1 $C
#element dispBeamColumnInt 8 9 2 $np 2 1 $C

fix 1 1 1 1
fix 2 0 0 0

# Set axial load 
pattern Plain 1 Constant {
	load 2 -85 0.0  0.0
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
	load 2 1.0 0.0 0.0;
} 
system BandGeneral;
constraints Transformation;
numberer Plain;

# Create a recorder to monitor nodal displacement and element forces
recorder Node -file nodeTop.out -node 2 -dof 1 disp
recorder Element -file elesX.out -time -ele 1 globalForce

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
algorithm ModifiedNewton -initial 
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
