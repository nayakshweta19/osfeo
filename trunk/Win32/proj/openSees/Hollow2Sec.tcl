wipe;
model  BasicBuilder  -ndm  2  -ndf  3 
node       1  +0.000000E+000  +0.000000E+000   -mass  0 0 0
#node       2  +0.000000E+000  +3.100000E+002  -mass  0 0 0
node       3  +0.000000E+000  +6.200000E+002   -mass  0 0 0
#node       4  +0.000000E+000  +9.300000E+002  -mass  0 0 0
node       5  +0.000000E+000  +1.240000E+003   -mass  0 0 0       
fix       1     1     1     1 
uniaxialMaterial Concrete01 1 -44.2 -0.002809355735456748 -35.36 -0.04 
uniaxialMaterial Concrete02 2 -35.0 -0.002499940215086195 -7.0 -0.04 0.1 3.5 1750.0
uniaxialMaterial Steel02 3 360.0 199948.00558280002 0.05 18.0 0.925 0.15 0.0 1.0 0.0 1.0 0.0
uniaxialMaterial Hysteretic 4 248.2113172752 0.0012413793103448274 322.67471245776 0.014 335.08527832152004 1.0 -248.2113172752 -0.0012413793103448274 -322.67471245776 -0.014 -335.08527832152004 -1.0 1.0 1.0 0.0 0.0 0.5
uniaxialMaterial ReinforcingSteel 5 300.0 405.0 200000.0 20000.0 0.008 0.016  -GABuck 5.0 2.0 0.0 0.0  -CMFatigue 0.26 0.506 0.389  -IsoHard 4.3 0.2  -MPCurveParams 0.333 18.0 4.0
uniaxialMaterial Steel01 6 325.0 199948.00558280002 0.05 0.0 1.0 0.0 1.0 
section Fiber 1 -GJ 1166694566960.3157 { 
layer straight 3 5 314.2 -225.0 -274.6 225.0 -274.6
layer straight 3 0 645.16 -374.6 -205.4 -250.4 -205.4
layer straight 3 2 314.2 250.4 274.6 250.4 205.4
patch quad 1 8 1 225.0 -180.0 -225.0 -180.0 -250.4 -205.4 250.4 -205.4
patch quad 2 6 12 -250.4 205.4 -374.6 205.4 -374.6 -205.4 -250.4 -205.4
layer straight 3 2 314.2 -374.6 274.6 -374.6 205.4
layer straight 3 2 314.2 -250.4 274.6 -250.4 205.4
patch quad 2 6 12 374.6 205.4 250.4 205.4 250.4 -205.4 374.6 -205.4
patch quad 2 16 8 250.4 274.6 -250.4 274.6 -250.4 205.4 250.4 205.4
patch quad 1 16 1 400.0 300.0 -400.0 300.0 -374.6 274.6 374.6 274.6
patch quad 1 24 1 400.0 300.0 374.6 274.6 374.6 -274.6 400.0 -300.0
layer straight 3 4 314.2 -250.4 180.0 -250.4 -180.0
layer straight 3 0 645.16 -250.4 274.6 -374.6 274.6
patch quad 1 24 1 -374.6 274.6 -400.0 300.0 -400.0 -300.0 -374.6 -274.6
layer straight 3 0 314.2 374.6 274.6 250.4 274.6
patch quad 1 1 12 250.4 205.4 225.0 180.0 225.0 -180.0 250.4 -205.4
layer straight 3 2 314.2 250.4 -205.4 250.4 -274.6
layer straight 3 2 314.2 -374.6 -205.4 -374.6 -274.6
layer straight 3 2 314.2 374.6 274.6 374.6 205.4
patch quad 1 1 12 -225.0 180.0 -250.4 205.4 -250.4 -205.4 -225.0 -180.0
patch quad 2 16 12 -250.4 274.6 -374.6 274.6 -374.6 205.4 -250.4 205.4
layer straight 3 2 314.2 -250.4 -205.4 -250.4 -274.6
layer straight 3 4 314.2 374.6 180.0 374.6 -180.0
patch quad 2 16 8 250.4 -205.4 -250.4 -205.4 -250.4 -274.6 250.4 -274.6
layer straight 3 5 314.2 -225.0 -205.4 225.0 -205.4
layer straight 3 4 314.2 250.4 180.0 250.4 -180.0
patch quad 2 6 12 374.6 274.6 250.4 274.6 250.4 205.4 374.6 205.4
layer straight 3 4 314.2 -374.6 180.0 -374.6 -180.0
layer straight 3 0 645.16 -250.4 205.4 -374.6 205.4
layer straight 3 0 645.16 250.4 -274.6 374.6 -274.6
layer straight 3 2 314.2 374.6 -205.4 374.6 -274.6
patch quad 1 8 1 250.4 205.4 -250.4 205.4 -225.0 180.0 225.0 180.0
patch quad 2 16 12 -250.4 -274.6 -250.4 -205.4 -374.6 -205.4 -374.6 -274.6
patch quad 1 16 1 374.6 -274.6 -374.6 -274.6 -400.0 -300.0 400.0 -300.0
layer straight 3 0 645.16 250.4 -205.4 374.6 -205.4
layer straight 3 5 314.2 225.0 274.6 -225.0 274.6
layer straight 3 0 645.16 374.6 205.4 250.4 205.4
layer straight 3 0 314.2 -374.6 -274.6 -250.4 -274.6
layer straight 3 5 314.2 225.0 205.4 -225.0 205.4
patch quad 2 16 12 374.6 -205.4 250.4 -205.4 250.4 -274.6 374.6 -274.6
};
section Fiber 2 -GJ 1085357445780.4637 { 
layer straight 4 5 50.3 -129.6 -155.0 129.6 -155.0
layer straight 4 0 50.3 -225.0 -85.0 -155.0 -85.0
layer straight 4 2 50.3 155.0 155.0 155.0 85.0
patch quad 1 8 1 130.0 -60.0 -130.0 -60.0 -155.0 -85.0 155.0 -85.0
patch quad 2 6 12 -155.0 85.0 -225.0 85.0 -225.0 -85.0 -155.0 -85.0
layer straight 4 2 50.3 -225.0 155.0 -225.0 85.0
layer straight 4 2 50.3 -155.0 155.0 -155.0 85.0
patch quad 2 6 12 225.0 85.0 155.0 85.0 155.0 -85.0 225.0 -85.0
patch quad 2 16 8 155.0 155.0 -155.0 155.0 -155.0 85.0 155.0 85.0
patch quad 1 16 1 250.0 180.0 -250.0 180.0 -225.0 155.0 225.0 155.0
patch quad 1 24 1 250.0 180.0 225.0 155.0 225.0 -155.0 250.0 -180.0
layer straight 4 1 50.3 -155.0 -59.6 -155.0 59.6
layer straight 4 0 50.3 -155.0 155.0 -225.0 155.0
patch quad 1 24 1 -225.0 155.0 -250.0 180.0 -250.0 -180.0 -225.0 -155.0
layer straight 4 0 50.3 225.0 155.0 155.0 155.0
patch quad 1 1 12 155.0 85.0 130.0 60.0 130.0 -60.0 155.0 -85.0
layer straight 4 2 50.3 155.0 -85.0 155.0 -155.0
layer straight 4 2 50.3 -225.0 -85.0 -225.0 -155.0
layer straight 4 2 50.3 225.0 155.0 225.0 85.0
patch quad 1 1 12 -130.0 60.0 -155.0 85.0 -155.0 -85.0 -130.0 -60.0
patch quad 2 16 12 -155.0 155.0 -225.0 155.0 -225.0 85.0 -155.0 85.0
layer straight 4 2 50.3 -155.0 -85.0 -155.0 -155.0
layer straight 4 1 50.3 225.0 -59.6 225.0 59.6
patch quad 2 16 8 155.0 -85.0 -155.0 -85.0 -155.0 -155.0 155.0 -155.0
layer straight 4 5 50.3 -129.6 -85.0 129.6 -85.0
layer straight 4 1 50.3 155.0 -59.6 155.0 59.6
patch quad 2 6 12 225.0 155.0 155.0 155.0 155.0 85.0 225.0 85.0
layer straight 4 1 50.3 -225.0 -59.6 -225.0 59.6
layer straight 4 0 50.3 -155.0 85.0 -225.0 85.0
layer straight 4 0 50.3 155.0 -155.0 225.0 -155.0
layer straight 4 2 50.3 225.0 -85.0 225.0 -155.0
patch quad 1 8 1 155.0 85.0 -155.0 85.0 -130.0 60.0 130.0 60.0
patch quad 2 16 12 -155.0 -155.0 -155.0 -85.0 -225.0 -85.0 -225.0 -155.0
patch quad 1 16 1 225.0 -155.0 -225.0 -155.0 -250.0 -180.0 250.0 -180.0
layer straight 4 0 50.3 155.0 -85.0 225.0 -85.0
layer straight 4 5 50.3 129.6 155.0 -129.6 155.0
layer straight 4 0 50.3 225.0 85.0 155.0 85.0
layer straight 4 0 50.3 -225.0 -155.0 -155.0 -155.0
layer straight 4 5 50.3 129.6 85.0 -129.6 85.0
patch quad 2 16 12 225.0 -85.0 155.0 -85.0 155.0 -155.0 225.0 -155.0
};
section Fiber 3 -GJ 1166694566960.3157 { 
patch quad 1 4 4 -400.0 25.0 -394.92 19.92 394.92 19.92 400.0 25.0
patch quad 1 4 4 400.0 25.0 394.92 19.92 394.92 -19.92 400.0 -25.0
patch quad 1 4 4 -400.0 25.0 -400.0 -25.0 -394.92 -19.92 -394.92 19.92
layer straight 3 1 201.1 394.92 19.92 394.92 -19.92
layer straight 3 1 645.16 5.079999999999984 19.92 -5.079999999999984 19.92
patch quad 1 8 8 -400.0 25.0 -310.0 25.0 -310.0 300.0 -400.0 300.0
layer straight 3 1 1256.6 -394.92 19.92 -394.92 -19.92
layer straight 3 1 645.16 -394.92 294.92 -394.92 30.08
layer straight 3 1 645.16 -315.08 294.92 -315.08 30.08
patch quad 2 16 16 -394.92 19.92 -394.92 -19.92 394.92 -19.92 394.92 19.92
patch quad 1 8 8 -400.0 -300.0 -310.0 -300.0 -310.0 -25.0 -400.0 -25.0
layer straight 3 1 314.2 389.84000000000003 19.92 389.84000000000003 -19.92
layer straight 3 1 645.16 -394.92 -30.08 -394.92 -294.92
patch quad 1 4 4 -394.92 -19.92 -400.0 -25.0 400.0 -25.0 394.92 -19.92
layer straight 3 1 645.16 -315.08 -30.08 -315.08 -294.92
layer straight 3 0 645.16 -389.84000000000003 19.92 -389.84000000000003 -19.92
layer straight 3 1 645.16 5.079999999999984 -19.92 -5.079999999999984 -19.92
};
section Fiber 4 -GJ 8331166899283.335 { 
patch quad 4 4 8 112.52199999999998 129.54 112.52199999999998 -129.54 134.61999999999998 -129.54 134.61999999999998 129.54
patch quad 4 8 4 -134.61999999999998 129.54 -134.61999999999998 -129.54 -112.52199999999998 -129.54 -112.52199999999998 129.54
patch quad 4 4 8 -112.52199999999998 -112.52199999999998 6.731 -112.52199999999998 -6.731 112.52199999999998 -6.731 112.52199999999998 6.731
};
section Fiber 5 -GJ 1166694566960.3157 { 
layer straight 6 5 50.3 -180.0 -150.0 180.0 -150.0
layer straight 6 0 50.3 -270.0 -90.0 -210.0 -90.0
layer straight 6 2 50.3 210.0 150.0 210.0 90.0
patch quad 1 12 1 180.0 -60.0 -180.0 -60.0 -210.0 -90.0 210.0 -90.0
patch quad 2 6 4 -210.0 90.0 -270.0 90.0 -270.0 -90.0 -210.0 -90.0
layer straight 6 2 50.3 -270.0 150.0 -270.0 90.0
layer straight 6 2 50.3 -210.0 150.0 -210.0 90.0
patch quad 2 6 4 270.0 90.0 210.0 90.0 210.0 -90.0 270.0 -90.0
patch quad 2 4 12 210.0 150.0 -210.0 150.0 -210.0 90.0 210.0 90.0
patch quad 1 18 1 300.0 180.0 -300.0 180.0 -270.0 150.0 270.0 150.0
patch quad 1 12 1 300.0 180.0 270.0 150.0 270.0 -150.0 300.0 -180.0
layer straight 6 1 50.3 -210.0 60.0 -210.0 -60.0
patch quad 1 12 1 -270.0 150.0 -300.0 180.0 -300.0 -180.0 -270.0 -150.0
patch quad 1 1 6 210.0 90.0 180.0 60.0 180.0 -60.0 210.0 -90.0
layer straight 6 2 50.3 210.0 -90.0 210.0 -150.0
layer straight 6 2 50.3 -270.0 -90.0 -270.0 -150.0
layer straight 6 2 50.3 270.0 150.0 270.0 90.0
patch quad 1 1 6 -180.0 60.0 -210.0 90.0 -210.0 -90.0 -180.0 -60.0
patch quad 2 4 4 -210.0 150.0 -270.0 150.0 -270.0 90.0 -210.0 90.0
layer straight 6 2 50.3 -210.0 -90.0 -210.0 -150.0
layer straight 6 1 50.3 270.0 60.0 270.0 -60.0
patch quad 2 4 12 210.0 -90.0 -210.0 -90.0 -210.0 -150.0 210.0 -150.0
layer straight 6 5 50.3 -180.0 -90.0 180.0 -90.0
layer straight 6 1 50.3 210.0 60.0 210.0 -60.0
patch quad 2 6 4 270.0 150.0 210.0 150.0 210.0 90.0 270.0 90.0
layer straight 6 1 50.3 -270.0 60.0 -270.0 -60.0
layer straight 6 2 50.3 270.0 -90.0 270.0 -150.0
patch quad 1 12 1 210.0 90.0 -210.0 90.0 -180.0 60.0 180.0 60.0
patch quad 2 4 4 -210.0 -150.0 -210.0 -90.0 -270.0 -90.0 -270.0 -150.0
patch quad 1 18 1 270.0 -150.0 -270.0 -150.0 -300.0 -180.0 300.0 -180.0
layer straight 6 5 50.3 180.0 150.0 -180.0 150.0
layer straight 6 0 50.3 270.0 90.0 210.0 90.0
layer straight 6 0 645.16 -270.0 -150.0 -210.0 -150.0
layer straight 6 5 50.3 180.0 90.0 -180.0 90.0
patch quad 2 4 4 270.0 -90.0 210.0 -90.0 210.0 -150.0 270.0 -150.0
};
section Fiber 6 -GJ 1166694566960.3157 { 
patch circ 1 36 24 0.0 0.0 1470 1500 0 360
layer circ 5 16 1256.6 0.0 0.0 830 0 360
patch circ 2 36 24 0.0 0.0 830 1470 0 360
};
geomTransf  PDelta       1;
element  dispBeamColumn       1       1       3     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
element  dispBeamColumn       2       3       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element  dispBeamColumn       3       3       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 
#element  dispBeamColumn       4       4       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 

timeSeries  Linear       1  -factor  +1.000000E+000 
initialize 
pattern  Plain       1  "Linear" { 
    load       5  +0.000000E+000  -2.800000E+005  +0.000000E+000 
} 
recorder  Node  -file  Node_DefoShape_Dsp.out  -time -nodeRange 1  5 -dof  1  2  3  disp 

constraints  Plain ;
test  EnergyIncr  +1.000000E-004    25     0     2 ;
integrator  LoadControl  +1.000000E+000 ;
algorithm  KrylovNewton ;
numberer  Plain ;
system  BandGeneral ;
analysis  Static ;
analyze     1 ;

loadConst -time 0.0 
remove recorders 
wipeAnalysis 
pattern  Plain       2  "Linear" { 
    load       5  1.000000E+005  +0.000000E+000   +0.000000E+000 
} 
#recorder Element -file hollow1_Col_1_se.dat -time -ele 1 section fiber 0. 150.0 stressStrain
#recorder Element -file hollow1_ele_p1_S.dat -time -ele 1 2 3 4 integrPoint 1 stresses
#recorder Element -file hollow1_ele_p1_E.dat -time -ele 1 2 3 4 integrPoint 1 strains
recorder Node -file hollow2_nodeReaction.dat -time -node 1 -dof 1 2 3 reaction;
recorder Node -file hollow2_node_disp_Output.dat -time -node 5 -dof 1 2 disp;


constraints Plain ; #Lagrange
set numbererTypeStatic Plain; numberer $numbererTypeStatic ;
set systemTypeStatic BandGeneral; system $systemTypeStatic ;
set TolStatic 1.0e-3; set Tol [expr $TolStatic]; set maxNumIterStatic 1000; set printFlagStatic 5;   
set testTypeStatic NormDispIncr ;	test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;

variable maxNumIterConvergeStatic 2000;	
variable printFlagConvergeStatic 0;
set Dincr [expr 1.0];	# displacement increment. you want this to be small, but not too small to slow analysis
set IDctrlNode 5; set IDctrlDOF 1;
set  algorithmTypeStatic KrylovNewton ; algorithm $algorithmTypeStatic;        
integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr;
set analysisTypeStatic Static; analysis $analysisTypeStatic ;
initialize
set iDmax { 5 10 15 20 25 30 35 40 45 50 55 60 65 70 } ; ##20 25 30 35 40 45 50 60 70 80 100 120 140 160 240 };  # vector of displacement-cycle peaks, in terms of storey drift ratio
set CycleType Full;  #  Full / Push / Half
set Ncycles 2;
source GeneratePeaks.tcl

foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep            
			puts "D1=$D1"
			set Dincr [expr $D1 - $D0]
			#puts "Dincr = $Dincr"
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			puts "Displacement Node: $IDctrlNode. Dof: $IDctrlDOF. Incr: $Dincr"
			#integrator ArcLength 1.0 0.1
			# ------------------first analyze command------------------------
			set ok [analyze 1]
			# -----------------if convergence failure-------------------------
			# max no. of iterations performed before "failure to converge" is ret'd
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent .."
				test NormDispIncr $Tol 2000 0;
				algorithm Newton -initial;
				set ok [analyze 1];
				test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;
				algorithm $algorithmTypeStatic;
				puts "Trying Newton with Initial Tangent failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Broyden ..";
				algorithm Broyden 40
				set ok [analyze 1 ]
				algorithm $algorithmTypeStatic
				puts "Trying Broyden failed to converge...";
			}
			if {$ok != 0} {
				puts "Trying Newton With LineSearch ..";
				algorithm NewtonLineSearch 0.8 
				set ok [analyze 1]
				algorithm $algorithmTypeStatic
				puts "Trying Newton With LineSearch failed to converge...";
			}
			#if {$ok != 0} {
			#	puts "Trying ArcLength ..";
			#	algorithm $algorithmTypeStatic 
			#	integrator ArcLength 1.0 0.1
			#	set ok [analyze 1]
			#	puts "Trying ArcLength With KrylovNewton failed to converge...";
			#}
			if {$ok != 0} {
				set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";
				set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] m]
				puts $putout
				return -1
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
			puts "$D0";
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
remove recorders;