wipe; # Based on hollow1.tcl
# Changes: 1. transfer into reinforcedSteel Materal model
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  2  -ndf  3 

# Define geometry 
# --------------- 

# Node    tag    xCrd    yCrd 
node       1  +0.000000E+000  +0.000000E+000 
#node       2  +0.000000E+000  +3.100000E+002 
#node       3  +0.000000E+000  +1.700000E+003 
#node       4  +0.000000E+000  +9.300000E+002 
node       5  +0.000000E+000  +3.40000E+003 

# Define Single Point Constraints 
# ------------------------------- 
# SPConstraint.tcl 

# SPC    tag    Dx    Dy    Rz 
fix       1     1     1     1 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.900000E+004  +0.000000E+000 

# Material "CoreConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
# C30
set Fc 14.3
set rFc 1.1;
set rFu 0.2;
set E [expr 4733.*sqrt($rFc*abs($Fc))]; #for psi: [expr 57000.*sqrt($rFc*abs($Fc))];
set fpc [expr -$rFc*abs($Fc)];
set epsc0 [expr 2.*$fpc/$E];
set fpcU [expr -$rFu*$rFc*abs($Fc)];
set epsU -0.02;
set lambda 0.1;
set fT  +0.1*$rFc*abs($Fc);
set Ets  $fT/0.002;

uniaxialMaterial Concrete02 2  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)] $lambda [expr abs(1.*$fT)] [expr abs(1.*$Ets)];
#uniaxialMaterial  Concrete02         -3.900000E+001  -2.000000E-003  -4.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

#uniaxialMaterial ConfinedConcrete01 $tag  $secType $fpc   $Ec -epscu $epscu  $nu     $L1  $phis  $S   $fyh   $Es0   $haRatio  $mu   $phiLon  -stRatio  $stRatio
#uniaxialMaterial ConfinedConcrete01    2     S1    $fpc  $E  -epscu -0.03    -nu 0.2 100.0 6.0 40.0 300.0 206000.0    0.00   200.0   8.0   -stRatio    0.85


# Material "CoverConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
set Fc 13.1;
set rFc 1.0;
set rFu 0.2;
set E [expr 4733.*sqrt($rFc*abs($Fc))]; #for psi: [expr 57000.*sqrt($rFc*abs($Fc))];
set fpc [expr -$rFc*abs($Fc)];
set epsc0 [expr 2.*$fpc/$E]
set fpcU [expr -$rFu*$rFc*abs($Fc)]
set epsU -0.02;
set lambda 0.1;
set fT  +0.1*$rFc*abs($Fc);
set Ets  $fT/0.002;
#uniaxialMaterial Concrete02 3  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)] $lambda [expr abs(1.*$fT)] [expr abs(1.*$Ets)]
#uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 
uniaxialMaterial  Concrete01 3  [expr -abs(1.*$fpc)] [expr -abs(1.*$epsc0)] [expr -abs(1.*$fpcU)] [expr -abs(1.*$epsU)]

# Material "reinforcedBar":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
#uniaxialMaterial  Steel02       4  +2.100000E+002  +2.000000E+005  +1.000000E-001  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "reinforcedSteel360":    matTag    Fy    Fu    E    Eh    epsh    epsu    <-GABuck  lsr  beta  r  gamma>    <-DMBuck  lsr  <alpha>>    <-CMFatigue  Cf  alpha  Cd>    <-IsoHard  <a1  <limit>>>    <-MPCurveParams  R1  R2  R3> 
#uniaxialMaterial  ReinforcingSteel       4  +3.600000E+002  +3.350000E+002  +2.000000E+005  +5.000000E+002  +5.000000E-002  +1.200000E-001  -GABuck  +6.000000E+000  +1.000000E+000  +4.000000E-001  +5.000000E-001  -DMBuck  +6.000000E+000  +1.000000E+000  -CMFatigue  +2.600000E-001  +5.060000E-001  +3.890000E-001  -IsoHard  +4.300000E+000  +1.000000E-002  -MPCurveParams  +3.330000E-001  +1.800000E+001  +4.000000E+000 
set MaterialTag 4;
set fY 200.;
set fU [expr 1.3*$fY];
set Es 1.95e5;
set Esh [expr 0.1*$Es];
set esh 0.002;
set eult 0.004;
set lsr1 6.0;
set beta 2.0;
set r 0.0; #0.4;
set gama 0.5;
set lsr2 8.0
set alpha1 0.95
set Cf 0.2;
set alpha2 0.5;
set Cd 0.6; #0.389;
set a1 0.3;
set limit 0.01; #0.01;
set R1 0.3;
set R2 5;
set R3 1;

uniaxialMaterial ReinforcingSteel $MaterialTag [expr 1.*$fY] [expr 1.*$fU] [expr 1.*$Es] [expr 1.*$Esh] [expr 1.*$esh] [expr 1.*$eult] \
                                               -GABuck [expr 1.*$lsr1] [expr 1.*$beta] [expr 1.*$r] [expr 1.*$gama] \
                                               -CMFatigue [expr 1.*$Cf] [expr 1.*$alpha2] [expr 1.*$Cd] \
                                               -IsoHard [expr 1.*$a1] [expr 1.*$limit] \
                                               -MPCurveParams [expr 1.*$R1] [expr 1.*$R2] [expr 1.*$R3]
#                                               -DMBuck [expr 1.*$lsr2] [expr 1.*$alpha1] \
#                                               -GABuck [expr 1.*$lsr1] [expr 1.*$beta] [expr 1.*$r] [expr 1.*$gama] \

# Define section(s) 
# ----------------- 
# Sections.tcl 

# Section "ElasticDefault":    secTag    E    A    Iz 
section  Elastic       1  +2.900000E+004  +1.800000E+002  +4.860000E+003 

# Section "Hollow01":    secTag 
section  Fiber       2  { 
    # PatchBox "Patch01":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3    12     1  +2.700000E+002  +1.800000E+002  +2.700000E+002  -1.800000E+002  +3.000000E+002  -1.800000E+002  +3.000000E+002  +1.800000E+002 
    patch  quad       3     1    18  -2.700000E+002  +1.800000E+002  -2.700000E+002  +1.500000E+002  +2.700000E+002  +1.500000E+002  +2.700000E+002  +1.800000E+002 
    patch  quad       3     1    18  -2.700000E+002  -1.500000E+002  -2.700000E+002  -1.800000E+002  +2.700000E+002  -1.800000E+002  +2.700000E+002  -1.500000E+002 
    patch  quad       3    12     1  -3.000000E+002  +1.800000E+002  -3.000000E+002  -1.800000E+002  -2.700000E+002  -1.800000E+002  -2.700000E+002  +1.800000E+002 
    # PatchBox "Patch02":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       2    10     2  +2.100000E+002  +1.500000E+002  +2.100000E+002  -1.500000E+002  +2.700000E+002  -1.500000E+002  +2.700000E+002  +1.500000E+002 
    patch  quad       2     2    14  -2.100000E+002  +1.500000E+002  -2.100000E+002  +9.000000E+001  +2.100000E+002  +9.000000E+001  +2.100000E+002  +1.500000E+002 
    patch  quad       2     2    14  -2.100000E+002  -9.000000E+001  -2.100000E+002  -1.500000E+002  +2.100000E+002  -1.500000E+002  +2.100000E+002  -9.000000E+001 
    patch  quad       2    10     2  -2.700000E+002  +1.500000E+002  -2.700000E+002  -1.500000E+002  -2.100000E+002  -1.500000E+002  -2.100000E+002  +1.500000E+002 
    # PatchBox "Patch03":    matTag    NSIJ    NSJK    Iy    Iz    Jy    Jz    Ky    Kz    Ly    Lz 
    patch  quad       3     6     1  +1.800000E+002  +9.000000E+001  +1.800000E+002  -9.000000E+001  +2.100000E+002  -9.000000E+001  +2.100000E+002  +9.000000E+001 
    patch  quad       3     2    12  -1.800000E+002  +9.000000E+001  -1.800000E+002  +6.000000E+001  +1.800000E+002  +6.000000E+001  +1.800000E+002  +9.000000E+001 
    patch  quad       3     2    12  -1.800000E+002  -6.000000E+001  -1.800000E+002  -9.000000E+001  +1.800000E+002  -9.000000E+001  +1.800000E+002  -6.000000E+001 
    patch  quad       3     6     1  -2.100000E+002  +9.000000E+001  -2.100000E+002  -9.000000E+001  -1.800000E+002  -9.000000E+001  -1.800000E+002  +9.000000E+001 
    # LayerStraight "Layer01":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +1.500000E+002  -2.700000E+002  +1.500000E+002 
    # LayerStraight "Layer02":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  +9.000000E+001  -2.700000E+002  +9.000000E+001 
    # LayerStraight "Layer03":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  -9.000000E+001  -2.700000E+002  -9.000000E+001 
    # LayerStraight "Layer04":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     8  +5.030000E+001  +2.700000E+002  -1.500000E+002  -2.700000E+002  -1.500000E+002 
    # LayerStraight "Layer05":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  +2.700000E+002  +3.000000E+001  +2.700000E+002  -3.000000E+001 
    # LayerStraight "Layer06":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  +2.100000E+002  +3.000000E+001  +2.100000E+002  -3.000000E+001 
    # LayerStraight "Layer07":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  -2.100000E+002  +3.000000E+001  -2.100000E+002  -3.000000E+001 
    # LayerStraight "Layer08":    matTag    numBar    areaBar    yStart    zStart    yEnd    zEnd 
    layer  straight       4     2  +5.030000E+001  -2.700000E+002  +3.000000E+001  -2.700000E+002  -3.000000E+001 
} 

# Define geometric transformation(s) 
# ---------------------------------- 
# GeoTran.tcl Corotational PDelta Linear
geomTransf  PDelta       1 

# Define element(s) 
# ----------------- 
# Element "Beam":            eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
element  dispBeamColumn       1       1       5     4     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  dispBeamColumn       2       3       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  dispBeamColumn       3       3       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 

# Element "Beam":    eleTag    NodeI    NodeJ    NIP    secTag    geoTranTag    <-mass massDens>    <-iter maxIters tol> 
#element  dispBeamColumn       4       4       5     5     2     1  -mass +4.665600E+005  -iter   20  +1.000000E-004 


# Define time series 
# ------------------ 
# TimeSeries.tcl 

# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# Get Initial Stiffness 
# --------------------- 
initialize 
pattern  Plain       1  "Linear" { 
    load       5  +0.000000E+000  -2.800000E+005  +0.000000E+000 
} 

recorder  Node  -file  Node_DefoShape_Dsp.out  -time -nodeRange 1  5 -dof  1  2  3  disp 
constraints  Plain 
test  EnergyIncr  +1.000000E-004    25     0     2 
integrator  LoadControl  +1.000000E+000 
algorithm  KrylovNewton 
numberer  Plain 
system  BandGeneral 
analysis  Static 
analyze     1 

# Reset for next analysis case 
# ---------------------------- 
loadConst -time 0.0 
remove recorders 
wipeAnalysis 

pattern  Plain       2  "Linear" { 
    # Load    nodeTag    LoadValues 
    load       5  1.000000E+005  +0.000000E+000   +0.000000E+000 
} 

recorder Node -file hollow2_nodeReaction.dat -time -node 1 -dof 1 2 3 reaction
recorder Node -file hollow2_node_disp_Output.dat -time -node 5 -dof 1 2 disp;

set Dincr [expr 1.0];	# displacement increment. you want this to be small, but not too small to slow analysis
set IDctrlNode 5;
set IDctrlDOF 1;

set RigidDiaphragm "OFF";
variable constraintsTypeStatic Plain;  # default;
if {  [info exists RigidDiaphragm] == 1} {
	if {$RigidDiaphragm=="ON"} {
		variable constraintsTypeStatic Lagrange;	#     for large model, try Transformation
	};	# if rigid diaphragm is on
};	# if rigid diaphragm exists
constraints $constraintsTypeStatic
set numbererTypeStatic Plain
numberer $numbererTypeStatic 
set systemTypeStatic BandGeneral;		# try UmfPack for large model  
system $systemTypeStatic 
variable TolStatic 1.0e-3;                        # Convergence Test: tolerance
variable Tol [expr $TolStatic];
variable maxNumIterStatic 1000;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
variable printFlagStatic 5;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
variable testTypeStatic NormDispIncr ;	# Convergence-test type: NormDispIncr, EnergyIncr
test $testTypeStatic $TolStatic $maxNumIterStatic $printFlagStatic;
variable maxNumIterConvergeStatic 2000;	
variable printFlagConvergeStatic 0;
variable algorithmTypeStatic KrylovNewton
algorithm $algorithmTypeStatic;        
integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr
set analysisTypeStatic Static
analysis $analysisTypeStatic 

# initialize in case we need to do an initial stiffness iteration
initialize

# Perform the Cycle pushover analysis
# characteristics of cyclic analysis	
#set iDmax { 5 10 15 20 25 30 35 40 45 50 55 60 65 70 } ; ##20 25 30 35 40 45 50 60 70 80 100 120 140 160 240 };  # vector of displacement-cycle peaks, in terms of storey drift ratio
set iDmax { 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160} ; ##20 25 30 35 40 45 50 60 70 80 100 120 140 160 240 };  # vector of displacement-cycle peaks, in terms of storey drift ratio
set CycleType Full;    # you can do Full / Push / Half cycles with the proc
set Ncycles 2;				# specify the number of cycles at each peak
#  ---------------------------------    perform Static Cyclic Displacements Analysis
source GeneratePeaks.tcl

# -- STATIC PUSHOVER/CYCLIC ANALYSIS

foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep            
			#puts "D1=$D1"
			set Dincr [expr $D1 - $D0]
			#puts "Dincr = $Dincr"
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			#puts "Displacement Node: $IDctrlNode. Dof: $IDctrlDOF. Incr: $Dincr"
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
			#puts "$D0";
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
remove recorders;