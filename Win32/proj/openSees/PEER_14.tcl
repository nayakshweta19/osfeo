# Test14.tcl: Rc joint subassembly
# Units: KN, M , KPa
# Arash, Apr 2003
# email:altntash@stanford.edu
#
# Time hystory alnalysis, displacement controlled
# Centerline model, elastic
#

set InputFile  PEER.txt;

# Number of integration points
set np 7
set NumSubElem	5
set AxialLoad	-650.0

# Define the model builder
model BasicBuilder -ndm 2 -ndf 3

# Define    nodes
#    tag	X							Y
node  101	1.83							0.0
node  199	0.23							0.0
node  201	0.0							1.07
node  299	0.0							0.25
node  301	-1.83							0.0
node  399	-0.23							0.0
node  401	0.0							-1.07
node  499	0.0							-0.25

for {set i 1} {$i < $NumSubElem} {incr i 1} {
	node  [expr 101+$i]	[expr 1.83 - $i*1.6/$NumSubElem]					0.0
	node  [expr 201+$i]	0.0					[expr 1.07 - $i*0.82/$NumSubElem]
	node  [expr 301+$i]	[expr -1.83 + $i*1.6/$NumSubElem]					0.0
	node	[expr 401+$i]	0.0					[expr -1.07 + $i*0.82/$NumSubElem]
}


# Single point constraints
#    node  DX    DY    RZ
fix	201	1	0	0
fix	401	1	1	0


# Define materials
#
#                          tag -f'c(ksi) -epsco  -f'cu      -epscu ratio ft       Ets
# Column cover
uniaxialMaterial Concrete02 1 -18.3254e3 -0.002   -0.0     -0.00683 0.1	0.0       10.0
# Column core
uniaxialMaterial Concrete02 2 -45.7177e3 -0.00273 -7.144e3 -0.06127 0.1	1.27914e3  10.0 
# Beam cover
uniaxialMaterial Concrete02 3 -18.3254e3 -0.002   -3.665e3 -0.00683 0.1	0.0         10.0 
# Beam core
uniaxialMaterial Concrete02 4 -33.7831e3 -0.00258 -6.757e3 -0.05271 0.1	1.27914e3   10.0 

# Steel model
#			tag   fy(kPa)   E(kPa)  b     R     c1     c2    a1   a2   a3   a4
uniaxialMaterial Steel02  5  420.298e3  2.03e8  0.02  18.5  0.925  0.15  0.0  1.0  0.0  1.0

source RCColumnsection.tcl
RCColumnsection	1	0.4572	0.4064	0.035		2	1	5	3	3.8794791e-4	2.8502296e-4	50	1	25	1

source RCbeamsection.tcl
RCbeamsection	2	0.508		0.4064	0.035		4	3	5	4	3.8794791e-4	4	1.9793261e-4	50	1

# Geometric transformation
geomTransf Linear 1

# Define elements
# Columns
#                           tag   ndI   ndJ    nPts  secID  transf
for {set i 1} {$i < $NumSubElem} {incr i 1} {
element dispBeamColumn	[expr 100+$i]	[expr 100+$i]	[expr 101+$i]	$np	1	1
element dispBeamColumn	[expr 300+$i]	[expr 300+$i]	[expr 301+$i]	$np	1	1
element dispBeamColumn	[expr 200+$i]	[expr 200+$i]	[expr 201+$i]	$np	2	1
element dispBeamColumn	[expr 400+$i]	[expr 400+$i]	[expr 401+$i]	$np	2	1
}

element dispBeamColumn	199	[expr 100+$NumSubElem]	199	$np	1	1
element dispBeamColumn	399	[expr 300+$NumSubElem]	399	$np	1	1
element dispBeamColumn	299	[expr 200+$NumSubElem]	299	$np	2	1
element dispBeamColumn	499	[expr 400+$NumSubElem]	499	$np	2	1

# Damage models for deteriorating material parameters
#			tag?	deltaU?	beta?	sigmaY?
damageModel ParkAng	1	0.1		0.15	308.754
damageModel ParkAng	2	0.1		0.15	215.493
damageModel ParkAng	3	0.1		0.15	218.813
damageModel ParkAng	4	0.15		0.15	357.52

# Deteriorating damage models
uniaxialMaterial PinchingDamage 11  99920.39	308.754	-308.754	0.04756	1.0		-0.01		0.1	-0.1	0.25	0.25   0.25   1   1   1   1
uniaxialMaterial PinchingDamage 12  113417.368	215.493	-218.813	0.13641	1.0	 	-0.01		0.1	-0.1	0.25	0.25   0.25   2   2   2   2
uniaxialMaterial PinchingDamage 13  75452.759	218.813	-215.493	0.13888	1.0	 	-0.01		0.1	-0.1	0.25	0.25   0.25   3   3   3   3
uniaxialMaterial PinchingDamage 14	59586.667	357.52	-357.52	0.07217	0.77736	-0.04266	0.02	-0.02	0.25	0.25   0.25   4   4   4   4

# joint with hinges, with damage calculations
#                tag  nd1 nd2 nd3 nd4  ndC  Mat1 mat2 mat3 mat4 matC lrg
element Joint2D  1    199 299 399 499   1     12   11   13   11   14   2	-damage	2	1	3	1	4


# Constant gravity loads
pattern Plain 1 Linear {
	#	node	FX	FY		MZ
	load	201	0.0	$AxialLoad	0.0 
}

initialize

test NormDispIncr 1.0e-8 10 1
algorithm Newton
integrator LoadControl .1 10 .1 .1
numberer RCM
constraints Penalty   1.0e12   1.0e12
system SparseGeneral -piv
analysis Static
analyze 10

loadConst -time 0.0
wipeAnalysis

set dt 1.0;
set factor 1.0;

pattern Plain 2 "Series -dt $dt -filePath $InputFile -factor $factor" {
	#  NODE DOF FACTOR
	sp  101   2    1.0
	sp  301   2   -1.0
}

# Create a recorder which writes to Node.out and prints
# the current load factor (pseudo-time) and dof 1 displacements at node 2 & 3
#recorder Node  -file ExternalNodes.out  disp  -node  101  301  -dof  2
recorder Element  101	199  -file  Element1.out  force
recorder Element  199  -file  Element1Section.out  -section 7 deformations

recorder Element  1         -file  Joint.out  defoANDforce
recorder Element  1         -file  JointDamage.out  damage
damageModel NormalizedPeak	5	0.01	-0.01	defo
recorder ElementDamage	199	-file	Element1Damage.out	-section	6	-dof	1	-damage	5


# Load control with variable
#                             
integrator LoadControl .05 100 .005 .5

# Convergence test
#                  tolerance maxIter displayCode
test NormDispIncr  1.0e-06     25         1

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Penalty   1.0e12   1.0e12

# System of equations solver
system SparseGeneral -piv

# Analysis for lateral load
analysis Static

# Perform the displacement history analysis
set CurrentStep    0
set MaxStep     10000
set ok             0

while {$ok == 0 && $CurrentStep < $MaxStep} {
	set ok [analyze 1 ]
	
	if {$ok != 0} {
		puts "Failed - Change the algorithm to linesearch"
		integrator ArcLength	0.05	1.0
		algorithm NewtonLineSearch    0.6
		set ok [analyze 1 ]
		test NormDispIncr  1.0e-06     25         1
		integrator LoadControl .05 100 .005 .5
		algorithm Newton
    	    }

	if {$ok != 0} {
		puts "Failed - Change the algorithm to linesearch, reduced tolerance to 1.0e-5"
		test NormDispIncr  1.0e-05     25         1
		algorithm NewtonLineSearch    0.6
		set ok [analyze 1 ]
		test NormDispIncr  1.0e-06     25         1
		algorithm Newton
    	    }

	if {$ok != 0} {
		puts "Failed - Change the algorithm to linesearch, reduced tolerance to 1.0e-4"
		test NormDispIncr  1.0e-04     25         1
		algorithm NewtonLineSearch    0.6
		set ok [analyze 1 ]
		test NormDispIncr  1.0e-06     25         1
		algorithm Newton
    	    }

	set CurrentStep [expr $CurrentStep + 1]
}

if { $ok != 0 } {
	puts "Sorry the analysis failed at step [expr $CurrentStep] "
}