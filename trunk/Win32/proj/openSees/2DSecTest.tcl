# --------------------------------------------------------------------------------------------------
# build a section
#		Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source LibUnits.tcl;			# define units

# MATERIAL parameters -------------------------------------------------------------------
set SecTagFlex 2;			# assign a tag number to the column flexural behavior
set SecTagAxial 3;			# assign a tag number to the column axial behavior	
set SecTag 1;			# assign a tag number to the column section tag

# COLUMN section
# calculated stiffness parameters
set EASec $Ubig;				# assign large value to axial stiffness
set MySec [expr 130000*$kip*$in];		# yield moment
set PhiYSec [expr 0.65e-4/$in];		# yield curvature
set EICrack [expr $MySec/$PhiYSec];		# cracked section inertia
set b 0.01 ;				# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
#uniaxialMaterial Steel01 $SecTagFlex $MySec $EICrack $b; 		# bilinear behavior for flexural moment-curvature
#uniaxialMaterial Elastic $SecTagAxial $EASec;			# this is not used as a material, this is an axial-force-strain response
#section Aggregator $SecTag $SecTagAxial P $SecTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)

# set fc fy E
set wfc 28.0;
set wfyv 560;
set wfyh1 450;
set wE 190000.0;
set rou1 0;
set rou2 0;
set rouv 0.0246;
set rouh1 0.00246; # #10@70

# tag fy E0 fpc rou
uniaxialMaterial SteelZ01 11 $wfyv $wE $wfc $rouv
uniaxialMaterial SteelZ01 12 $wfyh1 $wE $wfc $rouh1
# UniaxialMaterial: concreteZ01
# ConcreteZ01 tag f¡¯c ec0
uniaxialMaterial ConcreteL02 14 [expr -$wfc] -0.003
uniaxialMaterial ConcreteL02 15 [expr -$wfc] -0.003
set pi 3.141592654
set nDMatTag 21
# NDMaterial: FAFourSteelPCPlaneStress
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
#nDMaterial FAFourSteelPCPlaneStress $nDMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

#           RAFourSteetPCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?
#nDMaterial RAFourSteetPCPlaneStress $nDMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAFourSteelRCPlaneStress  (100)
#                                     tag rho p1 p2 s1 s2 c1 c2 angle1 angle2 angle3 angle4 roup1 roup2 rous1 rous2 fpc fpy fy E0 epsc0?
nDMaterial FAFourSteelRCPlaneStress $nDMatTag 0.0 11 12 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] [expr 0.0*$pi] [expr 0.5*$pi] $rou1 $rou2 $rouv $rouh1 $wfc 0.0 $wfyv $wE 0.003
#nDMaterial RAFourSteelRCPlaneStress $nDMatTag 0.0 11 12 11 12 14 15 [expr 1.0*$pi] [expr 0.966*$pi] [expr 0.5*$pi] [expr 0.0*$pi] $rou1 $rou2 $rouv $rouh1 0.0 0.0 $wfc $wfyv $wfyv $wE 0.003

# NDMaterial: FAReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlaneStress $nDMatTag 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

#           ReinforcedConcretePlaneStress matTag? rho? s1? s2? c1? c2? angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?
#nDMaterial ReinforcedConcretePlaneStress $nDMatTag 0.0 11 12 14 15 [expr 0.5*$pi]  [expr 0.0*$pi]  $rouv $rouh1  $wfc $wfyv $wE 0.003

# NDMaterial: CSMMRCPlaneStress  (100)
#                                     tag rho s1 s2 c1 c2 angle1 angle2 rous1 rous2 fpc  fy  E0 epsc0
#nDMaterial CSMMRCPlaneStress $nDMatTag 0.0 11 12 14 15 [expr 0.0*$pi] [expr 0.5*$pi] $rouv $rouh1 $wfc $wfyv $wE 0.002

set nDMatTag1 6;
# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D $nDMatTag1 2e4 0.2  0 
#           RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? angle1? angle2? rou1? rou2? pstrain? fpc? fyT? fy? E0? epsc0?

set nDBeamTag 10;
#nDmaterial for beam fiber 3D
nDMaterial BeamFiber $nDBeamTag $nDMatTag1;

set t1 35.0;
set NStrip1 1; # thickness 1
set t2 80.0;
set NStrip2 0; # thickness 2
set t3 42.0;
set NStrip3 1; # thickness 3
set np 4; # int. points
set C 0.4; # center of rotation 

#section definition 
#section WSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <shape?>
section WSection2d $SecTag $nDBeamTag  200 10 20 40 20 20 0.5
#section CSMMFiber2d $SecTag -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 21;
#	fiber2d  10 0 1 21;
#};
#section FiberInt 4 -NStrip $NStrip1 $t1 $NStrip2 $t2 $NStrip3 $t3 { ;
#	fiber2d -10 0 1 2;
#	fiber2d  10 0 1 2;
#};

# --------------------------------------------------------------------------------------------------
# Moment-Curvature analysis of section
#		Silvia Mazzoni & Frank McKenna, 2006
#

# define procedure
source MomentCurvature2D.tcl

# set AXIAL LOAD --------------------------------------------------------
set P [expr -10.0];	# + Tension, - Compression

# set maximum Curvature:
set Ku [expr 1.1];
set numIncr 100;	# Number of analysis increments to maximum curvature (default=100)
# Call the section analysis procedure
MomentCurvature2D $SecTag $P $Ku $numIncr