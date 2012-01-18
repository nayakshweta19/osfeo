# NTW1q1.tcl all shell element with quad3d element, no colomn elements, so there 3DOF/node
# the appended dofs there is fix in x-, y- axials
# use FAReinforcedConcrtePlaneStress material
wipe;
source Units&Constants.tcl;
######################## 
# Analysis-Sequence  1 #
######################## 

# Start of model generation 
# ========================= 
set LunitTXT "mm"
set FunitTXT "N"
# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  3  -ndf  3

# Define material(s) 
# ------------------ 
# Materials.tcl 

# Material "ElasticDefault":    matTag    E    eta  
uniaxialMaterial  Elastic       1  +2.000000E+005  +0.000000E+000 

# Material "CoreConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       2  -3.500000E+001  -2.000000E-003  -4.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "CoverConcrete":    matTag    fc'    epsc0    fcu'    epsu    lambda    ft    Ets 
uniaxialMaterial  Concrete02       3  -3.000000E+001  -2.000000E-003  -3.500000E+001  -1.400000E-002  +5.000000E-001  +6.000000E-001  +5.000000E+002 

# Material "reinforcedBar":    matTag    Fy    E    b    R0    cR1    cR2    <a1    a2    a3    a4>    <sig0> 
uniaxialMaterial  Steel02       4  +3.750000E+002  +2.000000E+005  +5.000000E-002  +1.850000E+001  +9.250000E-001  +1.500000E-001  +0.000000E+000  +1.000000E+000  +0.000000E+000  +1.000000E+000  +0.000000E+000 

# Material "Steel":    matTag    E    v    rho 
nDMaterial  ElasticIsotropic3D       5  +2.000000E+005  +2.500000E-001  +0.000000E+000 

# Material "TortionMat":    matTag    E    eta  
uniaxialMaterial  Elastic       6  +2.000000E+009  +0.000000E+000 

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------

# set fc fy E   ; #Units: L=mm, F=N, 
set wfc [expr     32.0];
set wfy [expr    430.0];
set wE  [expr    2.0e5];
set rou1 0.0059;
set rou2 0.0026;
set t [expr $in*6.0/$mm];

# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2

# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
uniaxialMaterial ConcreteZ02  13 [expr -$wfc] -0.0025  
uniaxialMaterial ConcreteZ02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025  
#uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 

set pi [expr 2.0*asin(1.0)];
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0 $eps0
nDMaterial FAReinforcedConcretePlaneStress 15  0.0 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.0025

#section PlateFiber secTag  ndMatTag  h
section  PlateFiber  2   15  $t

# Define section(s) 
# ----------------- 
# Section "ElasticDefault":    secTag    E    A    Iz  Iy  G  J 
section  Elastic       1  +2.000000E+005  +1.800000E+002  +4.860000E+003  +1.500000E+003  +1.115400E+004  +3.916000E+003 

# Section "PlateFiber":    secTag    matTag    h 
section  PlateFiber       3       5  +1.000000E-002 

# Define geometry 
# --------------- 
# generate the nodes and elements
set nxFlage 2;  # divied by 2, for middle wall
set nyFlage 4; # divided by 4, for 4 story
set flageNodeStart 1;
set flageEleStart 1;
#element  quad3d       1       1       11       15       5     $t  "PlaneStress"   15;
# block2d $nx $ny $e1 $n1 element (element arguments) {
set cmd "block2D $nxFlage $nyFlage $flageNodeStart $flageEleStart quad3d  \" $t  PlaneStress2D  5 \"  {
    1   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    2   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    3   [expr 5.4E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
    4   [expr 0.0E+0*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
}"
eval $cmd;

set nxWeb 1;  
set nyWeb 4;  # divided by 4, for 4 story
set webNodeStart [expr ($nxFlage+1)*($nyFlage+1)+1];
set webEleStart  [expr $nxFlage*$nyFlage+1];
#element  quad3d       9      16       6       10      20     $t  "PlaneStress"   15;
set cmd "block2D $nxWeb $nyWeb $webNodeStart $webEleStart quad3d  \" $t  PlaneStress2D  5 \"  {
    1   [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 0.00E+0*$in/$mm]
    2   [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 0.00E+0*$in/$mm]
    3   [expr 2.7E+1*$in/$mm]  [expr 0.0E+0*$in/$mm]  [expr 2.88E+2*$in/$mm]
    4   [expr 2.7E+1*$in/$mm]  [expr 7.8E+1*$in/$mm]  [expr 2.88E+2*$in/$mm]
}"
eval $cmd;

set nodeIDs [getNodeTags]
set endNodeTag [lindex $nodeIDs end]
set eleIDs [getEleTags]
set endEleTag [lindex $eleIDs end]

set fileModel [open "qmodel.ops" "w"]
puts $fileModel "\nmodel  BasicBuilder  -ndm  3  -ndf  3 \n"
foreach nid  $nodeIDs {
  puts $fileModel "node $nid [nodeCoord $nid]"
}

# ---------------------------------- 
# Define Single Point Constraints 
# ---------------------------------- 
for { set i 0 } { $i <= $nxWeb } { incr i 1 } {
  fix [expr $webNodeStart+$i] 1 1 1;
  puts $fileModel "fix [expr $webNodeStart+$i] 1 1 1;"
}
for { set i 0 } { $i <= $nxFlage } { incr i 1 } {
  fix [expr $flageNodeStart+$i] 1 1 1;
  puts $fileModel "fix [expr $flageNodeStart+$i] 1 1 1;"
}

if { $nyWeb != $nyFlage } {
  puts "ERROR, nyWeb != $nyFlage ! "
}
# Constraints; Define nodal masses
for { set i 1 } { $i <= $nyFlage} {incr i 1} {
  for { set j 0 } { $j <= $nxFlage} {incr j 1} {
    fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0
    puts $fileModel "fix [expr $flageNodeStart+($nxFlage+1)*$i+$j] 0 1 0"
    mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1000. 1000.
    puts $fileModel "mass [expr $flageNodeStart+($nxFlage+1)*$i+$j] 1000. 1000. 1000."
  }

  for { set j 0 } { $j <= $nxWeb} {incr j 1} {
    fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0
    puts $fileModel "fix [expr $webNodeStart+($nxWeb+1)*$i+$j] 1 0 0"
    mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1000. 1000. 1000.
    puts $fileModel "mass [expr $webNodeStart+($nxWeb+1)*$i+$j] 1000. 1000. 1000."
  }
  remove  sp [expr $webNodeStart+($nxWeb+1)*($i+1)-1] 1
  remove  sp [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 2
  equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3
  puts $fileModel "remove  sp [expr $webNodeStart+($nxWeb+1)*($i+1)-1] 1"
  puts $fileModel "remove  sp [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 2"
  puts $fileModel "equalDOF  [expr $webNodeStart+($nxWeb+1)*($i+1)-1] [expr $flageNodeStart+($nxFlage+1)*$i+$nxFlage/2] 1 2 3"

}

set eleIDs [getEleTags]
set endEleTag [lindex $eleIDs end]
foreach nid  $eleIDs {
  puts $fileModel "element quad $nid [eleNodes $nid]"
}
close $fileModel

# ---------------------------------- 
# Define time series 
# ---------------------------------- 
# TimeSeries "LinearDefault":    tsTag    cFactor 
timeSeries  Linear       1  -factor  +1.000000E+000 

# ---------------------------------- 
# Start of anaysis generation 
# ---------------------------------- 
# ---------------------------------- 
# Get Initial Stiffness 
# ---------------------------------- 
initialize 

# ---------------------------------- 
# Analysis: StaticDefaultCase 
# ---------------------------------- 
# ---------------------------------- 
# Define load pattern 
# ---------------------------------- 

# LoadPattern "PlainDefault":    patternTag    tsTag 
pattern  Plain       1       1  { 
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $nxFlage} {incr j 1} {
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004 
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004 
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }

    # SP    nodeTag    dofTag    DispValue 
 
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
}

# ---------------------------------- 
# Define recorder(s) 
# ---------------------------------- 

# Node Recorder "DefoShape":    fileName    <nodeTag>    dof    respType 
recorder  Node  -file  StaticDefaultCase_Node_DefoShape_Dsp.out  -time -nodeRange 1 $endNodeTag  -dof  1  2  3  disp 

# ---------------------------------- 
# Define analysis options 
# ---------------------------------- 

# AnalysisOptn "StaticDefault": Type: Static 
# ------------------------------------------ 
# Constraint Handler Lagrange
constraints  Plain;  
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

analyze     1 

# ---------------------------------- 
# Reset for next analysis case 
# ---------------------------------- 
loadConst -time 0.0 
remove recorders 
wipeAnalysis 

# ---------------------------------- 
# Analysis: EigenDefaultCase 
# ---------------------------------- 
# ---------------------------------- 
# Define recorder(s) 
# ---------------------------------- 

# Node Recorder "EigenVector":    fileName    <nodeTag>    dof    respType 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_1.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen1 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_2.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen2 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_3.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen3 
recorder  Node -file  EigenDefaultCase_Node_EigenVector_EigenVec_4.out  -time  -nodeRange 1 $endNodeTag  -dof  1  2  3  eigen4 
# ---------------------------------- 
# Define analysis options 
# ---------------------------------- 
# AnalysisOptn_2.tcl 
# ---------------------------------- 
# AnalysisOptn "EigenDefault": Type: Eigen 
# ---------------------------------- 
# Constraint Handler Lagrange
constraints  Transformation
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  Newton 
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

#set eigFID [open EigenDefaultCase_Node_EigenVector_EigenVal.out w] 
set lambda [eigen 15] 
#puts $eigFID $lambda
#close $eigFID 
record;
foreach lam $lambda {
  if { $lam > 0.0 } {
    lappend omega [expr sqrt($lam)]
    lappend f [expr sqrt($lam)/(2*$pi)]
    lappend T [expr (2*$pi)/sqrt($lam)]
    puts "Tn = [expr (2*$pi)/sqrt($lam)] s";
  } else {
    lappend omega -1.0
    lappend f "complex"
    lappend T "complex"
    puts "Tn = complex";
  }
}
# get values of eigenvectors for translational DOFs
#---------------------------------------------------
set fileEig [open "eigFile.eig" "w"]
for { set k 1 } { $k <= 15 } { incr k } {
  set tempStr "[lindex $T [expr $k-1]]"
  foreach tg $nodeIDs {
    lappend tempStr [lindex [nodeEigenvector $tg $k 1] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 2] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 3] 0]
  }
  puts $fileEig $tempStr;
}
close $fileEig 
# ---------------------------------- 
# Clean up 
# ---------------------------------- 
remove recorders;
#wipe 

# ---------------------------------- 
# Reset for next analysis case 
# ---------------------------------- 
loadConst -time 0.0 
wipeAnalysis 
# ---------------------------------- 
# Analysis: PushoverCase 
# ---------------------------------- 
# ---------------------------------- 
# Define time series 
# ---------------------------------- 

# TimeSeries "TimeSeriesX":    tsTag    dt    filePath    cFactor 
timeSeries  Path       2  -dt  +1.000000E-002  -filePath  TimeSeriesX.thf  -factor  +1.000000E-003 

# TimeSeries "TimeSeriesY":    tsTag    dt    filePath    cFactor 
timeSeries  Path       3  -dt  +1.000000E-002  -filePath  TimeSeriesY.thf  -factor  +1.000000E-003 

# LoadPattern_3.tcl 

# LoadPattern "PushoverPattern":    patternTag    tsTag 
##pattern  Plain       3    "Linear"  { 
##    # Load    nodeTag    LoadValues 
## 
##    # SP    nodeTag    dofTag    DispValue 
##    sp 5  1 1.0
##    sp 10 1 1.0
##    sp 15 1 1.0
##    sp 25 1 1.0
##    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
## 
##    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
##} 

pattern  Plain       2    "Linear"  { 
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $nxFlage} {incr j 1} {
      load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +1.000000E+000  0.000000E+000 
      #puts "load [expr $flageNodeStart+($nxFlage+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    for { set j 0 } { $j <= $nxWeb} {incr j 1} {
      load [expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +1.000000E+000  0.000000E+000 
      #puts "[expr $webNodeStart+($nxWeb+1)*$nyFlage+$j] +0.000000E+000  +0.000000E+000  -4.000000E+004 "
    }
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
} 

recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.out  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Node  -file  PushoverCase_Node_DefoShape_Dsp.dis  -time -nodeRange 1  $endNodeTag -dof  1  2  3  disp 
recorder  Element -file PushoverCase_Ele_All.stress -time -eleRange 1 $endEleTag stress
recorder  Element -file PushoverCase_Ele_All.strain -time -eleRange 1 $endEleTag strain

# ---------------------------------- 
# perform the analysis
# ---------------------------------- 
set IDctrlNode [expr $webNodeStart+($nxWeb+1)*$nyFlage];
set IDctrlDOF 2;
set Dincr 0.01;
# Constraint Handler Lagrange Transformation
constraints  Plain
# Convergence Test 
test  NormDispIncr  +1.000000E-003    25     0     2 
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  BandGeneral 
# Analysis Type 
analysis  Static 

# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;  # constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator
source GeneratePeaks.tcl;
#  ---------------------------------    perform Static Cyclic Displacements Analysis
set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";  # format for screen/file output of DONE/PROBLEM analysis
set iDmax "  1 3 5 10";
set CycleType Full;      # you can do Full / Push / Half cycles with the proc
set Ncycles 1;      # specify the number of cycles at each peak
set Fact 1.0;

foreach Dmax $iDmax {
  set iDstep [GeneratePeaks $Dmax $Dincr $CycleType $Fact];  # this proc is defined above
  for {set i 1} {$i <= $Ncycles} {incr i 1} {
    set zeroD 0
    set D0 0.0
    foreach Dstep $iDstep {
      set D1 $Dstep
      set Dincr [expr $D1 - $D0]
      integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
      analysis Static
      # ----------------------------------------------first analyze command------------------------
      set ok [analyze 1]
      # ----------------------------------------------if convergence failure-------------------------
      if {$ok != 0} {
        # if analysis fails, we try some other stuff
        # performance is slower inside this loop  global maxNumIterStatic;      # max no. of iterations performed before "failure to converge" is ret'd
        if {$ok != 0} {
          puts "Trying Newton with Initial Tangent .."
          test NormDispIncr   $Tol 2000 0
          algorithm Newton -initial
          set ok [analyze 1]
          test $testTypeStatic $TolStatic      $maxNumIterStatic    0
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          puts "Trying Broyden .."
          algorithm Broyden 8
          set ok [analyze 1 ]
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          puts "Trying NewtonWithLineSearch .."
          algorithm NewtonLineSearch 0.8 
          set ok [analyze 1]
          algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
          set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
          puts $putout
          return -1
        }; # end if
      }; # end if
      # -----------------------------------------------------------------------------------------------------
      set D0 $D1;      # move to next step
    }; # end Dstep
  };    # end i
};  # end of iDmaxCycl
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
  puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
  puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}

