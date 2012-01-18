puts " -------------------------2D Model------------------------" 

puts " -------------------------------------------------------------- " 

puts " RC Wall Section" 

puts " -------------------------------------------------------------- " 

# Define model builder 
# -------------------- 
model basic -ndm 2 -ndf 3 
set secTag 1 

#_______________________________________________________ 
#_______________________________________________________ 


# define UNITS ---------------------------------------------------------------------------- 
set in 1.; # define basic units -- output units 
set kip 1.; # define basic units -- output units 
set sec 1.; # define basic units -- output units 
set LunitTXT "inch"; # define basic-unit text for output 
set FunitTXT "kip"; # define basic-unit text for output 
set TunitTXT "sec"; # define basic-unit text for output 
set ft [expr 12.*$in]; # define engineering units 
set ksi [expr $kip/pow($in,2)]; 
set psi [expr $ksi/1000.]; 
set lbf [expr $psi*$in*$in]; # pounds force 
set pcf [expr $lbf/pow($ft,3)]; # pounds per cubic foot 
set psf [expr $lbf/pow($ft,3)]; # pounds per square foot 
set in2 [expr $in*$in]; # inch^2 
set in4 [expr $in*$in*$in*$in]; # inch^4 
set cm [expr $in/2.54]; # centimeter, needed for displacement input in MultipleSupport excitation 
set PI [expr 2*asin(1.0)]; # define constants 
set g [expr 32.2*$ft/pow($sec,2)]; # gravitational acceleration 
set Ubig 1.e10; # a really large number 
set Usmall [expr 1/$Ubig]; # a really small number 


#_______________________________________________________ 
#_______________________________________________________ 

# set up parameters for section and element definition 
#-------------------------------------------------- 
set IDcore 1; # ID tag for core concrete 
set IDcover 2; # ID tag for cover concrete 
set IDsteel 3; # ID tag for steel 

# Define materials for nonlinear section 
# ------------------------------------------ 
#CONCRETE------------------------------------------------------ 
# Confined concrete: (CORE CONCRETE) 
set kfc 1.3;	 # Ratio between confined and unconfined compressive strength 
set fc [expr -6*$ksi]; # CONCRETE Compressive Strength, ksi (+Tension, -Compression) 
set Ec [expr 57*$ksi*sqrt(-$fc/$psi)]; # Concrete Elastic Modulus 
set fc1C [expr $kfc*$fc]; # CONFINED concrete (mander model), maximum stress 
set eps1C [expr 2.*$fc1C/$Ec]; # strain at maximum stress 
set fc2C [expr 0.2*$fc]; # ultimate stress 
set eps2C [expr 2*$eps1C]; # strain at ultimate stress 

uniaxialMaterial Concrete01 $IDcore $fc1C $eps1C $fc2C $eps2C; 

# Unconfined concrete: (COVER CONCRETE) 
set fc1U $fc; # UNCONFINED concrete (todeschini parabolic model), maximum stress 
set eps1U -0.003; # strain at maximum stress 
set fc2U [expr 0.1*$fc]; # ultimate stress 
set eps2U -0.006; # strain at ultimate stress 

uniaxialMaterial Concrete01 $IDcover $fc1U $eps1U $fc2U $eps2U; 

# STEEL -------------------------------------------------------------- 
# Reinforcing steel 
set fy [expr 60*$ksi];	 # STEEL yield stress 
set Es	[expr 29000.*$ksi];	 # modulus of steel 
set Bs	0.01;	 # strain-hardening ratio 
set R0 18;	 # control the transition from elastic to plastic branches 
set cR1 0.925;	 # control the transition from elastic to plastic branches 
set cR2 0.15;	 # control the transition from elastic to plastic branches 

uniaxialMaterial Steel02 $IDsteel $fy $Es $Bs $R0 $cR1 $cR2; 

# Define rebar diameter & area per bar size 

set db3 [expr 0.375*$in] 
set db4 [expr 0.5*$in] 
set db5 [expr 0.625*$in] 
set db6 [expr 0.75*$in] 
set db7 [expr 0.875*$in] 
set db8 [expr 1.0*$in] 
set db9 [expr 1.128*$in] 
set db10 [expr 1.27*$in] 
set db11 [expr 1.41*$in] 

set Ab3 [expr 0.11*$in2] 
set Ab4 [expr 0.2*$in2] 
set Ab5 [expr 0.31*$in2] 
set Ab6 [expr 0.44*$in2] 
set Ab7 [expr 0.6*$in2] 
set Ab8 [expr 0.79*$in2] 
set Ab9 [expr 1.0*$in2] 
set Ab10 [expr 1.27*$in2] 
set Ab11 [expr 1.56*$in2] 

# Define cross-section geometry for section 
# ------------------------------------------ 

puts "Section Dimensions" 
# set some paramaters 
set bf [expr 210*$in]; # flange width 
set h [expr 138*$in]; # nominal depth 
set bw [expr 30*$in]; # web thickness 
set tf [expr 24*$in]; # flange thickness 
set bbf [expr 74*$in];	 # width of boundary element 1 into flange 
set tbe1 [expr 50*$in];	 # depth of boundary element 1 in to flange end of web 
set tbe2 [expr 32*$in];	 # depth of boundary element 2 into far end of web 
set cover [expr 2$in];	 # Clear concrete cover over rebar 

puts "bf = $bf in" 
puts "h = $h in" 
puts "bw = $bw in" 
puts "tf = $tf in" 

# some variables derived from the parameters 
set y1 [expr $tbe1-$tf] 
set y2 [expr $h-$tf] 
set y3 [expr $y2-$tbe2] 
set z1 [expr $bf-$bw] 
set z2 [expr $bbf-$bw] 

puts "y1 = $y1 in" 
puts "y2 = $y2 in" 
puts "y3 = $y2 in" 
puts "z1 = $z1 in" 
puts "z2 = $z2 in" 

#Gross area of the section 
set Ag [expr $z1*$tf+$h*$bw] 
puts "Ag = $Ag in^2" 

#------------------------------------------------------ 
#Boundary Element definitions 
#	 __________ 
#	| 1________|	
#	| |a 
#	| |	
#	| |	 1 is location of boundary element 1 (1a is in web section, 1b is in flange section) 
#	| |	 2 is location of boundary element 2 
#	| |	
#	|2|	
#	|_|	 a is the location of the origin (0,0) to define for wall patches (inside corner of web and flange) 
# 
#	 y 
#	 | 
#	 | 
#	 z____| 

#Define Steel in Boundary Element 1a (web section) 

#set rho1a 0.01;	 # steel ratio in boundary element 1 (in flange/web) 
#set As1a [expr $rho1a*$bw*$tbe1];	 # Area of steel in BE given rho for the BE 

# 2 rows of bars with (13) #10 bars per row 
set nb1a 13;	 # Number of bars per row (y direction) 

#Define Steel in Boundary Element 1b (flange section) 

#set rho1b 0.01;	 # steel ratio in boundary element 1 (in flange/web) 
#set As1b [expr $rho1b*$z2*$tf];	 # Area of steel in BE given rho for the BE 

# 2 rows of bars with (5) #10 bars per row 
set nb1b 5;	 # Number of bars per row (z direction) 

#Define Steel in Boundary Element 2 (in web end) 

#set rho2 0.01;	 # steel ratio in boundary element 2 (web end) 
#set As2 [expr $rho2*$bw*$tbe2];	 # Area of steel in BE given rho for the BE 

# 2 rows of bar with (6) #10 bars per row 
set nb2 6;	 # Number of bars per row (y direction) 

#----------------------------------------------------------------- 
#Section model set up and discritization 

#Define discritization size of concrete patches 
set nfdb 20; # number of fibers along boundary element depth 
set nftb 20; # number of fibers along boundary element thickness 
set nfbf 20; # number of fibers along flange/web width 
set nftf 20; # number of fibers along flange/web thickness 

# Discretize the section with layers and patches 
section Fiber $secTag { 

# Define boundary elements (confined) 
# ID, nfIJ, nfJK, yI, zI, yJ, zJ, yK, zK, yL, zL 
patch quad $IDcore $nftb $nfdb -$y2 $bw -$y2 0 -$y3 0 -$y3 $bw;	 #BE 2 
patch quad $IDcore $nftb $nfdb -$y1 $bw -$y1 0 $tf 0 $tf $bw;	 #BE 1a 
patch quad $IDcore $nftb $nfdb 0 0 0 -$z2 $tf -$z2 $tf 0;	 #BE 1b 
patch quad $IDcore $nftb $nfdb 0 0 0 -$z2 $tf -$z2 $tf 0;	 #BE 1b 

# Define web/flange elements (unconfined) 
patch quad $IDcover $nfbf $nftf -$y3 $bw -$y3 0 -$y1 0 -$y1 $bw 
patch quad $IDcover $nfbf $nftf 0 -$z2 0 -$z1 $tf -$z1 $tf -$z2 


# Define Bottom boundary element reinforcement 
# ID, #bars, Abars, y start, z start, y end, z end 
layer straight $IDsteel $nb1a $Ab10 [expr $tf-$cover] $cover -$y1 $cover; #BE 1a Row 1 
layer straight $IDsteel $nb1a $Ab10 [expr $tf-$cover] [expr $bw-$cover] -$y1 [expr $bw-$cover]; #BE 1a Row 2 
layer straight $IDsteel $nb1b $Ab10 $cover 0 $cover -$z2;	 #BE 1b row 1 
layer straight $IDsteel $nb1b $Ab10 [expr $tf-$cover] 0 [expr $tf-$cover] -$z2;	 #BE 1b row 2 
layer straight $IDsteel $nb2 $Ab10 -$y3 [expr $bw-$cover] [expr -$y2+$cover] [expr $bw-$cover];	 #BE 2 row 1 
layer straight $IDsteel $nb2 $Ab10 -$y3 $cover [expr -$y2+$cover] $cover;	 #BE 2 row 2 
} 



#_______________________________________________________ 
#_______________________________________________________ 


# set AXIAL LOAD -------------------------------------------------------- 
set P [expr 0.1*$fc*$Ag*$kip]; # + Tension, - Compression (NOTE: $f'c is negitive, as previously defined) 
puts "Axial Load on Section = $P kip" 

# set Maximum Curvature: 
set Ku [expr 0.001/$in]; 
set numIncr 100; # Number of analysis increments to maximum curvature (default=100) 
puts "Number of analysis increments: $numIncr 
" 

#_______________________________________________________ 
#_______________________________________________________ 


# ANALYSIS of section -------------------------------------------------- 

# Define two nodes at (0,0) 
node 1001 0.0 0.0 
node 1002 0.0 0.0 

# Fix all degrees of freedom except axial and bending 
fix 1001 1 1 1 
fix 1002 0 1 0 

# Define element 
# tag ndI ndJ secTag 
element zeroLengthSection 2001 1001 1002 $secTag -orient 1 0 0 0 1 0 

# Create Moment Curvature recorder ----> output moment (col 1) & curvature (col 2) 
recorder Node -xml Mphi.out -time -node 1002 -dof 3 disp; 

#Create Stress-Strain recorder ---> col1 = Moment (k-in), col2 = stress (ksi), col3 = strain (in/in) 
recorder Element -xml SSFlangeSteel1.out -time -ele 2001 section fiber $tf $cover $IDsteel stressStrain;	# Top steel recorder 1 
recorder Element -xml SSFlangeSteel2.out -time -ele 2001 section fiber $tf -$z1 $IDsteel stressStrain;	 # Top steel recorder 2 
recorder Element -xml SSFlangeConc1.out -time -ele 2001 section fiber $tf $cover $IDcore stressStrain;	# Top concrete recorder 1 
recorder Element -xml SSFlangeConc2.out -time -ele 2001 section fiber $tf -$z1 $IDcore stressStrain;	 # Top concrete recorder 2 
recorder Element -xml SSWebSteel.out -time -ele 2001 section fiber -$y2 $cover $IDsteel stressStrain;	 # bottom steel recorder 
recorder Element -xml SSWebConc.out -time -ele 2001 section fiber -$y2 $cover $IDcore stressStrain;	 # bottom concrete recorder 

# Define constant axial load 
pattern Plain 3001 Constant { 
load 1002 $P 0.0 0.0 
} 

# Define analysis parameters 
integrator LoadControl 0 1 0 0 
system SparseGeneral -piv; # Overkill, but may need the pivoting! 
test EnergyIncr 1.0e-1 100 
numberer Plain 
constraints Plain 
algorithm Newton 
analysis Static 

# Do one analysis for constant axial load 
analyze 1 

# Define reference moment 
pattern Plain 3002 "Linear" { 
load 1002 0.0 0.0 1.0 
} 

# Compute curvature increment 
set dK [expr $Ku/$numIncr] 

# Use displacement control at node 1002 for section analysis, dof 3 
integrator DisplacementControl 1002 3 $dK 1 $dK $dK 

# Do the section analysis 
set ok [analyze $numIncr] 
puts "Analysis Run" 

# ----------------------------------------------if convergence failure------------------------- 
set IDctrlNode 1002 
set IDctrlDOF 3 
set Dmax $Ku 
set Dincr $dK 
set TolStatic 1.e-2; 
set testTypeStatic EnergyIncr 
set maxNumIterStatic 6 
set algorithmTypeStatic Newton 
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Curv=%.4f /%s"; # format for screen/file output of DONE/PROBLEM analysis 

if {$ok != 0} { 
# if analysis fails, we try some other stuff, performance is slower inside this loop 
set Dstep 0.0; 
set ok 0 
while {$Dstep <= 1.0 && $ok == 0} { 
set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ] 
set Dstep [expr $controlDisp/$Dmax] 
set ok [analyze 1]; # this will return zero if no convergence problems were encountered 
if {$ok != 0} {; # reduce step size if still fails to converge 
set Nk 4; # reduce step size 
set DincrReduced [expr $Dincr/$Nk]; 
integrator DisplacementControl $IDctrlNode $IDctrlDOF $DincrReduced 
for {set ik 1} {$ik <=$Nk} {incr ik 1} { 
set ok [analyze 1]; # this will return zero if no convergence problems were encountered 
if {$ok != 0} { 
# if analysis fails, we try some other stuff 
# performance is slower inside this loop global maxNumIterStatic; # max no. of iterations performed before "failure to converge" is ret'd 
puts "Trying Newton with Initial Tangent .." 
test NormDispIncr $TolStatic 2000 0 
algorithm Newton -initial 
set ok [analyze 1] 
test $testTypeStatic $TolStatic $maxNumIterStatic 0 
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
if {$ok != 0} {; # stop if still fails to converge 
puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT] 
return -1 
}; # end if 
}; # end for 
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr; # bring back to original increment 
}; # end if 
}; # end while loop 
}; # end if ok !0 
# ----------------------------------------------------------------------------------------------------- 
global LunitTXT; # load time-unit text 
if { [info exists LunitTXT] != 1} {set LunitTXT "Length"}; # set blank if it has not been defined previously. 
if {$ok != 0 } { 
puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT] 
} else { 
puts [format $fmt1 "DONE" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT] 
} 

puts "Analysis Finished" 
puts " -------------------------------------------------------------- "

