##################################################################
## Units&Constants.tcl#
## Defining units and constants#
## Li Ning#
##2004/10/26#
##################################################################

# ======================================
# Outputs All With International Units
# ======================================

# Constant
set PI [expr 4*atan(1.)]
set g  9.81
set U 1.e20; # a really large number
set u [expr 1/$U]; # a really small number

# Basic units
set m   1.
set sec 1.
set N   1.

# Other units
# Angle
set rad 1.
set deg [expr $PI/180.*$rad]

# Length
set mm  0.001
set mm2 [expr $mm*$mm]
set mm4 [expr $mm2*$mm2]
# Force
set kN  [expr 1000.*$N]

# Pressure
set MPa [expr $N/pow($mm,2)]
##################################################################
# ===============================================
# Units in US
# ===============================================
set in [expr 25.4*$mm]; 			# define basic units
set kip [expr 0.4536*1000*$g*$N]; 		# define basic units
set kips $kip;
set ft [expr 12.*$in]; 				# define engineering units
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set in2 [expr $in*$in]; 			# inch^2
set in4 [expr $in2*$in2]; 			# inch^4
set cm [expr $in/2.54];				# centimeter, needed for displacement input in MultipleSupport excitation
set lbf [expr $kip/1000];
set pcf [expr $lbf/$ft/$ft/$ft];
set psf [expr $lbf/$ft/$ft];






