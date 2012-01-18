proc GeneratePeaks {Dmax {DincrStatic 0.01} {CycleType "Full"} {Fact 1} } {;	# generate incremental disps for Dmax
	###########################################################################
	## GeneratePeaks $Dmax $DincrStatic $CycleType $Fact 
	###########################################################################
	# generate incremental disps for Dmax
	# this proc creates a file which defines a vector then executes the file to return the vector of disp. increments
	# by Silvia Mazzoni, 2006
	# input variables
	#	$Dmax	: peak displacement (can be + or negative)
	#	$DincrStatic	: displacement increment (optional, default=0.01, independently of units)
	#	$CycleType	: Full (0->+peak), Half (0->+peak->0), Full (0->+peak->0->-peak->0)   (optional, def=Full)
	#	$Fact	: scaling factor (optional, default=1)
	#	$iDstepFileName	: file name where displacement history is stored temporarily, until next disp. peak
	# output variable
	#	$iDstep	: vector of displacement increments
	set outFileID [open tmpDsteps.tcl w]
	set Disp $DincrStatic;
	puts $outFileID "set iDstep {\n0.0";
	puts $outFileID $Disp; #puts $outFileID $Disp;	# open vector definition and some 0
	set Dmax [expr $Dmax*$Fact];	# scale value
	if {$Dmax<0} {;  # avoid the divide by zero
		set dx [expr -$DincrStatic]
	} else {
		set dx $DincrStatic;
	}
	set NstepsPeak [expr int(abs($Dmax)/$DincrStatic)]
	for {set i 1} {$i < $NstepsPeak} {incr i 1} {;		# zero to one
		set Disp [expr $Disp + $dx]
		puts $outFileID $Disp;			# write to file
	}
	if {$CycleType !="Push"} {
		for {set i 1} {$i <= $NstepsPeak} {incr i 1} {;		# one to zero
			set Disp [expr $Disp - $dx]
			puts $outFileID $Disp;			# write to file
		}
		if {$CycleType !="HalfCycle"} {
			for {set i 1} {$i <= $NstepsPeak} {incr i 1} {;		# zero to minus one
				set Disp [expr $Disp - $dx]
				puts $outFileID $Disp;			# write to file
			}
			for {set i 1} {$i < $NstepsPeak} {incr i 1} {;		# minus one to zero
				set Disp [expr $Disp + $dx]
				puts $outFileID $Disp;			# write to file
			}
		}
	}
	puts $outFileID " 0.0\n}";		# close vector definition
	close $outFileID
	source tmpDsteps.tcl;		# source tcl file to define entire vector
	return $iDstep
}