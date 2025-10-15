set outfile [open Hnum_CaNi.dat w]

set sel1 [atomselect top "segname A"]
set sel2 [atomselect top "segname C"]
set nf [molinfo top get numframes]

for { set i 0 } { $i <= $nf } { incr i } { 
	$sel1 frame $i
	$sel2 frame $i
	set nhb1 [llength [lindex [measure hbonds 3.5 30 $sel1 $sel2] 0]]
	set nhb2 [llength [lindex [measure hbonds 3.5 30 $sel2 $sel1] 0]]
	set nhb [expr $nhb1 + $nhb2]
	puts $outfile "$i $nhb"
	# puts $outfile "\n"
}
close $outfile
