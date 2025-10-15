
set fil [open rgyr_lig-2.dat w] 
                                   
set nf [molinfo top get numframes]
set sel [atomselect top "segname C and noh"]

# rmsd calculation loop
for {set i 0} {$i < $nf } { incr i } {
	$sel frame $i
	puts $fil "[measure rgyr $sel weight mass]"
}
close $fil


