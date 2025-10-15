     set outfile [open rmsf_lig.txt w]
     set sel [atomselect top "segname C or segname A"]

set nf [molinfo top get numframes]
set ref [atomselect top "segname C or segname A" frame 0]
set comp [atomselect top "segname C or segname A"]
set all [atomselect top all]

for {set i 1} {$i<$nf} {incr i} {
    $comp frame $i
    $all frame $i
    $all move [measure fit $comp $ref]
}

	set alpharmsf [measure rmsf $sel first 0 last -1]
	set res [$sel get resid]

     for { set i 0 } { $i < [llength $res] } { incr i } { 
	puts $outfile "[lindex $res $i]\t[lindex $alpharmsf $i]"
     }

     close $outfile 