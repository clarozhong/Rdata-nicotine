 set sel [atomselect top "protein and serial 2494 or chain C and serial 3825"]
     
     set outfile [open dist_ILNi.dat w]
     
     set nf [molinfo top get numframes]

     for { set i 0 } { $i <= $nf } { incr i } { 
          $sel frame $i
          $sel update
    
    set coords [$sel get {x y z}]
    set atom1 [lindex $coords 0]
    set atom2 [lindex $coords 1]
    
    set dist {}
    
    lappend dist [veclength [vecsub $atom2 $atom1]]
    puts $outfile $dist
    }
    close $outfile

