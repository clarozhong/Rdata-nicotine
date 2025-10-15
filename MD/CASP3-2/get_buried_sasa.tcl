set outfile [open sasa_CaNi_w.dat w]

set nf [molinfo top get numframes] 
set pro [atomselect top all] 
set C [atomselect top "segname C"]
set A [atomselect top protein]


puts $outfile "Frames\tsa_pro\tsa_A\tsa_C\tsa_buried"

for { set i 0 } { $i < $nf } { incr i } { 
$pro frame $i 
$A frame $i
$C frame $i

set sapro [measure sasa 1.4 $pro]
set saA [measure sasa 1.4 $A] 
set saC [measure sasa 1.4 $C]
set saburied [expr $saA+$saC-$sapro]

puts $outfile "$i\t$sapro\t$saA\t$saC\t$saburied" 
} 

close $outfile