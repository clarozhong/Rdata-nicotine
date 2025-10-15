set outfile [open "rmsd_lig-2.dat" w]
set nf [molinfo top get numframes]
set frame0 [atomselect top "segname C" frame 0]
set sel [atomselect top "segname C"]
set body0 [atomselect top "segname C" frame 0]
set body [atomselect top "segname C"]
set all [atomselect top all]

# Debugging output
puts "Number of selected atoms in frame 0: [$frame0 num]"
puts "Number of selected atoms in current frame: [$sel num]"

# rmsd calculation loop
for {set i 1} {$i < $nf} {incr i} {
    $sel frame $i
    $body frame $i
    $all frame $i
    $all move [measure fit $sel $frame0]
    puts $outfile "$i\t[measure rmsd $body $body0]"
}
close $outfile