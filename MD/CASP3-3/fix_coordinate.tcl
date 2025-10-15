# % $Id: rmsd.tcl,v 1.3 2005/02/18 17:56:49 mbach Exp $

                                             

  set numframe [molinfo top get numframes]
  #set fir [atomselect top all frame first]
  #set b1 [measure minmax $fir]
  #set b11 [lindex $b1 0]
  #set b12 [lindex $b1 1]
  #set b2 [vecsub $b12 $b11]
  #set a [lindex $b2 0]
  #set b [lindex $b2 1]
  #set c [lindex $b2 2]
  set a 86.622614795  
  set b 77.2850064275  
  set c 89.8027887538   
  set sel [atomselect top "segname A"]
  $sel frame 0
  set coord [measure center $sel]

    set nx 0; set ny 0; set nz 0     

  for {set t 1} {$t < $numframe} {incr t} {
    $sel frame $t
    set oldcoord $coord
    set coord [measure center $sel] 
  
    set dx [lindex [vecsub $coord $oldcoord] 0]
    if {$dx > [expr "$a / 2"]} {incr nx -1}
    if {$dx < [expr "-1 * $a / 2"]} {incr nx}
    set dy [lindex [vecsub $coord $oldcoord] 1]
    if {$dy > [expr "$b / 2"]} {incr ny -1}
    if {$dy < [expr "-1 * $b / 2"]} {incr ny}
    set dz [lindex [vecsub $coord $oldcoord] 2]
    if {$dz > [expr "$c / 2"]} {incr nz -1}
    if {$dz < [expr "-1 * $c / 2"]} {incr nz}
    $sel moveby [list [expr "$nx * $a"] [expr "$ny * $b"] [expr "$nz * $c"]]
  }

 
  set sel [atomselect top "segname C"]
  $sel frame 0
  set coord [measure center $sel]


    set nx 0; set ny 0; set nz 0
  for {set t 1} {$t < $numframe} {incr t} {
    $sel frame $t
    set oldcoord $coord
    set coord [measure center $sel] 

    set dx [lindex [vecsub $coord $oldcoord] 0]
    if {$dx > [expr "$a / 2"]} {incr nx -1}
    if {$dx < [expr "-1 * $a / 2"]} {incr nx}
    set dy [lindex [vecsub $coord $oldcoord] 1]
    if {$dy > [expr "$b / 2"]} {incr ny -1}
    if {$dy < [expr "-1 * $b / 2"]} {incr ny}
    set dz [lindex [vecsub $coord $oldcoord] 2]
    if {$dz > [expr "$c / 2"]} {incr nz -1}
    if {$dz < [expr "-1 * $c / 2"]} {incr nz}
    $sel moveby [list [expr "$nx * $a"] [expr "$ny * $b"] [expr "$nz * $c"]]
  }



  

