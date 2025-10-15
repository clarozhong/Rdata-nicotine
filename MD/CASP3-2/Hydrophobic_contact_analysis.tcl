proc hydrophobic_contact_occupancy { {dis 4.0} {sel1 protein} {sel2 none} {outfile stdout} {mol top} } {
  if {! [string compare $dis help]} {
     puts "Usage: hydrophobic_contact_occupancy distance selection1 selection2 mol"
     return
  }

  if {! [string compare $mol top]} {
     set mol [molinfo top]
  }

  if {[string compare $outfile stdout]} {
     set outfile [open $outfile w];
  }

  set hbondallframes {}
  set hbondcount {}

  set sumofall 0
  set sumof60 0
  set sumof50 0
  set sumof10 0
  set numof60 0
  set numof50 0
  set numof10 0

  # Select hydrophobic residues for K3 and B2
  set sel1 [atomselect $mol "segname A"]
  set sel2 [atomselect $mol "segname C"]

  set framenumberbackup [molinfo $mol get frame]
  set numberofframes [molinfo $mol get numframes] 

  for { set i 0 } { $i < $numberofframes } { incr i } {
      molinfo $mol set frame $i

      set hbondsingleframe [measure contacts $dis $sel1 $sel2]

      for { set j 0 } { $j < [llength [lindex $hbondsingleframe 0] ] } { incr j } {
          set newhbond {}
          lappend newhbond [lindex $hbondsingleframe 0 $j ]
          lappend newhbond [lindex $hbondsingleframe 1 $j ] 
          set hbondexist [lsearch $hbondallframes $newhbond]
          if { $hbondexist == -1 } {
             lappend hbondallframes $newhbond
             lappend hbondcount 1
          } else {
             lset hbondcount $hbondexist [expr { [lindex $hbondcount $hbondexist] + 1 } ]
          }
      }
  }

  for { set i 0 } { $i < [llength $hbondallframes] } { incr i } {
      set donor [atomselect $mol "index [lindex $hbondallframes $i 0]"]
      set acceptor [atomselect $mol "index [lindex $hbondallframes $i 1]"]
      set occupancy [expr {100*[lindex $hbondcount $i]/($numberofframes+0.0) } ]

      set sumofall [expr { $sumofall + $occupancy} ]

      if { $occupancy >= 50 } { 
         set sumof50 [expr { $sumof50 + $occupancy}]
         set numof50 [expr { $numof50 + 1}]
      }

      if { $occupancy >= 60 } { 
         set sumof60 [expr { $sumof60 + $occupancy}]
         set numof60 [expr { $numof60 + 1}]
      }

      if { $occupancy >= 10 } { 
         set sumof10 [expr { $sumof10 + $occupancy}]
         set numof10 [expr { $numof10 + 1}]
      }

      if {$i == 0 } {
         puts $outfile "donor \t acceptor \t occupancy(%)"
      }
      puts $outfile [format "%s-%s%i-%s-%i \t %s-%s%i-%s-%i \t %.2f" [$donor get segname] [$donor get resname] [$donor get resid] [$donor get name] [$donor get index] [$acceptor get segname] [$acceptor get resname] [$acceptor get resid] [$acceptor get name] [$acceptor get index] $occupancy]
  }

  if {[string compare $outfile stdout]} {
     close $outfile
  }

  molinfo $mol set frame $framenumberbackup
}

########################################################################
set pro [atomselect top "segname A"]
set lig [atomselect top "segname C"]
hydrophobic_contact_occupancy 3.5 $pro $lig hydrophobic_complex.dat
