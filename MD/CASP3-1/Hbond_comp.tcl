##计算Occupancy的平均值，计算>50%的平均值，计算>10%的平均值
############################################################################

proc hbond_occupancy_5 { {dis 3.5} {ang 30} {sel1 protein} {sel2 none} {outfile stdout} {mol top} } {
  if {! [string compare $dis help]} {
     puts "Usage: hbond_occupancy distance angle selection1 selection2 mol"
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
##############用于计算Occupancy均值的变量
  	set sumofall 0
	set sumof60 0
	set sumof50 0
	set sumof10 0
	set numof60 0
	set numof50 0
	set numof10 0


  if { ! [string compare $sel1 protein]} {
     set sel1 [atomselect $mol $sel1]
  }

  set framenumberbackup [molinfo $mol get frame]
  set numberofframes [molinfo $mol get numframes] 

  for { set i 0 } { $i < $numberofframes } { incr i } {
      molinfo $mol set frame $i

      if {[string compare $sel2 none]} {
         set hbondsingleframe [measure hbonds $dis $ang $sel1 $sel2]
		# 反过来计算
		set hbondsingleframe1 [measure hbonds $dis $ang $sel2 $sel1]
		
		
      } else {
         set hbondsingleframe [measure hbonds $dis $ang $sel1]

      }

      for { set j 0 } { $j < [llength [lindex $hbondsingleframe 0] ] } { incr j } {
          set newhbond {}
          lappend newhbond [lindex $hbondsingleframe 0 $j ]
          lappend newhbond [lindex $hbondsingleframe 1 $j ] 
          #lappend newhbond [lindex $hbondsingleframe 2 $j ]  #屏蔽此行可将供体受体相同但氢原子不同的情况归于同一氢键，但同时也屏蔽掉了氢原子的信息
          set hbondexist [lsearch $hbondallframes $newhbond]
          if { $hbondexist == -1 } {
             lappend hbondallframes $newhbond
             lappend hbondcount 1
          } else {
             lset hbondcount $hbondexist [expr { [lindex $hbondcount $hbondexist] + 1 } ]
          }
      }

	# 将hbondsingleframe1中的值整合到$hbondallframes中
       for { set j 0 } { $j < [llength [lindex $hbondsingleframe1 0] ] } { incr j } {
	   set newhbond {}
          lappend newhbond [lindex $hbondsingleframe1 0 $j ]
          lappend newhbond [lindex $hbondsingleframe1 1 $j ] 
          lappend newhbond [lindex $hbondsingleframe1 2 $j ]
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
	#####屏蔽掉了氢原子的信息
	#set hydrogen [atomselect $mol "index [lindex $hbondallframes $i 2]"]
      set occupancy [expr {100*[lindex $hbondcount $i]/($numberofframes+0.0) } ]

	#############################################################################################################################
	#########对Occupancy进行分类（all、>50%、>10%）、求和
	
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
	##############################################################################################################################

      if { ([$donor get name] == "N") || ([$donor get name] == "CA") || ([$donor get name] == "C") || ([$donor get name] == "O") } {
         set donortype "Main"
      } else {
         set donortype "Side"
      } 
      if { ([$acceptor get name] == "N") || ([$acceptor get name] == "CA") || ([$acceptor get name] == "C") || ([$acceptor get name] == "O") } {
         set acceptortype "Main"
      } else {
         set acceptortype "Side"
      } 

      if {$i == 0 } {
         puts $outfile "donor \t acceptor \t occupancy(%)"
      }
      #####puts $outfile [format "%s%i-%s \t %s%i-%s \t %.2f%%" [$donor get resname] [$donor get resid] $donortype [$acceptor get resname] [$acceptor get resid] $acceptortype $occupancy]
	#输出原子index和name
	puts $outfile [format "%s-%s%i-%s-%i \t %s-%s%i-%s-%i \t %.2f" [$donor get segname] [$donor get resname] [$donor get resid] [$donor get name] [$donor get index] [$acceptor get segname] [$acceptor get resname] [$acceptor get resid] [$acceptor get name] [$acceptor get index] $occupancy]
  }

######################################################################################################################################
#计算、输出各类型Occupancy的平均值和数目

	#puts $outfile "\n"
	#puts $outfile [format "%.2f \t %.2f \t %.2f" $sumofall $sumof50 $sumof10]
	#puts $outfile "type \t average \t number"
	#puts $outfile [format "%s \t %.2f%% \t %i" "all" [expr $sumofall/([llength $hbondallframes])] [llength $hbondallframes]]
	#puts $outfile [format "%s \t %.2f%% \t %i" ">60%" [expr $sumof50/$numof60] $numof60]
	#puts $outfile [format "%s \t %.2f%% \t %i" ">50%" [expr $sumof50/$numof50] $numof50]
	#puts $outfile [format "%s \t %.2f%% \t %i" ">10%" [expr $sumof10/$numof10] $numof10]
#######################################################################################################################################

  if {[string compare $outfile stdout]} {
     close $outfile
  }

  molinfo $mol set frame $framenumberbackup
}

########################################################################
set A [atomselect top "protein and segname A and name \"\[ONH\].*\""]
set C [atomselect top "segname C and name \"\[ONH\].*\""]
hbond_occupancy_5 3.5 30 $A $C Hbond_pcv5-3_EM.dat


#######################################################################
