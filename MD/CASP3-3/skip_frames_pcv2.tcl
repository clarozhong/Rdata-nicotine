mol load psf CaNi_ion.psf



animate read dcd  equ1-1.dcd skip 20 waitfor all



set pro [atomselect top "segname A or segname C"]



animate write dcd  equ1-1_skip20.dcd sel $pro



quit
