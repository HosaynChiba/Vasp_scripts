# Gnuplot script file for plotting data in file "*.dat"
# This file is called   estilo.p

#=======================================================
#     Line wide point size and font aspects
#=======================================================
#
#     Specify the line wide and dot size

  linw=1.6
  poin=1
  fsize=14

#======================================================
#          Output settings terminals
#=====================================================
#
#    Type of file you want to print the plot
#    eps is the most recomended
#    Default: Shows it on your screen

#set term post eps size 7.00 cm ,4.00 cm enhanced "Times-Roman, 14"
#set output "outputfile.eps"

#======================================================
#    Standard initialization to create the plot
#    some tags are for latin-american speakers
#======================================================

 set termopt enhanced         # Allows to put ^super/_{sub} indexes
 set   autoscale              # scale axes automatically
 unset log                    # remove any log-scaling
 unset label                  # remove any previous labels
 set xtic auto                # set xtics automatically
 set ytic auto                # set ytics automatically
 set encoding iso_8859_1      # Allow to put accents 
# set grid                     # Set the type of grid

#======================================================
#     SPECIFY PLOT RANGE
#======================================================

 set xr [-4:4]
 set xtics 1
 set yr [-100:100]
 set xzeroaxis lt -1
 set yzeroaxis lt 0 lw 1.5
 set xtics nomirror
 set ytics nomirror
#========================================================
#    KEYWORDS to control the PLOT LABELS
#========================================================
#    There are two ways to use the predertermine positions
#    or to indicate the x,y position using 
#    set key nobox at x, y

 set key nobox at 3.8, 90      #default inside the plot

#Para poner las leyendas en dif posiciones gnuplot utiliza
# top/bottom/center  left/right/center box/nobox
# set key center right box


#========================================================
#        PLOT INSIDE THE TOTAL CANVAS
#========================================================
#     Controls the position of the plot on the canvas
#     you can add extra "air" on the sides
#            l-> left t-> top r-> range

# set lmargin                          
# set tmargin 
# set rmargin
# set bmargin 5


#=======================================================
#   TITLES AND XY TIC LABELS
#=======================================================
#      You can modify the font and size but
#      is not necesary because you specify that on 
#      the outprint format

 set title  "Density of state of Fe_{16}Te_{32}"
 set xlabel "Energy [eV]"
 set ylabel "Density of States [states/eV]" 
 set label "" at 2.0,80

#=======================================================
#       Line Styles 
#=======================================================

 set style line 5  lt 1 lw 1.5 pt 5  ps poin  lc rgb "gray"
 set style line 11 lt 1 lw 1.5 pt 11 ps poin  lc rgb "black"
 set style line 10 lt 1 lw 1.5 pt 10 ps poin  lc rgb "red" 
 set style line 12 lt 1 lw 1.3 pt 12 ps poin  lc rgb "green"
 set style line 13 lt 1 lw 1.5 pt 12 ps poin  lc rgb "blue"

#======================================================
#      Add information and arrows to the plot
#======================================================

# set label "cruce evitado" at 1.3,-0.7

#Label can be "nohead"
set arrow from -3.2, 50 to -3.2, 70 ls 11
set arrow from -3.2, -50 to -3.2, -70 ls 11


#======================================================
#           Fit data to equations 
#======================================================
#    
# f1(x)= (a1*x)/(a1*x+1)
# a1=0
# fit f1(x) 'langmuir.dat' using 1:2 via a1



#======================================================
#          Output settings terminals
#=====================================================
#
#    Type of file you want to print the plot
#    eps is the most recomended 
#    Default: Shows it on your screen

set term pngcairo size 1200 , 800 enhanced font "Times-New-Roman, 24"
set output "plot_tot.png"


#=======================================================
#             Plot instructions
#=======================================================

#    filedatax.dat is your data file
#    using tag means which versus which column you are plotting
#    title mean the label title on the key box

# plot    "b3lyp_tio2.dat"  using (($1*27.21)+5.48):($4/27.21) title 'Total DOS' with filledcurve y1=0 ls 5, \
#         "b3lyp_tio2.dat"  using (($1*27.21)+5.48):($2/27.21) title 'Ti (p)' with lines ls 11, \
#         "b3lyp_tio2.dat"  using (($1*27.21)+5.48):($3/27.21) title 'O (d)' with lines ls 10;
plot "tdos.dat" using 1:2 title 'TOT' with lines ls 11,\
     "tdos.dat" using 1:3 title '' with lines ls 11,\
     "PDOS_Fe_UP.dat" using 1:11 title 'Fe' with lines ls 10,\
     "PDOS_Fe_DW.dat" using 1:11 title '' with lines ls 10,\
     "PDOS_Te_UP.dat" using 1:11 title 'Te' with lines ls 13,\
     "PDOS_Te_DW.dat" using 1:11 title '' with lines ls 13;
     

 unset style line
 reset
