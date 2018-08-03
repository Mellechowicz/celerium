reset
set term epslatex standalone color size 4, 4
set loadpath './'
set out './wykres.tex'

set multiplot

NXPLOTS = 1
NYPLOTS = 1

XSEP = 0.0
YSEP = 0.00

TMARGIN = 0.98
BMARGIN = 0.1
LMARGIN = 0.1
RMARGIN = 0.8

PLOTHEIGHT = (TMARGIN-BMARGIN - (NXPLOTS - 1)*YSEP )/NYPLOTS
PLOTWIDTH = (RMARGIN-LMARGIN - (NYPLOTS - 1)*XSEP )/NXPLOTS 

left_m(i, j) = LMARGIN + i*(PLOTWIDTH + XSEP)
right_m(i, j) = left_m(i, j) + PLOTWIDTH
top_m(i, j) = TMARGIN - j*(PLOTHEIGHT + YSEP)
bottom_m(i, j) = top_m(i, j) - PLOTHEIGHT




file = 'data.dat'

####################################

i = 0
j = 0

set tmargin at screen top_m(i, j)
set bmargin at screen bottom_m(i, j)
set lmargin at screen left_m(i, j)
set rmargin at screen right_m(i, j)

set key width -4
set key vertical maxrows 2
set key at graph 0.90, 0.98
set key samplen 1.2
set key spacing 1.3
unset key


#set ylabel '$d^2/d^2_\mathrm{HF}$' offset 1, 0, 0
#set xlabel 'Site intex' offset 0, 0.5, 0

#set xrange [0.99 : 40.001]
#set yrange [0.4 : 0.65]
#set xtics (1, 10, 20, 30, 40)
#set ytics -10, 0.1, 10
#set mxtics 2
#set mytics 2

#set style line 100 lt 4 lw 2 pt 5 lc rgb 'black' dt 2


set cbrange [-0.3 : 0.3]
set palette defined ( 0 "blue", 1 "white", 2 "red")

set view map
set size square
set pm3d
#set dgrid3d 20, 20
#set isosamples 10, 10

#plot file palette w points ti ''

splot file palette pt 5 ps 0.7







unset multiplot

set out 
