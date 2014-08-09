set object 1 circle center 0, 0, 0 size 0.7 arc [0:360]
set object 1 fc rgb "green" 
set object 1 fc rgb "green" fillstyle empty border
set object 2 circle size 1 fc rgb "navy" fillstyle empty border
set view equal xy
unset xtics; unset ytics
set yrange [-1.1:1.1]
set size ratio 1 1
set xzeroaxis linetype 0
set yzeroaxis linetype 0
set yrange [-1.2:1.2]
unset border
set title "Budget de l'Etat"
set yrange [2:2]
set xlabel "PC1"
set ylabel "PC2"
set key off

plot [-1:1] [-1:1] [-1:1] "cercle-budget" using 2:3:1 with labels font "Times,8"

