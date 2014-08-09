reset
set noxtics
set noytics
set title "ACP :  Depenses de l'Etat - 1872-1971" 
set key off
set xlabel "PC1"
set ylabel "PC2"
plot [-4:3.5] [-3.5:3]  "pc.dat" index 0 using 4:3:1 with labels font "Times, 7"
