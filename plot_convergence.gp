set terminal pngcairo 
set output 'convergence.png'

set title "Historique de Convergence de l'Algorithme de Richardson"
set xlabel "Itérations"
set ylabel "Erreur (Résidu)"
set grid

set logscale y 

plot "bin/RESVEC.dat" using 1 with linespoints title "Erreur par Itération"

