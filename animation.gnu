reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
unset key
set border 0
set view 342,200
set zrange [0:8]
set title "Evolution of diffusion equation over time."
set zlabel "Phi(x,y,t)"
set ylabel "y"
set xlabel "x"
set grid
#set dgrid3d 30,30, 1
#set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
#set style data line
cmdStr = "ls -l out* | wc -l"
numFiles = system(cmdStr)
stats "out0.txt"
do for [ii=1:STATS_blocks] {
	    splot for [i=0:numFiles-1] 'out'.i.'.txt' index (ii-1) using 1:2:3
		pause 0.1
} 
