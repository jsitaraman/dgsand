for ds in 1.25 1.0 0.625 0.5 0.3125 0.25;do
  sed 's/0.3125/'"$ds"'/g' makeStructGrid_p1.m.temp > makeStructGrid_p1.m
  octave makeStructGrid_p1.m
  ../tools/a.out grid1.dat
  ../tools/a.out grid2.dat
  ../src/dgsand input.dgsand1 input.dgsand2 "$ds" > p2log_"$ds"
  ../src/dgsand input.dgsandp1 input.dgsandp2 "$ds" > p1log_"$ds"
done
