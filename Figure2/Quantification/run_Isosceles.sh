
echo "strat at: `date +%T`"
old_time=`date +%s`

Rscript Isosceles.R 

new_time=`date +%s`
echo "end at: `date +%T`"
echo "running time: $(($new_time - $old_time))s"
