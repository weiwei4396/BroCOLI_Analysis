
echo "strat at: `date +%T`"
old_time=`date +%s`

Rscript Bambu.R 

new_time=`date +%s`
echo "running time: $(($new_time - $old_time))s"
echo "end at: `date +%T`"
