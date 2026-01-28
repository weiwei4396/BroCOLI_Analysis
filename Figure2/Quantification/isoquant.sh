echo "strat at: `date +%T`"
old_time=`date +%s`

isoquant.py --reference /data/workdir/panw/Data/SIRV_Set2/SIRV.fasta \
  --genedb /data/workdir/panw/Data/SIRV_Set2/SIRV_100.gtf \
  --bam /data/workdir/panw/Data/SIRV_Set2/mini230_sorted.bam \
  --data_type nanopore -o /data/workdir/panw/Data/SIRV_Set2/Quantification/Isoquant

new_time=`date +%s`
echo "running time: $(($new_time - $old_time))s"
echo "end at: `date +%T`"



