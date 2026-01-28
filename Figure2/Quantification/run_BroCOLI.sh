
SAM=/data/workdir/panw/Data/SIRV_Set2/mini230_sorted.sam
fa=/data/workdir/panw/Data/SIRV_Set2/SIRV.fasta
GTF=/data/workdir/panw/Data/SIRV_Set2/SIRV_100.gtf

g++ -std=c++11 -pthread /data/workdir/panw/C_code/source/BroCOLI_bulk_20251217.cpp -I /data/workdir/panw/C_code/source/ -o BroCOLI

echo "strat at: `date +%T`"
old_time=`date +%s`

./BroCOLI -t 1 -r 0 -s $SAM -f $fa -g $GTF -o /data/workdir/panw/Data/SIRV_Set2/Quantification/BroCOLI/OUT

new_time=`date +%s`
echo "end at: `date +%T`"
echo "running time: $(($new_time - $old_time))s"

