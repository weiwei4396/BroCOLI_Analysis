echo "strat at: `date +%T`"
old_time=`date +%s`

# SIRV direct RNA sample1 annotations=100%
#/data/workdir/panw/softwares/stringtie/stringtie-2.2.1/stringtie -L -G /data/workdir/panw/SIRV_Set2/SIRV_100.gtf -e --rf -o /data/workdir/panw/SIRV_Set2/Quantification/Stringtie/updated.gtf /data/workdir/panw/SIRV_Set2/samples_123_sorted.bam

/data/workdir/panw/softwares/stringtie/stringtie-2.2.1/stringtie -L -G /data/workdir/panw/Data/SIRV_Set2/SIRV_100.gtf --rf /data/workdir/panw/Data/SIRV_Set2/mini230_sorted.bam > counts_transcript.txt

new_time=`date +%s`
echo "running time: $(($new_time - $old_time))s"
echo "end at: `date +%T`"
