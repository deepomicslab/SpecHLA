dir=/mnt/d/HLAPro_backup/HLAPro/com_base_error/
bash $dir/run.blast.baseerror.sh $1
perl $dir/merge.novel.baseerror1.pl
perl $dir/merge.novel.baseerror2.pl
perl $dir/count_fre2.pl

