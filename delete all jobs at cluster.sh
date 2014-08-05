# delete all jobs
qstat | awk '$3 ~ "michal" {cmd="qdel " $1; system(cmd); close(cmd)}'