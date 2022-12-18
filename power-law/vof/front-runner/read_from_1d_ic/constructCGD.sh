awk 'BEGIN{ ORS=" " }{ print $1 }END{ print "\n" }' < t100_u >> u.cgd
awk '{ print $2 }' < t100_u >> u.cgd
