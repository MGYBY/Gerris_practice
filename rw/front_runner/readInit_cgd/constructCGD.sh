awk 'BEGIN{ ORS=" " }{ print $1 }END{ print "\n" }' < t_25.4.txt >> h.cgd
awk '{ print $2 }' < t_25.4.txt >> h.cgd


awk 'BEGIN{ ORS=" " }{ print $1 }END{ print "\n" }' < t_25.4.txt >> hu.cgd
awk '{ print $3 }' < t_25.4.txt >> hu.cgd
