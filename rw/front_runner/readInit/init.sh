# count the line number of the water depth field
lines=`wc -l initH | awk '{print $1}'`

# initial water depth file
# x y H
# awk -v lines=$lines '
# BEGIN{ for (x = 1; x <= 32; x += 1) 
# BEGIN{
#   print lines " 0 0"
# } {
# print $1 " " (x-0.50)*(3.0/32.0) " " $3;
# }
# }' #< initH | delaunay > initH.gts

awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  for (x = 1; x <= 32; x += 1) print $1 " " (x-0.50)*(3.0/32.0) " " $3;
}' < initH | delaunay > test.gts