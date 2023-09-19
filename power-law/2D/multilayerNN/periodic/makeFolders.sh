preString=lp_
mkdir lp_1 lp_2
for i in $(seq 4 4 20)
do
    mkdir ${preString}$i
done
