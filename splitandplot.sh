for i in $(seq -f "%03g" 1 202)
do
   cat ../mu1e5-100by100by300.sum | head -n $i | tail -n 1 > temporary.raw
   /groups/astro/ianpaga/anaconda3/bin/python plotaframe.py $i
done

