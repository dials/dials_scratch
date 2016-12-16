
for t in "1:100" "4:20" "5:20" "6:200"
do
  i=$(echo $t | cut -d":" -f1)
  n=$(echo $t | cut -d":" -f2)
  for j in $(seq $n)
  do
    result=$(dials.python test_reflection_data.py ~/Desktop/profile_modelling_experiment/data/${i}_${j}_indexed.pickle  ~/Desktop/profile_modelling_experiment/data/${i}_${j}_experiments.json)
    echo $i $j $result
  done
done
