INPUT="80486dx2-large.jpg"
COUNT=2000
for i in $(seq -w 1 $COUNT)
do
  OUTPUT="temp${i}.jpg"
  #echo ${INPUT} ${OUTPUT}
  ./jFloaty.exe -i ${INPUT} -o ${OUTPUT} > /dev/null
  md5sum ${OUTPUT}
  INPUT=${OUTPUT}
done
