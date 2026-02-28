#!/usr/bin/bash
set -e

INPUT="Tyler.png"
NPROC=`nproc`
JFLOATY="../../jFloaty.exe"

if [[ ! -x ${JFLOATY} ]]; then
    echo "${JFLOATY} is not executable or not present. Please run make in project directory."
    exit -1
fi

echo "Processing images using up to ${NPROC} jobs"

# Function definition (must come before it is called)
process_q() {
  local q=$1
  local s=$2
  local g=$3
  local d=$4
  JPG="${d}/temp${q}.jpg"
  DJBMP="${d}/temp${q}_djpeg.bmp"
  JFPNG="${d}/temp${q}_jfloaty.png"
  DJPNG="${d}/temp${q}_djpeg.png"
  ../../jFloaty.exe -i ${INPUT} -q ${q} -subsample ${s} -d ${g} -o ${JPG} -i ${JPG} -m ${g} -clip -o16 ${JFPNG} 2>&1 > /dev/null
  djpeg -bmp ${JPG} > ${DJBMP}
  ../../jFloaty.exe -i ${DJBMP} -m ${g} -clip -o16 ${DJPNG}
  rm ${DJBMP}
  echo ${JFPNG} ${DJPNG}
}

for s in $(seq 0 1); do
for g in $(seq 1 4 5); do

  if [ ${s} == "0" ]; then
    yuv="444"
  else
    yuv="420"
  fi
  DIR="gain${g}_${yuv}"

  mkdir -p ${DIR}

  ../../jFloaty.exe -i ${INPUT} -d ${g} -dither8 -m ${g} -o "${DIR}/temp_dither8.png" -clip 2>&1 > /dev/null &
  ../../jFloaty.exe -i ${INPUT} -d ${g} -floor8  -m ${g} -o "${DIR}/temp_floor8.png"  -clip 2>&1 > /dev/null &
  ../../jFloaty.exe -i ${INPUT} -d ${g} -round8  -m ${g} -o "${DIR}/temp_round8.png"  -clip 2>&1 > /dev/null &

  for i in $(seq -w 100 -1 1); do
    process_q ${i} ${s} ${g} ${DIR} &
    while [ `jobs | wc -l` -ge $NPROC ]
    do
      sleep 1
    done
  done

  wait

  echo "Image processing complete. Generating dssim report."
  dssim ${INPUT} ${DIR}/temp*_*.png | tee ${DIR}/report.txt

  echo "DSSIM complete. Plotting"
  python plotter.py ${DIR}
done
done
