#!/bin/bash
#$ -l mem=24G,time=1:: -S /bin/bash -N PLEIO -j y 
echo $NSLOTS
echo 'NUMSLOTS'
#ps -eo nlwp | tail -n +2 | awk '{ num_threads += $1 } END { print num_threads }'
wkd='/ifs/scratch/msph/eigen/hl3565/01_MTGB/codes/source/pleio'
${wkd}/pleio.py --metain ${wkd}/example/input.txt.gz --sg ${wkd}/example/sg.txt.gz --ce ${wkd}/example/re.txt.gz --nis 10000 --parallel --create --out ${wkd}/test
