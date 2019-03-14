#!/bin/ksh

set -x
# set -o noclobber is used to avoid writing in existing "content" files
#    This was needed for parallel execution of this script by the r1 prep_lbcs_loop task
set -o noclobber

rep_in=${1}
wild=${2}
rep_ou=${3}
nthreads=${4}

model_inputs ()
{
set -x
date
INPUT=$1
DEST=$2
lnk_dest=$3
bname=$(basename ${INPUT})
if [ "${lnk_dest}" != "${bname}" ] ; then lnk_dest=${lnk_dest}/${bname} ; fi
list_A=$(r.fstliste.new -izfst $INPUT -typvar "A" |cut -d ":" -f 11 | sort -u )
list_P=$(r.fstliste.new -izfst $INPUT -typvar "P" |cut -d ":" -f 11 | sort -u )
list_I=$(r.fstliste.new -izfst $INPUT -typvar "I" |cut -d ":" -f 11 | sort -u )

for i in $(echo "${list_A} ${list_P} ${list_I}" | tr ' ' '\n' | sort -u | tr '\n' ' ') ; do
  valid=$(echo $i | cut -c1-8).$(echo $i | cut -c9-14)
  dir=${DEST}/VALID_${valid}
  mkdir -p ${dir}
  ln -sf ../../${lnk_dest} ${dir}
done
}

if [ "${wild}" == "@NIL@" ] ; then
  wild=
fi

if [ -d $rep_in ] ; then

  cd $rep_in
  count=0

  for file in $(ls -1 ${wild}) ; do
    if [ ! -d $file ] ; then
      count=$(( count + 1 ))
      list=${TASK_WORK}/lis_$(basename $file)_$$.lis
      model_inputs $file ${rep_ou} MODEL_inrep 1> $list 2>&1 &
    fi
    if [ $count -eq ${nthreads} ] ; then
      wait
      count=0
    fi
  done
  wait

else
  model_inputs $rep_in ${rep_ou} ANALYSIS
fi

for i in ${rep_ou}/VALID_*  ; do
  if [ -d $i ] ; then
    cd $i
    cnt=$(ls -1 * 2> /dev/null | wc -l)
    set +e
    cat > ${TMPDIR}/content$$ <<EOF
$(echo $cnt)
$(ls -1 * 2> /dev/null)
EOF
    mv ${TMPDIR}/content$$ content
    set -e
  fi
done
