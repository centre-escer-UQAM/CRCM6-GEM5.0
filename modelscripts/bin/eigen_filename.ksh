#!/bin/ksh
#
nmlfile=$1

if [ -n "${nmlfile}" ] ; then
if [ -e ${nmlfile} ] ; then
   GRDTYP=$(rpy.nml_get -f ${nmlfile} grid/grd_typ_s 2> /dev/null | sed "s/'//g")
   if [ -z "${GRDTYP}" ] ; then exit ; fi

   if [ "$GRDTYP" == "GY" ] ; then
      dim=$(rpy.nml_get -f ${nmlfile} grid/grd_nj      2> /dev/null)
      if [ -z "${dim}" ] ; then
         dim=$(rpy.nml_get -f ${nmlfile} grid/grd_ni   2> /dev/null)
      fi
      OVERLAP=$(rpy.nml_get -f ${nmlfile} grid/Grd_overlap 2> /dev/null)
      if [ -z "${OVERLAP}" ] ; then OVERLAP=0.0 ; fi
      echo ${GRDTYP} ${dim} ${OVERLAP} > $TMPDIR/gem_grid$$
   else
      GNI=$(rpy.nml_get -f ${nmlfile} grid/grd_ni   2> /dev/null)
      DX=$(rpy.nml_get -f ${nmlfile} grid/grd_dx   2> /dev/null)
      LATR=$(rpy.nml_get -f ${nmlfile} grid/Grd_latr 2> /dev/null)
      if [ -z "${LATR}" ] ; then LATR=0.0 ; fi
      LONR=$(rpy.nml_get -f ${nmlfile} grid/Grd_lonr 2> /dev/null)
      if [ -z "${LONR}" ] ; then LONR=180.0 ; fi
      echo ${GRDTYP} ${GNI} ${DX} ${LATR} ${LONR} > $TMPDIR/gem_grid$$
   fi

   fn=eigenv_v1_$(md5sum $TMPDIR/gem_grid$$ | cut -d" " -f 1)
   /bin/rm -f $TMPDIR/gem_grid$$
   echo $fn
fi
fi
