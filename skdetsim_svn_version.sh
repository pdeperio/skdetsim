#!/bin/sh
if test -e skdetsim_svn_version.h 
then
  rm skdetsim_svn_version.h
  echo replacing skdetsim_svn_version.h
else
  echo creating skdetsim_svn_version.h
fi
printf "       CHARACTER*16 SKDETSIM_SVN_VERSION /'" >> skdetsim_svn_version.h
svnversion -n . >> skdetsim_svn_version.h
printf "'/\n" >> skdetsim_svn_version.h
touch skidentify.F
