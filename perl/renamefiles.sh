#!/bin/sh
# INCOMPLETE!!!
# Rename the bunch of files created inside folder run0? using the
# ADAF code to some other preferred name.
# Example:
# 	Say you run this code inside folder run01 and call this code as
#	renamefiles.sh test
#	It will do the following renaming of files:
# 	out_01 spec_01 spec_01_ spec_01_ssd
#	out_teste spec_teste spec_teste_ spec_teste_ssd

# Gets name of current directory
currentdir=`pwd`



          if test-commands; then
            consequent-commands;
          [elif more-test-commands; then
            more-consequents;]
          [else alternate-consequents;]
          fi

# Depending on the folder in which you are (run01, run02 etc), take the
# appropriate action
if pwd | grep run01 ; then
  cp out_01 


if [ $currentdir=" ] ; then

for X in `ls *.mp3`
    do
    TARGET=`basename $X .mp3`
    mpg123 -s $X | sox -t raw -r 44100 -s -w -c 2 - ${TARGET}.wav
    done
