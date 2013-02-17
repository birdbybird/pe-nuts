#! /bin/bash

cmd="gcc -Wall -pedantic -std=c99 rudin-shapiro.c bitbitter.c -o pe384"

if [ "$1" == "-log" ]
then
  echo "logging"
  $cmd > compile_log 2>&1
  # ( $cmd 2>&1 ) > compile_log
else
  $cmd
fi
