#!/bin/bash

echo "Cleaning old files..."
if [ -e out.log ]; then
  mv -v out.log out.log.bak
fi
if [ -e server.out ]; then
  mv -v server.out server.out.bak
fi
if [ -e clients.out ]; then
  mv -v clients.out clients.out.bak
fi
if [ -e auto-rockstar.cfg ]; then
  rm -v auto-rockstar.cfg
fi
if [ $(ls halos/* 2> /dev/null | wc -l) != "0" ]; then
  rm -rv halos/*
fi

echo "Submitting run script..."
echo "qsub run_rockstar.pbs"
qsub run_rockstar.pbs

