#!/bin/bash
rm -f *.o
rm -f *.mod
gfortran -ffree-line-length-none -c constant.f90
gfortran -ffree-line-length-none -c func.f90
gfortran -ffree-line-length-none -c pick_event.f90
gfortran *.o -o pickevent
