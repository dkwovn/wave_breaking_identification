#!/bin/bash
rm -f *.o
rm -f *.mod
gfortran -ffree-line-length-none -c constant.f90
gfortran -ffree-line-length-none -c drct.f90
gfortran -ffree-line-length-none -c detect.f90
gfortran -ffree-line-length-none -c reverse.f90
gfortran -ffree-line-length-none -c main.f90
gfortran *.o -o mkct

