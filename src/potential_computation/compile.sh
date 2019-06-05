#!/bin/bash

gfortran simulation_module.f90 io_module.f90 pbc_displacement_module.f90 energy_module.f90 potential.f90 -o potential.exe

rm *mod

