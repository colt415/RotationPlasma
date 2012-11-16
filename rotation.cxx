//Rotation plasma with simplified 2-fluid equation

#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//2D initial profiles
Field2D Ni0, phi0, rho0;

//3D evolving fields
Field3D rho, Ni, phi;




