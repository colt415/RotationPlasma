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
Field3D Ni;

//3D total values
Field3D Nit, phi, rho;

//E x B drift velocity
Vector3D vE;

//Flags
int phi_flags;
//Field routines, inverse laplacian
int solve_phi_tridag(Field3D &r,Field3D &p,int flags);

int physics_init(bool restarting){
  
/*  mesh->get(Ni0,"Ni0");
  mesh->get(rho0,"rho0");
  mesh->get(phi0,"Phi0");
*/
  
  GRID_LOAD(Ni0);
  GRID_LOAD(phi0);
  GRID_LOAD(rho0);

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("2fluid");
  OPTION(options, phi_flags,   0);

  bout_solve(rho,"rho");
  bout_solve(Ni,"Ni");
  
  phi=phi0;
//  rho=0.0;  
  if(!restarting){
   // phi+=phi0;
    rho+=rho0;
  }
  
  //Add other variables to be dumped to file
  dump.add(phi,"phi",1);
  
  dump.add(Ni0,"Ni0",0);
  dump.add(phi0,"phi0",0);
  dump.add(rho0,"rho0",0);
  
  return 0;
}


//just define a macro for V_E dot Grad
#define vE_Grad(f,p) (b0xGrad_dot_Grad(p,f))

int physics_run(BoutReal t){
  
  solve_phi_tridag(rho,phi,phi_flags);

  mesh->communicate(Ni,rho,phi);

  //Nit=Ni0;

  //Density equation
  ddt(Ni)=0.0;
  ddt(Ni)-=vE_Grad(Ni0,phi);
  ddt(Ni)-=vE_Grad(Ni,phi);
  
  //Vorticity
/*  
  vE.x=DDZ(phi)*sqrt(mesh->g33);
  vE.y=0.0;
  vE.z=-DDX(phi)*sqrt(mesh->g11);
  vE.covariant=false;
  ddt(rho)=0.0;
  ddt(rho)+=Div(V_dot_Grad(vE,Grad(phi)));
  ddt(rho)+=V_dot_Grad(V_dot_Grad(vE,Grad(phi)),Ni)/Ni0;
  ddt(rho)+=V_dot_Grad(V_dot_Grad(vE,Grad(phi)),Ni0)/Ni0;
*/
  ddt(rho)=0.0;
  ddt(rho)+=Ni*vE_Grad(rho,phi);
  
  return 0;
}

#include <invert_laplace.hxx>

// Performs inversion of rho (r) to get phi (p)
int solve_phi_tridag(Field3D &r, Field3D &p, int flags)
{
  //output.write("Solving phi: %e, %e -> %e\n", max(abs(r)), min(Ni0), max(abs(r/Ni0)));

  if(invert_laplace(r, p, flags, NULL)) {
    return 1;
  }

  //Field3D pertPi = Ti*Ni0 + Ni*Ti0;
  //p -= pertPi/Ni0;
  return(0);
}

