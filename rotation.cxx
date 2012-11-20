//Rotation plasma with simplified 2-fluid equation

#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//2D initial profiles
Field2D Ni0, phi0, rho0;

//3D evolving fields
Field3D Ni,phi,rho;

//3D total values
Field3D Nit, phit, rhot;

//E x B drift velocity
Vector3D vE;

//Flags
int phi_flags;
//Field routines, inverse laplacian
//int solve_phi_tridag(Field3D &r,Field3D &p,int flags);

//bool bout_exb,arakawa;

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
//  OPTION(options, bout_exb,    false);
//  OPTION(options, arakawa,    false);
 
  bout_solve(rho,"rho");
  bout_solve(Ni,"Ni");
  
/*  phi=phi0;
//  rho=0.0;  
  if(!restarting){
   // phi+=phi0;
    rho+=rho0;
  }*/
  
  //Add other variables to be dumped to file
  dump.add(phi,"phi",1);
  
  dump.add(Ni0,"Ni0",0);
  dump.add(phi0,"phi0",0);
  dump.add(rho0,"rho0",0);
  
  return 0;
}


//just define a macro for V_E dot Grad
#define vE_Grad(f,p) b0xGrad_dot_Grad(p,f)
// Routines for ExB terms (end of this file)
/*const Field2D vE_Grad(const Field2D &f, const Field2D &p);
const Field3D vE_Grad(const Field2D &f, const Field3D &p);
const Field3D vE_Grad(const Field3D &f, const Field2D &p);
const Field3D vE_Grad(const Field3D &f, const Field3D &p);
*/

int physics_run(BoutReal t){
  
//  solve_phi_tridag(rho,phi,phi_flags);
  phi=invert_laplace(rho/Ni0,phi_flags,NULL);
//  phi.applyBoundary();

  mesh->communicate(Ni,rho,phi);

  //Nit=Ni0;

  //Density equation
  ddt(Ni)=0.0;
  ddt(Ni)-=vE_Grad(Ni0,phi);
  ddt(Ni)-=vE_Grad(Ni,phi0);
  ddt(Ni)-=vE_Grad(Ni0,phi0);
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
/* 
  vE.x=DDZ(phi+phi0)*sqrt(mesh->g33);
  vE.y=0.0;
  vE.z=-DDX(phi+phi0)*sqrt(mesh->g11);
  vE.covariant=false;
  vE.applyBoundary();
  mesh->communicate(vE);
*/
  ddt(rho)=0.0;
//  ddt(rho)+=V_dot_Grad(vE,Grad(phi))*Grad(Ni)/Ni0;
//  ddt(rho)+=DDX(Ni)*vE_Grad(DDX(phi),phi)/Ni0;
//  ddt(rho)+=DDZ(Ni)*sqrt(mesh->g33)*vE_Grad(DDZ(phi)*sqrt(mesh->g33),phi)/Ni0;
  ddt(rho)-=vE_Grad(rho,phi);
  ddt(rho)-=vE_Grad(rho,phi0);
//  ddt(rho)-=vE_Grad(rho0,phi0); 

  ddt(rho)-=(mesh->g11)*sqrt(mesh->g11)*sqrt(mesh->g33)*DDX(Ni0)*(DDZ(phi)*D2DX2(phi0)-DDX(phi0)*D2DXDZ(phi));

  ddt(rho)-=(mesh->g33)*sqrt(mesh->g33)*sqrt(mesh->g11)*DDZ(Ni)*(DDZ(phi)*D2DXDZ(phi)-DDX(phi0)*D2DZ2(phi));
//  ddt(rho)+=Div(V_dot_Grad(vE,Grad(phi)));

  return 0;
}


/*// Performs inversion of rho (r) to get phi (p)
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
*/
/*
const Field2D vE_Grad(const Field2D &f, const Field2D &p)
{
  Field2D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f);
  }
  return result;
}

const Field3D vE_Grad(const Field2D &f, const Field3D &p)
{
  Field3D result;
  if(arakawa) {
    // Arakawa scheme for perpendicular flow. Here as a test
    
    result.allocate();
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(p)*DDX(f) - DDX(p)*DDZ(f)
          BoutReal Jpp = 0.25*( (p[jx][jy][jzp] - p[jx][jy][jzm])*
                                (f[jx+1][jy] - f[jx-1][jy]) -
                                (p[jx+1][jy][jz] - p[jx-1][jy][jz])*
                                (f[jx][jy] - f[jx][jy]) )
            / (mesh->dx[jx][jy] * mesh->dz);

          // J+x
          BoutReal Jpx = 0.25*( f[jx+1][jy]*(p[jx+1][jy][jzp]-p[jx+1][jy][jzm]) -
                                f[jx-1][jy]*(p[jx-1][jy][jzp]-p[jx-1][jy][jzm]) -
                                f[jx][jy]*(p[jx+1][jy][jzp]-p[jx-1][jy][jzp]) +
                                f[jx][jy]*(p[jx+1][jy][jzm]-p[jx-1][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          // Jx+
          BoutReal Jxp = 0.25*( f[jx+1][jy]*(p[jx][jy][jzp]-p[jx+1][jy][jz]) -
                                f[jx-1][jy]*(p[jx-1][jy][jz]-p[jx][jy][jzm]) -
                                f[jx-1][jy]*(p[jx][jy][jzp]-p[jx-1][jy][jz]) +
                                f[jx+1][jy]*(p[jx+1][jy][jz]-p[jx][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          
          result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
        }
    
  }else if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f);
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field2D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f);
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field3D &p)
{
  Field3D result;
  if(arakawa) {
    // Arakawa scheme for perpendicular flow. Here as a test
    
    result.allocate();
    
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(p)*DDX(f) - DDX(p)*DDZ(f)
          BoutReal Jpp = 0.25*( (p[jx][jy][jzp] - p[jx][jy][jzm])*
                                (f[jx+1][jy][jz] - f[jx-1][jy][jz]) -
                                (p[jx+1][jy][jz] - p[jx-1][jy][jz])*
                                (f[jx][jy][jzp] - f[jx][jy][jzm]) )
            / (mesh->dx[jx][jy] * mesh->dz);

          // J+x
          BoutReal Jpx = 0.25*( f[jx+1][jy][jz]*(p[jx+1][jy][jzp]-p[jx+1][jy][jzm]) -
                                f[jx-1][jy][jz]*(p[jx-1][jy][jzp]-p[jx-1][jy][jzm]) -
                                f[jx][jy][jzp]*(p[jx+1][jy][jzp]-p[jx-1][jy][jzp]) +
                                f[jx][jy][jzm]*(p[jx+1][jy][jzm]-p[jx-1][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          // Jx+
          BoutReal Jxp = 0.25*( f[jx+1][jy][jzp]*(p[jx][jy][jzp]-p[jx+1][jy][jz]) -
                                f[jx-1][jy][jzm]*(p[jx-1][jy][jz]-p[jx][jy][jzm]) -
                                f[jx-1][jy][jzp]*(p[jx][jy][jzp]-p[jx-1][jy][jz]) +
                                f[jx+1][jy][jzm]*(p[jx+1][jy][jz]-p[jx][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          
          result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
        }
    
  }else if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f) + VDDZ(-DDX(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f);
  }
  return result;
}
*/
