/*
* Terms in energy budget
* -------------------
* Calculate the energy conversion term following Kang (2010) 
* 
*/ 
#include "energy.h"
#include "phys.h"
#include "grid.h"
#include "mynetcdf.h"
#include "initialization.h"
#include "state.h"

// Local function
static void VelocityDecomposition(REAL *Uc_tmp, REAL *Vc_tmp, REAL *depth_tmp, REAL **uc_prime_tmp, REAL **vc_prime_tmp, gridT *grid, physT *phys);
static void KineticEnergy(REAL *Ek0, REAL **Ek_prime, REAL **Ek0_prime, REAL *Uc_tmp, REAL *Vc_tmp, REAL **uc_prime_tmp, REAL **vc_prime_tmp, gridT *grid);

/*
 * Function: VelocityDecomposition
 * Usage: VelocityDecomposition(grid,phys);
 * ---------------------------------
 * Computes the barotropic and baroclinic velocities Uc, Vc and uc_prime, vc_prime
 *
 */
static void VelocityDecomposition(REAL *Uc_tmp, REAL *Vc_tmp, REAL *depth_tmp, REAL **uc_prime_tmp, REAL **vc_prime_tmp, gridT *grid, physT *phys){
  int i, j, k;

    // Uc_tmp=0;
    // Vc_tmp=0;
    // depth_tmp=0;

    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        depth_tmp[i] += grid->dzz[i][k];
        Uc_tmp[i] += phys->uc[i][k]*grid->dzz[i][k];
        Vc_tmp[i] += phys->vc[i][k]*grid->dzz[i][k];
    }

    for(i=0;i<grid->Nc;i++) {
      Uc_tmp[i] = Uc_tmp[i]/depth_tmp[i];
      Vc_tmp[i] = Vc_tmp[i]/depth_tmp[i];
    }

    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i]+1;k++){
        uc_prime_tmp[i][k] = phys->uc[i][k]-Uc_tmp[i];
        vc_prime_tmp[i][k] = phys->vc[i][k]-Vc_tmp[i];
      }
    }
 
}

/*
 * Function: KineticEnergy
 * Usage: KineticEnergy(grid,phys);
 * ---------------------------------
 * Computes the barotropic and baroclinic kinetic energy Ek0, Ek_prime, and cross term Ek0_prime
 *
 */
static void KineticEnergy(REAL *Ek0, REAL **Ek_prime, REAL **Ek0_prime, REAL *Uc_tmp, REAL *Vc_tmp, REAL **uc_prime_tmp, REAL **vc_prime_tmp, gridT *grid) {
  int i, k;
  // REAL *Uc_tmp, *Vc_tmp, **uc_prime_tmp, **vc_prime_tmp;

  // VelocityDecomposition(Uc_tmp, Vc_tmp, depth_tmp, uc_prime_tmp, vc_prime_tmp, grid, phys);
  for(i=0;i<grid->Nc;i++) {
    Ek0[i]=RHO0*(pow(Uc_tmp[i],Uc_tmp[i])+pow(Vc_tmp[i],Vc_tmp[i]))/2;
    for(k=0;k<grid->Nk[i]+1;k++){
      Ek_prime[i][k]=RHO0*(pow(uc_prime_tmp[i][k],uc_prime_tmp[i][k])+pow(vc_prime_tmp[i][k],vc_prime_tmp[i][k]))/2;
      Ek0_prime[i][k]=RHO0*(pow(Uc_tmp[i],uc_prime_tmp[i][k])+pow(Vc_tmp[i],vc_prime_tmp[i][k]))/2;
    }
  }     
 
}

/*
 * Function: BarotropicW
 * Usage: BarotropicW(phys->W,grid,phys,comm);
 * -----------------------------------------------------------------------------------------------
 * Computes the barotropic W following (7) in Kang and Fringer (2012)
 *
 */
 void BarotropicW(REAL **W, gridT *grid, physT *phys, MPI_Comm comm) {
  int i, iptr, k, nf, ne;
  REAL UH_face, depth_face, height_face;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i]+1;k++)
      W[i][k] = 0;

    // for each face
    for(nf=0;nf<grid->nfaces[i];nf++) {
      // get the edge pointer
      ne = grid->face[i*grid->maxfaces+nf];
      
      // compute total water depth for each face
      depth_face = 0;
      // compute integrated velocity for each face
      UH_face = 0;
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
        depth_face += grid->dzf[ne][k];
        UH_face += phys->u[ne][k]*grid->dzf[ne][k];
      }

      // water colume height z+d for each face
      height_face = 0;

      // compute W from the horizontal divergence of barotropic horizontal velocities
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
        height_face+=grid->dzf[ne][k]; 
        if (height_face>0){
          W[i][k]-=height_face*UH_face*
            grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/depth_face/grid->Ac[i];
        }
      }
    }

  }
}


/*
 * Function: Conversion
 * Usage: Conversion(phys->C1,phys->C2,phys->W,grid,phys,prop,comm,myproc);
 * -----------------------------------------------------------------------------------------------
 * Computes the conversion term C following (17) in Kang and Fringer (2012)
 *
 */
void Conversion(REAL *C1, REAL *C2, REAL **W, gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc) {
  int i, j, k;
  REAL z, **T_initial, **S_initial, **p_initial;
  REAL **rho_initial, **rho_prime;
  REAL **dqdz;


  T_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  S_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  p_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  rho_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  rho_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  dqdz = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");

  for(i=0;i<grid->Nc;i++) {
    T_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    S_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    p_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    rho_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    rho_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    dqdz[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  for(i=0;i<grid->Nc;i++) {
      z = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        z-=grid->dz[k]/2;
        T_initial[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);
        S_initial[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z,prop);
        p_initial[i][k]=RHO0*prop->grav*z;
        z-=grid->dz[k]/2;
      }
  }

  // initial density
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      rho_initial[i][k]=StateEquation(prop,S_initial[i][k],T_initial[i][k],p_initial[i][k]);
    }
  }

  // compute pertubation density rho_prime by extracting phys->rho with initial density
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      rho_prime[i][k] = phys->rho[i][k] - rho_initial[i][k];
    }
  }

  // compute dq/dz, store gradients at same location of W: top and bottom of each cell
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      dqdz[i][k]=2.0*(phys->q[i][k-1]-phys->q[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
    }
    dqdz[i][grid->ctop[i]]=dqdz[i][grid->ctop[i]+1];
    // dqdz[i][grid->Nk[i]]=dqdz[i][grid->Nk[i]-1];
  }


  // SunFree(T_initial,grid->Nkmax*grid->Nc*sizeof(REAL),"InitializePhyiscalVariables");

  // compute depth integrated conversion term
  for(i=0;i<grid->Nc;i++) {
    C1[i] = 0;
    C2[i] = 0;
    // for(k=grid->ctop[i]+1;k<grid->Nk[i]+1;k++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      C1[i]+=RHO0*rho_prime[i][k]*prop->grav*W[i][k]*grid->dzz[i][k]; 
      C2[i]-=RHO0*dqdz[i][k]*W[i][k]*grid->dzz[i][k];
    }
  }


}


/*
 * Function: EnergyFlux
 * Usage: EnergyFlux(phys->C1,phys->C2,phys->W,grid,phys,prop,comm,myproc);
 * -----------------------------------------------------------------------------------------------
 * Computes the Energy Flux term F_0 and F_prime following (15) and (16) in Kang and Fringer (2012)
 * Note that (see ../analysis/SUNTANS_island_energybudget.ipynb)
 * 1) The integration for depth-integrated energy flux calculation did not include eta, need to fix this! (S.Tan, 05/03/2023)
 * 2) TO UPDATE: Ek' does not include w (S.Tan, 05/03/2023)
 *
 */
void EnergyFlux(REAL *Fx_0, REAL *Fy_0, REAL *Fx_prime, REAL *Fy_prime, gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc) {
  int i, j, k;
  REAL *Uc_tmp, *Vc_tmp, *depth_tmp, **uc_prime_tmp, **vc_prime_tmp;
  REAL *Ek0, **Ek_prime, **Ek0_prime;

  Uc_tmp = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  Vc_tmp = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  depth_tmp = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  uc_prime_tmp = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  vc_prime_tmp = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");

  Ek0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  Ek_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  Ek0_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocatePhysicalVariables");

  for(i=0;i<grid->Nc;i++) {
    uc_prime_tmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    vc_prime_tmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    Ek_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    Ek0_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  for(i=0;i<grid->Nc;i++) {
    Uc_tmp[i] = 0;
    Vc_tmp[i] = 0;
    depth_tmp[i] = 0;
    Ek0[i] = 0;
    for(k=0;k<grid->Nk[i]+1;k++){
      uc_prime_tmp[i][k] = 0;
      vc_prime_tmp[i][k] = 0;
      Ek_prime[i][k] = 0;
      Ek0_prime[i][k] = 0;
    }
  }

  VelocityDecomposition(Uc_tmp, Vc_tmp, depth_tmp, uc_prime_tmp, vc_prime_tmp, grid, phys);
  KineticEnergy(Ek0, Ek_prime, Ek0_prime, Uc_tmp, Vc_tmp, uc_prime_tmp, vc_prime_tmp, grid);

  for(i=0;i<grid->Nc;i++) {
      Fx_0[i] = Uc_tmp[i];
      Fy_0[i] = Vc_tmp[i];
      Fx_prime[i] = uc_prime_tmp[i][0];
      Fy_prime[i] = vc_prime_tmp[i][0];
  }

}