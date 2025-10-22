/*
* Terms in energy budget
* -------------------
* Calculate the energy conversion term following Kang (2010) 
* 
*/ 
#include "energy.h"
// #include "phys.h"
// #include "grid.h"
// #include "initialization.h"
// #include "state.h"

// Local function
static void VelocityDecomposition(energyT *energy, gridT *grid, physT *phys);
static void DensityDecomposition(energyT *energy, gridT *grid, physT *phys, propT *prop);
static void qGradient(energyT *energy, gridT *grid, physT *phys);
static void KineticEnergy(energyT *energy, gridT *grid, physT *phys);
static void PotentialEnergy(energyT *energy, gridT *grid, physT *phys, propT *prop);

/*
 * Function: AllocateEnergyVariables()
 * ------------------------------------
 * Allocate memory to the energy variable arrays
 *
 */
void AllocateEnergyVariables(gridT *grid, energyT **energy, propT *prop){
  int i, k;

  // allocate energy structure
  *energy = (energyT *)SunMalloc(sizeof(energyT),"AllocateEnergyVariables");

  // Allocate 2D arrays
  (*energy)->Uc = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->Vc = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->depth = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ek0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ep0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  // Allocate 3D arrays
  (*energy)->uc_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->vc_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->T_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->S_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->p_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->rho_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->rho_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->dqdz = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ek_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->Ek0_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ep_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<grid->Nc;i++) {
    (*energy)->uc_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->vc_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->T_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->S_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->p_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->rho_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->rho_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->dqdz[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->Ek_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->Ek0_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->Ep_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
  }

}

/*
 * Function: ZeroEnergyVariables() 
 * ---------------------------------
 * Zero the energy arrays
 */
void ZeroEnergyVariables(gridT *grid, energyT *energy, propT *prop ,MPI_Comm comm){
  int i,k;

  for(i=0;i<grid->Nc;i++) {
    energy->Uc[i]=0;
    energy->Vc[i]=0;
    energy->depth[i]=0;
    energy->Ek0[i]=0;
    energy->Ep0[i]=0;
    for(k=0;k<grid->Nk[i];k++) {
      energy->uc_prime[i][k] = 0;
      energy->vc_prime[i][k] = 0;
      energy->T_initial[i][k] = 0;
      energy->S_initial[i][k] = 0;
      energy->p_initial[i][k] = 0;
      energy->rho_initial[i][k] = 0;
      energy->rho_prime[i][k] = 0;
      energy->dqdz[i][k] = 0;
      energy->Ek_prime[i][k] = 0;
      energy->Ek0_prime[i][k] = 0;
      energy->Ep_prime[i][k] = 0;
    }
  }

}

/*
 * Function: FreeEnergyVariables()
 * ---------------------------------------------
 * This function frees all space allocated in AllocateEnergyVariables
 *
 */
void FreeEnergyVariables(gridT *grid, energyT *energy, propT *prop){
  int i;

  // free all the arrays over depth for cell-oriented
  for(i=0;i<grid->Nc;i++) {
    free(energy->uc_prime[i]);
    free(energy->vc_prime[i]);
    free(energy->T_initial[i]);
    free(energy->S_initial[i]);
    free(energy->p_initial[i]);
    free(energy->rho_initial[i]);
    free(energy->rho_prime[i]);
    free(energy->dqdz[i]);
    free(energy->Ek_prime[i]);
    free(energy->Ek0_prime[i]);
    free(energy->Ep_prime[i]);
  }
  free(energy->Uc);
  free(energy->Vc);
  free(energy->depth);
  free(energy->Ek0);
  free(energy->Ep0);
  free(energy->uc_prime);
  free(energy->vc_prime);
  free(energy->T_initial);
  free(energy->S_initial);
  free(energy->p_initial);
  free(energy->rho_initial);
  free(energy->rho_prime);
  free(energy->dqdz);  
  free(energy->Ek_prime);
  free(energy->Ek0_prime);
  free(energy->Ep_prime);

  free(energy);
}

/*
 * Function: VelocityDecomposition()
 * ---------------------------------
 * Computes the barotropic and baroclinic velocities Uc, Vc and uc_prime, vc_prime
 * Note: The integration for U and V calculation did not include eta, need to fix this! (S.Tan, 05/03/2023)
 */
static void VelocityDecomposition(energyT *energy, gridT *grid, physT *phys){
  int i, k;
  REAL ddz, Udz, Vdz;

    for(i=0;i<grid->Nc;i++) {
      ddz = 0;
      Udz = 0;
      Vdz = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
        ddz += grid->dzz[i][k];
        Udz += phys->uc[i][k]*grid->dzz[i][k];
        Vdz += phys->vc[i][k]*grid->dzz[i][k];
      }
      energy->depth[i] = ddz;
      energy->Uc[i] = Udz/ddz;
      energy->Vc[i] = Vdz/ddz;
    }

    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i];k++){//for(k=0;k<grid->Nk[i]+1;k++){
        energy->uc_prime[i][k] = phys->uc[i][k]-energy->Uc[i];
        energy->vc_prime[i][k] = phys->vc[i][k]-energy->Vc[i];
      }
    }
 
}

/*
 * Function: DensityDecomposition()
 * ---------------------------------
 * Decompose density rho(x,y,z,t): total density into
- rho0: constant reference density
- rho_b(z): background density
- rho_prime(x,y,z,t): deviation, may be interpreted as perturbation density due to wave motions
 * 
 */
static void DensityDecomposition(energyT *energy, gridT *grid, physT *phys, propT *prop){
  int i, k;
  REAL z;

  for(i=0;i<grid->Nc;i++) {
      z = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        z-=grid->dz[k]/2;
        energy->T_initial[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);
        energy->S_initial[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z,prop);
        energy->p_initial[i][k]=RHO0*prop->grav*z;
        z-=grid->dz[k]/2;
      }
  }

  for(i=0;i<grid->Nc;i++) {
    // for(k=grid->ctop[i];k<grid->Nk[i];k++) {
    for(k=0;k<grid->Nk[i];k++){
      // initial density
      energy->rho_initial[i][k]=StateEquation(prop,energy->S_initial[i][k],energy->T_initial[i][k],energy->p_initial[i][k]);
      // compute pertubation density rho_prime by extracting phys->rho with initial density
      energy->rho_prime[i][k]=phys->rho[i][k]-energy->rho_initial[i][k];
    }
  }

}

/*
 * Function: qGradient()
 * ---------------------------------
 * Compute vertical gradient for q: dqdz
 * 
 */
static void qGradient(energyT *energy, gridT *grid, physT *phys){
  int i, k;

  // compute dq/dz, store gradients at same location of W: top and bottom of each cell
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      energy->dqdz[i][k]=2.0*(phys->q[i][k-1]-phys->q[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
    }
    energy->dqdz[i][grid->ctop[i]]=energy->dqdz[i][grid->ctop[i]+1];
    // energy->dqdz[i][grid->Nk[i]]=energy->dqdz[i][grid->Nk[i]-1];
  }
 
}


/*
 * Function: KineticEnergy()
 * ---------------------------------
 * Computes the barotropic and baroclinic kinetic energy Ek0, Ek_prime, and cross term Ek0_prime
 * 
 */
static void KineticEnergy(energyT *energy, gridT *grid, physT *phys) {
  int i, k;

  for(i=0;i<grid->Nc;i++) {
    energy->Ek0[i]=RHO0*(pow(energy->Uc[i],2)+pow(energy->Vc[i],2))/2;
    for(k=0;k<grid->Nk[i];k++){
      energy->Ek_prime[i][k]=RHO0*(pow(energy->uc_prime[i][k],2)+pow(energy->vc_prime[i][k],2))/2;
      energy->Ek0_prime[i][k]=RHO0*(energy->Uc[i]*energy->uc_prime[i][k]+energy->Vc[i]*energy->vc_prime[i][k])/2;
    }
  }     
 
}

/*
 * Function: PotentialEnergy()
 * ---------------------------------
 * Computes the perturbation potential energy due to surface elevation Ep0, available potential energy (APE) Ep_prime
 *
 */
static void PotentialEnergy(energyT *energy, gridT *grid, physT *phys, propT *prop) {
  int i, k;

  for(i=0;i<grid->Nc;i++) {
    energy->Ep0[i]=RHO0*prop->grav*pow(phys->h[i],2)/2;
    for(k=0;k<grid->Nk[i];k++){
      energy->Ep_prime[i][k]=prop->grav;
    }
  }     
 
}

/*
 * Function: BarotropicW()
 * Usage: BarotropicW(phys->W,grid,phys,comm);
 * -----------------------------------------------------------------------------------------------
 * Computes the barotropic W following (7) in Kang and Fringer (2012)
 *
 */
 void BarotropicW(REAL **W, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, iptr, k, nf, ne;
  REAL UH_face, depth_face, height_face;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i]+1;k++){
      W[i][k] = 0;
    }
      
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
          W[i][k]-=height_face*UH_face*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/depth_face/grid->Ac[i];
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
void Conversion(REAL *C1, REAL *C2, REAL **W, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, j, k;
  REAL z, dC1, dC2;

  // DensityDecomposition(energy, grid, phys, prop);
  qGradient(energy, grid, phys);

  // compute depth integrated conversion term
  for(i=0;i<grid->Nc;i++) {
    C1[i] = 0;
    C2[i] = 0;
    // for(k=grid->ctop[i]+1;k<grid->Nk[i]+1;k++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      C1[i]+=RHO0*energy->rho_prime[i][k]*phys->W[i][k]*prop->grav*grid->dzz[i][k]; 
      C2[i]-=RHO0*energy->dqdz[i][k]*phys->W[i][k]*grid->dzz[i][k]; 
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
void EnergyFlux(REAL *Fx_0, REAL *Fy_0, REAL *Fx_prime, REAL *Fy_prime, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, j, k;

  VelocityDecomposition(energy, grid, phys);
  KineticEnergy(energy, grid, phys);
  PotentialEnergy(energy, grid, phys, prop);

  for(i=0;i<grid->Nc;i++) {
      Fx_0[i] = energy->Uc[i];
      Fy_0[i] = energy->depth[i];
      Fx_prime[i] = energy->uc_prime[i][0];
      Fy_prime[i] = energy->Ep0[i];
  }

}