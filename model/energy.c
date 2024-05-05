/*
* Terms in energy budget
* -------------------
* Calculate the energy conversion term following Kang (2010) 
* 
*/ 
#include "energy.h"
#include "initialization.h"
#include "state.h"

// Local function
static void VelocityDecomposition(energyT *energy, gridT *grid, physT *phys);
static void DensityDecomposition(energyT *energy, gridT *grid, physT *phys, propT *prop);
static void PressureDecomposition(energyT *energy, gridT *grid, physT *phys, propT *prop);
static void qGradient(energyT *energy, gridT *grid, physT *phys);
static void KineticEnergy(energyT *energy, gridT *grid, physT *phys);
static void PotentialEnergy(energyT *energy, gridT *grid, physT *phys, propT *prop);
static void BarotropicW(gridT *grid, energyT *energy, physT *phys, propT *prop);
static void swap(REAL xp, REAL yp);
void printArray(REAL* arr, int size);


/*
 * Function: AllocateEnergyVariables()
 * ------------------------------------
 * Allocate memory to the energy variable arrays
 *
 */
void AllocateEnergyVariables(gridT *grid, energyT **energy, propT *prop){
  int i, k;
  int Ntmp, Nkmax;

  // allocate energy structure
  *energy = (energyT *)SunMalloc(sizeof(energyT),"AllocateEnergyVariables");

  // Allocate 2D arrays
  (*energy)->Uc = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->Vc = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->depth = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ek0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ep0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  // Allocate 3D arrays
  (*energy)->W = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->uc_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->vc_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->T_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->S_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->p_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->rho_initial = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->rho_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->p0 = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->p_b = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->p_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->dqdz = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ek_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->Ek0_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  (*energy)->Ep_prime = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateEnergyVariables");

  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<grid->Nc;i++) {
    (*energy)->W[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->uc_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->vc_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->T_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->S_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->p_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->rho_initial[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->rho_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->p0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->p_b[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->p_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->dqdz[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->Ek_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
    (*energy)->Ek0_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");

    (*energy)->Ep_prime[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateEnergyVariables");
  }

  // for testinterp1()
  Ntmp = 2;
  Nkmax = 50;
  (*energy)->xtmp = (REAL *)SunMalloc(Ntmp*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->ytmp = (REAL *)SunMalloc(Ntmp*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->xitmp = (REAL *)SunMalloc(Nkmax*sizeof(REAL *),"AllocateEnergyVariables");
  (*energy)->yitmp = (REAL *)SunMalloc(Nkmax*sizeof(REAL *),"AllocateEnergyVariables");
}


/*
 * Function: ZeroEnergyVariables() 
 * ---------------------------------
 * Zero the energy arrays
 */
void ZeroEnergyVariables(gridT *grid, energyT *energy, propT *prop ,MPI_Comm comm){
  int i,k;
  int Ntmp, Nkmax;

  for(i=0;i<grid->Nc;i++) {
    energy->Uc[i]=0;
    energy->Vc[i]=0;
    energy->depth[i]=0;
    energy->Ek0[i]=0;
    energy->Ep0[i]=0;
    for(k=0;k<grid->Nk[i];k++) {
      energy->W[i][k] = 0;
      energy->uc_prime[i][k] = 0;
      energy->vc_prime[i][k] = 0;
      energy->T_initial[i][k] = 0;
      energy->S_initial[i][k] = 0;
      energy->p_initial[i][k] = 0;
      energy->rho_initial[i][k] = 0;
      energy->rho_prime[i][k] = 0;
      energy->p0[i][k] = 0;
      energy->p_b[i][k] = 0;
      energy->p_prime[i][k] = 0;
      energy->dqdz[i][k] = 0;
      energy->Ek_prime[i][k] = 0;
      energy->Ek0_prime[i][k] = 0;
      energy->Ep_prime[i][k] = 0;
    }
  }

  // for testinterp1()
  Ntmp = 2;
  Nkmax = 50;
  for(k=0;k<Ntmp;k++) {
    energy->xtmp[k]=0;
    energy->ytmp[k]=0;
  }
  for(k=0;k<Nkmax;k++) {
    energy->xitmp[k]=0;
    energy->yitmp[k]=0;
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
    free(energy->W[i]);
    free(energy->uc_prime[i]);
    free(energy->vc_prime[i]);
    free(energy->T_initial[i]);
    free(energy->S_initial[i]);
    free(energy->p_initial[i]);
    free(energy->rho_initial[i]);
    free(energy->rho_prime[i]);
    free(energy->p0[i]);
    free(energy->p_b[i]);
    free(energy->p_prime[i]);
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
  free(energy->W);
  free(energy->uc_prime);
  free(energy->vc_prime);
  free(energy->T_initial);
  free(energy->S_initial);
  free(energy->p_initial);
  free(energy->rho_initial);
  free(energy->rho_prime);
  free(energy->p0);
  free(energy->p_b);
  free(energy->p_prime);
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
 *
 * Note: The integration for U and V calculation did not include eta, need to fix this! (S.Tan, 05/03/2023), fixed (S.Tan, 11/06/2023)
 -- dzz include eta, so fixed when integrate from k=0 (dzz[0,:,2750]+eta[1,2750]=dzz[1,:,2750])
 */
static void VelocityDecomposition(energyT *energy, gridT *grid, physT *phys){
  int i, k;
  REAL ddz, Udz, Vdz;

    for(i=0;i<grid->Nc;i++) {
      ddz = 0;
      Udz = 0;
      Vdz = 0;
      for(k=0;k<grid->Nk[i];k++){
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
      z=phys->h[i];
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        z-=0.5*grid->dzz[i][k];
        energy->T_initial[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);
        energy->S_initial[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z,prop);
        energy->p_initial[i][k]=-RHO0*prop->grav*z;
        z-=0.5*grid->dzz[i][k];
      }
  }


  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      // initial density
      energy->rho_initial[i][k]=StateEquation(prop,energy->S_initial[i][k],energy->T_initial[i][k],energy->p_initial[i][k]);
      // compute pertubation density rho_prime by extracting phys->rho with initial density
      energy->rho_prime[i][k]=phys->rho[i][k]-energy->rho_initial[i][k];
    }
  }

}


/*
 * Function: PressureDecomposition()
 * ---------------------------------
 * Decompose pressure p: total pressure into p = ph + q: hydrostatic and non-hydrostatic pressures 
- ph: total hydrostatic pressure
- p0: constant reference pressure
- p_b: background pressure
- p_prime: deviation
 *
 * Note: The integration for p_b and p_prime calculation did not include eta, need to fix this! (S.Tan, 05/03/2023), fixed (S.Tan, 11/06/2023)
 * Now I'm fixing it, z[0] follows free surface, is not necessarily -dz[0]/2, z[-1] follows topo, not nesessary z[-2]-dz[-1]/2 (S.Tan, 11/06/2023)
 */ 

static void PressureDecomposition(energyT *energy, gridT *grid, physT *phys, propT *prop){
  int i, k;
  REAL z;

  for(i=0;i<grid->Nc;i++) {
    z=phys->h[i];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z-=0.5*grid->dzz[i][k];
      energy->p0[i][k]=RHO0*prop->grav*(phys->h[i]-z);
      if(k==grid->ctop[i]){
        energy->p_b[i][k]=RHO0*energy->rho_initial[i][k]*prop->grav*(phys->h[i]-z);
        energy->p_prime[i][k]=RHO0*energy->rho_prime[i][k]*prop->grav*(phys->h[i]-z);
      } else {
        energy->p_b[i][k]=energy->p_b[i][k-1]+RHO0*energy->rho_initial[i][k-1]*prop->grav*grid->dzz[i][k-1]/2+RHO0*energy->rho_initial[i][k]*prop->grav*grid->dzz[i][k]/2;
        energy->p_prime[i][k]=energy->p_prime[i][k-1]+RHO0*energy->rho_prime[i][k-1]*prop->grav*grid->dzz[i][k-1]/2+RHO0*energy->rho_prime[i][k]*prop->grav*grid->dzz[i][k]/2;
      }
      z-=0.5*grid->dzz[i][k];
    }
  }

}


/*
 * Function: qGradient()
 * ---------------------------------
 * Compute vertical gradient for q: dqdz
 * calculation inspired by WPredictor in phys.c
 */
static void qGradient(energyT *energy, gridT *grid, physT *phys){
  int i, k;

  // compute dq/dz, store gradients at same location of W: top and bottom of each cell, Nk+1
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      energy->dqdz[i][k]=2.0*(phys->q[i][k-1]-phys->q[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
    }
    energy->dqdz[i][grid->ctop[i]]=-2.0*phys->q[i][grid->ctop[i]]/grid->dzz[i][grid->ctop[i]];
  }
 
}


/*
 * Function: KineticEnergy()
 * ---------------------------------
 * Computes the barotropic and baroclinic kinetic energy Ek0, Ek_prime, and cross term Ek0_prime
 * Ek' now include w (S.Tan, 10/03/2023), note w is at top and bottom of each cell, Ek' is at center of cell
 */
static void KineticEnergy(energyT *energy, gridT *grid, physT *phys) {
  int i, k;

  for(i=0;i<grid->Nc;i++) {
    energy->Ek0[i]=RHO0*(pow(energy->Uc[i],2)+pow(energy->Vc[i],2))/2;
    for(k=0;k<grid->Nk[i];k++){
      energy->Ek_prime[i][k]=RHO0*(pow(energy->uc_prime[i][k],2)+pow(energy->vc_prime[i][k],2)+pow(0.5*(phys->w[i][k]+phys->w[i][k+1]),2))/2;
      energy->Ek0_prime[i][k]=RHO0*(energy->Uc[i]*energy->uc_prime[i][k]+energy->Vc[i]*energy->vc_prime[i][k])/2;
    }
  }     
 
}

// https://www.geeksforgeeks.org/c-program-to-sort-an-array-in-ascending-order/
// C program to sort the array in an
// ascending order using selection sort
static void swap(REAL xp, REAL yp)
{
    REAL temp = xp;
    xp = yp;
    yp = temp;
}


// Function to perform Selection Sort
static void selectionSort(REAL* arr, int n)
{
    int i, j, min_idx;
    REAL temp1, temp2;
  
    // One by one move boundary of
    // unsorted subarray
    for (i = 0; i < n - 1; i++) {
        // Find the minimum element in
        // unsorted array
        min_idx = i;
        for (j = i + 1; j < n; j++){
          if (arr[j] < arr[min_idx]){
            min_idx = j;
          }
        }
  
        // Swap the found minimum element
        // with the first element
        temp1 = arr[min_idx];
        temp2 = arr[i];
        arr[min_idx] = temp2;
        arr[i] = temp1;
    }
}


/*
 * Function: PotentialEnergy()
 * ---------------------------------
 * Computes the perturbation potential energy due to surface elevation Ep0, available potential energy (APE) Ep_prime
 *
 * Note: the current Ep_prime is the triangle (APE3) instead of arc-triangle (APE2) in Kang et al. (2010) (S.Tan somewhere spring, 2023)
 * in this new update, free surface is taken into consideration, and area of arc-triangle is calculated (S.Tan, 11/07/2023)
 * Note: now I use the Taylor expansion from Kang et al. (2010), note that at the deepest point, it calculate the area of triange as if there's rho_initial beneath it (S.Tan, 11/27/2023)
 */
static void PotentialEnergy(energyT *energy, gridT *grid, physT *phys, propT *prop) {
  int i, j, k, min_idx;
  REAL *z_i, *z_p, temp, *N2, *N2z; 
  REAL a1 = 23.36;
  REAL a2 = 3.13;
  REAL a3 = -44.12;
  REAL a4 = 293.12;

  N2 = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"ComputeZ");
  N2z= (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"ComputeZ");
  z_i = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"ComputeZ");
  z_p = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"ComputeZ");


  for(i=0;i<grid->Nc;i++) {
    energy->Ep0[i]=RHO0*prop->grav*pow(phys->h[i],2)/2; // this is depth integrated value!!!

    // compute N2, store gradients at same location of eho: center of each cell, Nk+1

    // for(i=0;i<grid->Nc;i++) {
    //   for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
    //     energy->N2[i][k]=-prop->grav/RHO0*(energy->rho_initial[i][k-1]+energy->rho_initial[i][k+1]-2*energy->rho_initial[i][k])/(grid->dzz[i][k]**2);
    //   }
    //   energy->N2[i][grid->ctop[i]]=-prop->grav/RHO0*(-2.0*energy->rho_initial[i][grid->ctop[i]]/grid->dzz[i][grid->ctop[i]]);
    //   energy->N2[i][grid->Nk[i]]=-prop->grav/RHO0*(-2.0*energy->rho_initial[i][grid->Nk[i]]/grid->dzz[i][grid->Nk[i]]);
    // }

    // initilize variables, fill with 0
    for(k=0;k<grid->Nkmax;k++){
      N2[k]=0;
      N2z[k]=0;
    }

    z_i[0] = -grid->dz[0]*0.5;
    for(k=1;k<grid->Nkmax;k++){
      z_i[k] = z_i[k-1]-grid->dz[k];
    }

    // z for perturbed rho profile, follow free surface and bottom
    z_p[grid->ctop[i]]=phys->h[i]-grid->dzz[i][grid->ctop[i]]/2;
    N2[grid->ctop[i]]=prop->grav*prop->gamma*a1*exp(-(-z_p[grid->ctop[i]]+a3)/a4)/a4;
    N2z[grid->ctop[i]]=prop->grav*prop->gamma*a1*exp(-(-z_p[grid->ctop[i]]+a3)/a4)/a4/a4;
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      z_p[k]=z_p[k-1]-(grid->dzz[i][k-1]+grid->dzz[i][k])/2;
      N2[k]=prop->grav*prop->gamma*a1*exp(-(-z_p[k]+a3)/a4)/a4;
      N2z[k]=prop->grav*prop->gamma*a1*exp(-(-z_p[k]+a3)/a4)/a4/a4;
    }

    // integration arc-triangle (APE2)
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      energy->Ep_prime[i][k]=pow(prop->grav,2)*pow((energy->rho_prime[i][k]*RHO0),2)/2/RHO0/N2[k]+pow(prop->grav,3)*N2z[k]*pow((energy->rho_prime[i][k]*RHO0),3)/6/pow(RHO0,2)/pow(N2[k],3);
    }
  }     
 
}


/*
 * Function: BarotropicW()
 * -----------------------------------------------------------------------------------------------
 * Computes the barotropic W following (7) in Kang and Fringer (2012), at surface and bottom of each cell, Nk+1
 *
 */
 static void BarotropicW(gridT *grid, energyT *energy, physT *phys, propT *prop) {
  int i, iptr, k, nf, ne;
  REAL UH_face, depth_face, height_face;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i]+1;k++){
      energy->W[i][k] = 0;
    }
    // no flux into bottom cells
    energy->W[i][grid->Nk[i]] = 0;
      
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
        // set the flux vertical area is explicit
        if (k<grid->Nke[ne] && height_face>0){
          energy->W[i][k]-=height_face*UH_face*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/depth_face/grid->Ac[i];
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
 * Note that W and dqdz are on surface and bottom of each cell, C is calculated at center of cell
 */
void Conversion(REAL *C1, REAL *C2, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, k;

  DensityDecomposition(energy, grid, phys, prop);
  qGradient(energy, grid, phys);
  BarotropicW(grid, energy, phys, prop);

  // compute depth integrated conversion term
  for(i=0;i<grid->Nc;i++) {
    C1[i] = 0;
    C2[i] = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      C1[i]+=0.5*RHO0*energy->rho_prime[i][k]*(energy->W[i][k]+energy->W[i][k+1])*prop->grav*grid->dzz[i][k]; 
      C2[i]-=0.5*RHO0*(energy->dqdz[i][k]*energy->W[i][k]+energy->dqdz[i][k+1]*energy->W[i][k+1])*grid->dzz[i][k]; 
    }
  }

}


/*
 * Function: EnergyFlux
 * Usage: EnergyFlux(phys->C1,phys->C2,phys->W,grid,phys,prop,comm,myproc);
 * -----------------------------------------------------------------------------------------------
 * Computes the Energy Flux term F_0 and F_prime following (15) and (16) in Kang and Fringer (2012)
 * Note that (see ../analysis/SUNTANS_island_energybudget.ipynb)
 * 1) The integration for depth-integrated energy flux calculation did not include eta, need to fix this! (S.Tan, 05/03/2023), fixed - integrate with dzz (S.Tan, 11/07/2023)
 * 2) TO UPDATE: Ek' does not include w (S.Tan, 05/03/2023) - resolved (S.Tan, 10/03/2023)
 *
 */
void EnergyFlux(REAL *Fx_0, REAL *Fy_0, REAL *Fx_prime, REAL *Fy_prime, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, j, k;
  REAL Fx_0_1, Fx_0_2, Fx_0_3, Fx_0_4;
  REAL Fy_0_1, Fy_0_2, Fy_0_3, Fy_0_4;
  REAL Fx_prime_1, Fx_prime_2, Fx_prime_3, Fx_prime_4, Fx_prime_5;
  REAL Fy_prime_1, Fy_prime_2, Fy_prime_3, Fy_prime_4, Fy_prime_5;

  VelocityDecomposition(energy, grid, phys);
  KineticEnergy(energy, grid, phys);
  DensityDecomposition(energy, grid, phys, prop);
  PotentialEnergy(energy, grid, phys, prop);
  PressureDecomposition(energy, grid, phys, prop);

  // barotropic terms
  for(i=0;i<grid->Nc;i++) {
    Fx_0_1 = 0; 
    Fx_0_2 = RHO0*prop->grav*energy->Uc[i]*energy->depth[i]*phys->h[i];
    Fx_0_3 = 0; 
    Fx_0_4 = 0; 
    Fy_0_1 = 0; 
    Fy_0_2 = RHO0*prop->grav*energy->Vc[i]*energy->depth[i]*phys->h[i];
    Fy_0_3 = 0; 
    Fy_0_4 = 0; 
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      Fx_0_1+=energy->Ek0[i]*grid->dzz[i][k]; 
      Fx_0_3+=energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fx_0_4+=RHO0*phys->q[i][k]*grid->dzz[i][k]; 
      Fy_0_1+=energy->Ek0[i]*grid->dzz[i][k]; 
      Fy_0_3+=energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fy_0_4+=RHO0*phys->q[i][k]*grid->dzz[i][k]; 
    }
   Fx_0_1*=energy->Uc[i];
   Fx_0_3*=energy->Uc[i];
   Fx_0_4*=energy->Uc[i];
   Fy_0_1*=energy->Vc[i];
   Fy_0_3*=energy->Vc[i];
   Fy_0_4*=energy->Vc[i];

   Fx_0[i] = Fx_0_1+Fx_0_2+Fx_0_3+Fx_0_4;
   Fy_0[i] = Fy_0_1+Fy_0_2+Fy_0_3+Fy_0_4;
  }

  // baroclinic terms
  for(i=0;i<grid->Nc;i++) {
    Fx_prime_1 = 0;
    Fx_prime_2 = 0; 
    Fx_prime_3 = 0; 
    Fx_prime_4 = 0; 
    Fx_prime_5 = 0;
    Fy_prime_1 = 0;
    Fy_prime_2 = 0; 
    Fy_prime_3 = 0; 
    Fy_prime_4 = 0; 
    Fy_prime_5 = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      Fx_prime_1+=phys->uc[i][k]*energy->Ek_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_2+=phys->uc[i][k]*energy->Ek0_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_3+=phys->uc[i][k]*energy->Ep_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_4+=energy->uc_prime[i][k]*energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_5+=energy->uc_prime[i][k]*RHO0*phys->q[i][k]*grid->dzz[i][k]; 
      Fy_prime_1+=phys->vc[i][k]*energy->Ek_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_2+=phys->vc[i][k]*energy->Ek0_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_3+=phys->vc[i][k]*energy->Ep_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_4+=energy->vc_prime[i][k]*energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_5+=energy->vc_prime[i][k]*RHO0*phys->q[i][k]*grid->dzz[i][k]; 
    }
   Fx_prime[i] = Fx_prime_1+Fx_prime_2+Fx_prime_3+Fx_prime_4+Fx_prime_5;
   Fy_prime[i] = Fy_prime_1+Fy_prime_2+Fy_prime_3+Fy_prime_4+Fy_prime_5;
  }

}


/*
 * Function: EnergyFluxDecompose
 * Usage: EnergyFluxDecompose(phys->C1,phys->C2,phys->W,grid,phys,prop,comm,myproc);
 * -----------------------------------------------------------------------------------------------
 * Computes the Energy Flux term F_0 and F_prime following (15) and (16) in Kang and Fringer (2012)
 * Note that (see ../analysis/SUNTANS_island_energybudget.ipynb)
 * 1) The integration for depth-integrated energy flux calculation did not include eta, need to fix this! (S.Tan, 05/03/2023)
 * 2) TO UPDATE: Ek' does not include w (S.Tan, 05/03/2023)
 *
 */
void EnergyFluxDecompose(REAL *Fx_01, REAL *Fx_02, REAL *Fx_03, REAL *Fx_04, REAL *Fy_01, REAL *Fy_02, REAL *Fy_03, REAL *Fy_04, REAL *Fx_prime1, REAL *Fx_prime2, REAL *Fx_prime3, REAL *Fx_prime4, REAL *Fx_prime5, REAL *Fy_prime1, REAL *Fy_prime2, REAL *Fy_prime3, REAL *Fy_prime4, REAL *Fy_prime5, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm) {
  int i, j, k;
  REAL Fx_0_1, Fx_0_2, Fx_0_3, Fx_0_4;
  REAL Fy_0_1, Fy_0_2, Fy_0_3, Fy_0_4;
  REAL Fx_prime_1, Fx_prime_2, Fx_prime_3, Fx_prime_4, Fx_prime_5;
  REAL Fy_prime_1, Fy_prime_2, Fy_prime_3, Fy_prime_4, Fy_prime_5;

  VelocityDecomposition(energy, grid, phys);
  KineticEnergy(energy, grid, phys);
  DensityDecomposition(energy, grid, phys, prop);
  PotentialEnergy(energy, grid, phys, prop);
  PressureDecomposition(energy, grid, phys, prop);

  // barotropic terms
  for(i=0;i<grid->Nc;i++) {
    Fx_0_1 = 0; 
    Fx_0_2 = RHO0*prop->grav*energy->Uc[i]*energy->depth[i]*phys->h[i];
    Fx_0_3 = 0; 
    Fx_0_4 = 0; 
    Fy_0_1 = 0; 
    Fy_0_2 = RHO0*prop->grav*energy->Vc[i]*energy->depth[i]*phys->h[i];
    Fy_0_3 = 0; 
    Fy_0_4 = 0; 
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      Fx_0_1+=energy->Ek0[i]*grid->dzz[i][k]; 
      Fx_0_3+=energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fx_0_4+=RHO0*phys->q[i][k]*grid->dzz[i][k]; 
      Fy_0_1+=energy->Ek0[i]*grid->dzz[i][k]; 
      Fy_0_3+=energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fy_0_4+=RHO0*phys->q[i][k]*grid->dzz[i][k]; 
    }
   Fx_0_1*=energy->Uc[i];
   Fx_0_3*=energy->Uc[i];
   Fx_0_4*=energy->Uc[i];
   Fy_0_1*=energy->Vc[i];
   Fy_0_3*=energy->Vc[i];
   Fy_0_4*=energy->Vc[i];

   Fx_01[i] = Fx_0_1;
   Fx_02[i] = Fx_0_2; 
   Fx_03[i] = Fx_0_3;
   Fx_04[i] = Fx_0_4; 
   Fy_01[i] = Fy_0_1;
   Fy_02[i] = Fy_0_2; 
   Fy_03[i] = Fy_0_3;
   Fy_04[i] = Fy_0_4; 
  }

  // baroclinic terms
  for(i=0;i<grid->Nc;i++) {
    Fx_prime_1 = 0;
    Fx_prime_2 = 0; 
    Fx_prime_3 = 0; 
    Fx_prime_4 = 0; 
    Fx_prime_5 = 0;
    Fy_prime_1 = 0;
    Fy_prime_2 = 0; 
    Fy_prime_3 = 0; 
    Fy_prime_4 = 0; 
    Fy_prime_5 = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      Fx_prime_1+=phys->uc[i][k]*energy->Ek_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_2+=phys->uc[i][k]*energy->Ek0_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_3+=phys->uc[i][k]*energy->Ep_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_4+=energy->uc_prime[i][k]*energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fx_prime_5+=energy->uc_prime[i][k]*RHO0*phys->q[i][k]*grid->dzz[i][k]; 
      Fy_prime_1+=phys->vc[i][k]*energy->Ek_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_2+=phys->vc[i][k]*energy->Ek0_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_3+=phys->vc[i][k]*energy->Ep_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_4+=energy->vc_prime[i][k]*energy->p_prime[i][k]*grid->dzz[i][k]; 
      Fy_prime_5+=energy->vc_prime[i][k]*RHO0*phys->q[i][k]*grid->dzz[i][k]; 
    }
   Fx_prime1[i] = Fx_prime_1;
   Fx_prime2[i] = Fx_prime_2;
   Fx_prime3[i] = Fx_prime_3;
   Fx_prime4[i] = Fx_prime_4;
   Fx_prime5[i] = Fx_prime_5;
   Fy_prime1[i] = Fy_prime_1;
   Fy_prime2[i] = Fy_prime_2;
   Fy_prime3[i] = Fy_prime_3;
   Fy_prime4[i] = Fy_prime_4;
   Fy_prime5[i] = Fy_prime_5;
  }

}


/*
 * Function: BottomDrag()
 * ---------------------------------
 * Compute bottm drag D_0 and D_prime from Cd and velocities at the bottom cell
 */
void BottomDrag(REAL *D_0, REAL *D_prime, energyT *energy, gridT *grid, physT *phys, propT *prop, MPI_Comm comm){
  int i, k;

  VelocityDecomposition(energy, grid, phys);

  for(i=0;i<grid->Nc;i++) {
    k=grid->Nk[i]-1; 
    D_0[i] = RHO0*prop->CdB*sqrt(pow(phys->uc[i][k],2)+pow(phys->vc[i][k],2))*(phys->uc[i][k]*energy->Uc[i]+phys->vc[i][k]*energy->Vc[i]);
    D_prime[i] = RHO0*prop->CdB*sqrt(pow(phys->uc[i][k],2)+pow(phys->vc[i][k],2))*(phys->uc[i][k]*energy->uc_prime[i][k]+phys->vc[i][k]*energy->vc_prime[i][k]+pow(0.5*(phys->w[i][k]+phys->w[i][k+1]),2));
  }

}


/*
 * Function: EnergyDiff()
 * ---------------------------------
 * Compute depth-integrated Ek0, Ep0, Ek_prime, Ep_prime
 */
void EnergyDiff(energyT *energy, gridT *grid, physT *phys, propT *prop, MPI_Comm comm){
  int i, k;
  REAL E1, E3, E4;

  for(i=0;i<grid->Nc;i++) {
    E1 = 0;
    E3 = 0;
    E4 = 0;
    phys->dEk0[i] = 0;
    phys->dEk_prime[i] = 0;
    phys->dEp_prime[i] = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      E1+=energy->Ek0[i]*grid->dzz[i][k]; 
      E3+=energy->Ek_prime[i][k]*grid->dzz[i][k]; 
      E4+=energy->Ep_prime[i][k]*grid->dzz[i][k]; 
    }
    phys->dEk0[i] = E1;
    phys->dEp0[i] = energy->Ep0[i]; // energy->Ep0 is already depth integration (Kang 2010 eq. 5.53)
    phys->dEk_prime[i] = E3;
    phys->dEp_prime[i] = E4;
  }

}


/*
 * Function: testinterp1()
 * ---------------------------------
 * this function is to test my interp1 in util.c, it works fine
 */
void testinterp1(REAL *xout, REAL *yout, energyT *energy, gridT *grid, physT *phys, propT *prop){
  int i, j, k, Ntmp, Nkmax;

  Ntmp = 2;
  Nkmax = 50;

  energy->xtmp[0]=-100;
  energy->xtmp[1]=100;
  energy->ytmp[0]=-10;
  energy->ytmp[1]=10;

  for(i=0;i<Nkmax;i++) {
    energy->xitmp[i]=0;
    energy->yitmp[i]=0;
  }
  for(i=1;i<Nkmax;i++) {
    energy->xitmp[i]=-i;
  }

  interp1(energy->xtmp, energy->ytmp, Ntmp, energy->xitmp, energy->yitmp, Nkmax);

  for(i=0;i<grid->Nc;i++) {
    if(i<Nkmax){
      xout[i]=energy->xitmp[i];
      yout[i]=energy->yitmp[i];
    }else{
      xout[i]=0;
      yout[i]=0;
    }
  }

}