/*
 * Energy header file
 *
 */

#ifndef _energy_h
#define _energy_h

#include "phys.h"
#include "grid.h"
#include "memory.h"
#include "initialization.h"
#include "state.h"
#include "met.h"

/*
 * Main energy variable struct.
 */
typedef struct _energyT {
  
  REAL *depth;

  // decomposed velocities
  REAL *Uc;
  REAL *Vc;
  REAL **uc_prime;
  REAL **vc_prime;

  // decomposed densities, & related
  REAL **T_initial;
  REAL **S_initial;
  REAL **p_initial;
  REAL **rho_initial;
  REAL **rho_prime;

  // pressure gradients
  REAL **dqdz;

  // energy terms
  REAL *Ek0;
  REAL *Ep0;
  REAL **Ek_prime;
  REAL **Ek0_prime;
  REAL **Ep_prime;

} energyT;

/* *** Public Functions *** */
void AllocateEnergyVariables(gridT *grid, energyT **energy, propT *prop);
void ZeroEnergyVariables(gridT *grid, energyT *energy, propT *prop, MPI_Comm comm);
void BarotropicW(REAL **W, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);
void Conversion(REAL *C1, REAL *C2, REAL **W, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);
void EnergyFlux(REAL *Fx_0, REAL *Fy_0, REAL *Fx_prime, REAL *Fy_prime, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);

#endif
