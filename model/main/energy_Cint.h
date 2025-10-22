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
#include "util.h"

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

  // decomposed pressures
  REAL **p0;
  REAL **p_b;
  REAL **p_prime;

  // pressure gradients
  REAL **dqdz;

  // barotropic vertical velocity
  REAL **W;

  // energy terms
  REAL *Ek0;
  REAL *Ep0;
  REAL **Ek_prime;
  REAL **Ek0_prime;
  REAL **Ep_prime;

  // for testinterp1()
  REAL *xtmp;
  REAL *ytmp;
  REAL *xitmp;
  REAL *yitmp;  

} energyT;

/* *** Public Functions *** */
void AllocateEnergyVariables(gridT *grid, energyT **energy, propT *prop);
void ZeroEnergyVariables(gridT *grid, energyT *energy, propT *prop, MPI_Comm comm);
void Conversion(REAL *C1, REAL *C2, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);
void EnergyFlux(REAL *Fx_0, REAL *Fy_0, REAL *Fx_prime, REAL *Fy_prime, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);
void EnergyFluxDecompose(REAL *Fx_01, REAL *Fx_02, REAL *Fx_03, REAL *Fx_04, REAL *Fy_01, REAL *Fy_02, REAL *Fy_03, REAL *Fy_04, REAL *Fx_prime1, REAL *Fx_prime2, REAL *Fx_prime3, REAL *Fx_prime4, REAL *Fx_prime5, REAL *Fy_prime1, REAL *Fy_prime2, REAL *Fy_prime3, REAL *Fy_prime4, REAL *Fy_prime5, gridT *grid, energyT *energy, physT *phys, propT *prop, MPI_Comm comm);
void EnergyDiff(energyT *energy, gridT *grid, physT *phys, propT *prop, MPI_Comm comm);
void testinterp1(REAL *xout, REAL *yout, energyT *energy, gridT *grid, physT *phys, propT *prop);

#endif
