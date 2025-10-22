/*
 * File: phys.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for phys.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _phys_h
#define _phys_h

#include "suntans.h"
#include "grid.h"
#include "fileio.h"

/*
 * Enumerated type definitions
 *
 */

// effectively flags for use in interpolation schemes
typedef enum _interpolation {
  nRT1, nRT2,
  tRT1, tRT2,
  QUAD, PEROT, LSQ
} interpolation;

/*
 * Main physical variable struct.
 *
 */
typedef struct _physT {
  REAL **u;
  REAL **uc;
  REAL **vc;
  REAL **wc, **wf;

  /*  new variables for nodal and tangential velocities */
  // definitions follow from Wang et al 2011
  // nRT1[Np][Nk][numpcneighs] so that for each node there is a 
  // value for each cell neighbor (non-unique to a node) and 
  // note that the number of cell neighbors varies based on node
  REAL ***nRT1u;
  REAL ***nRT1v;
  // nRT2[Np][Nk] has a unique value for each node since it is 
  // the area-weigted average of all the nRT1 values 
  // around the node
  REAL **nRT2u;
  REAL **nRT2v;
  // tRT1[Ne][Nk] is the area-weighted average of the nRT1 values
  // for the neighboring cells of an edge
  REAL **tRT1;
  // tRT2[Ne][Nk] is the area-weighted average of the nRT2 values
  // for each node at the end of the edge
  REAL **tRT2;

  REAL **uold;
  REAL **vold;
  REAL *D;
  REAL **w;
  REAL **q;
  REAL **qc;
  REAL **s;
  REAL **T;
  REAL **s0;
  REAL **rho;
  REAL *h;
  REAL *hcorr;
  unsigned char *active;

  REAL **boundary_u;
  REAL **boundary_v;
  REAL **boundary_w;
  REAL **boundary_s;
  REAL **boundary_T;
  REAL **boundary_tmp;
  REAL **boundary_rho;
  REAL *boundary_h;
  REAL *boundary_flag;

  REAL **nu_tv;
  REAL **kappa_tv;
  REAL **nu_lax;
  REAL *tau_T;
  REAL *tau_B;
  REAL *CdT;
  REAL *CdB;
  REAL *z0B;
  REAL *z0T;
  REAL **qT;
  REAL **lT;

  REAL mass;
  REAL mass0;
  REAL volume;
  REAL volume0;
  REAL Ep;
  REAL Ep0;
  REAL Ek;
  REAL Eflux1;
  REAL Eflux2;
  REAL Eflux3;
  REAL Eflux4;
  REAL smin;
  REAL smax;

  REAL *htmp;
  REAL *hold;
  REAL *htmp2;
  REAL *htmp3;
  REAL *hcoef;
  REAL *hfcoef;
  REAL **stmp;
  REAL **stmp2;
  REAL **stmp3;
  REAL **utmp;
  REAL **utmp2;
  REAL **ut;
  REAL **Cn_R;
  REAL **Cn_T;
  REAL **Cn_U;
  REAL **Cn_U2; //AB3
  REAL **Cn_W;
  REAL **Cn_W2; //AB3
  REAL **Cn_q;
  REAL **Cn_l;
  REAL **wnew;
  REAL **wtmp;
  REAL **wtmp2;
  REAL **qtmp;
  REAL **user_def_nc;
  REAL **user_def_ne;
  //GLS Turbulence variables
  REAL **TP;
  REAL **TB;
  REAL **TD;

  REAL *ap;
  REAL *am;
  REAL *bp;
  REAL *bm;

  REAL *a;
  REAL *b;
  REAL *c;
  REAL *d;

  // Horizontal facial scalar
  REAL **SfHp;
  REAL **SfHm;

  //Define variables for TVD schemes
  REAL *Cp;
  REAL *Cm;
  REAL *rp;
  REAL *rm;
  REAL *wp;
  REAL *wm;
  REAL **gradSx;
  REAL **gradSy; 
  
  //Variables for heat flux model
  REAL *Tsurf;
  REAL *dT;
  REAL **Ttmp;

  // Variables for netcdf write
  REAL *tmpvar;
  REAL *tmpvarW;
  REAL *tmpvarE;
  REAL *nctemp;

  // Least squares fit arrays
  REAL **A;
  REAL **Apr;
  REAL **AT;
  REAL *bpr;
  
  // pressure
  REAL **p_prime;  
  REAL **p_hyd;  

  // energy temporal integration
  REAL *C1tmp;  
  REAL *C2tmp;  

  REAL **C1tmp_;  
  REAL **C2tmp_;  

  REAL *Dtmp1;  
  REAL *Dtmp2; 

  REAL *F1tmp;  
  REAL *F2tmp;   
  REAL *F3tmp;  
  REAL *F4tmp; 

  REAL *F1tmp1;  
  REAL *F1tmp2;  
  REAL *F1tmp3;  
  REAL *F1tmp4;  

  REAL *F2tmp1;   
  REAL *F2tmp2;   
  REAL *F2tmp3;   
  REAL *F2tmp4;   

  REAL *F3tmp1;  
  REAL *F3tmp2;  
  REAL *F3tmp3;  
  REAL *F3tmp4;  
  REAL *F3tmp5;  

  REAL *F4tmp1; 
  REAL *F4tmp2; 
  REAL *F4tmp3; 
  REAL *F4tmp4; 
  REAL *F4tmp5; 

  // energy
  REAL *dEk0;
  REAL *dEp0;
  REAL *dEk_prime;
  REAL *dEp_prime;
  REAL *C1;  
  REAL *C2;  
  REAL *C1_int;  
  REAL *C2_int;  
  REAL **C1_;  
  REAL **C2_;  
  REAL **C1_int_;  
  REAL **C2_int_;  
  REAL *D_0;  
  REAL *D_prime;  
  REAL *D_0_int;  
  REAL *D_prime_int;  

  REAL *Fx_0;
  REAL *Fy_0;
  REAL *Fx_prime;
  REAL *Fy_prime;
  REAL *Fx_0_int;
  REAL *Fy_0_int;
  REAL *Fx_prime_int;
  REAL *Fy_prime_int;

  REAL *Fx_01;
  REAL *Fx_02;
  REAL *Fx_03;
  REAL *Fx_04;

  REAL *Fy_01;
  REAL *Fy_02;
  REAL *Fy_03;
  REAL *Fy_04;

  REAL *Fx_prime1;
  REAL *Fx_prime2;
  REAL *Fx_prime3;
  REAL *Fx_prime4;
  REAL *Fx_prime5;

  REAL *Fy_prime1;
  REAL *Fy_prime2;
  REAL *Fy_prime3;
  REAL *Fy_prime4;
  REAL *Fy_prime5;

  REAL *Fx_01_int;
  REAL *Fx_02_int;
  REAL *Fx_03_int;
  REAL *Fx_04_int;

  REAL *Fy_01_int;
  REAL *Fy_02_int;
  REAL *Fy_03_int;
  REAL *Fy_04_int;

  REAL *Fx_prime1_int;
  REAL *Fx_prime2_int;
  REAL *Fx_prime3_int;
  REAL *Fx_prime4_int;
  REAL *Fx_prime5_int;

  REAL *Fy_prime1_int;
  REAL *Fy_prime2_int;
  REAL *Fy_prime3_int;
  REAL *Fy_prime4_int;
  REAL *Fy_prime5_int;

} physT;

/*
 * Main property struct.
 *
 */
typedef struct _propT {
  REAL dt, Cmax, rtime, amp, omega, flux, timescale, theta0, theta, thetaM, 
       thetaS, thetaB, nu, nu_H, tau_T, z0T, CdT, z0B, CdB, CdW, relax, epsilon, qepsilon, resnorm, 
       dzsmall, beta, kappa_s, kappa_sH, gamma, kappa_T, kappa_TH, grav, Coriolis_f, CmaxU, CmaxW, 
       laxWendroff_Vertical, latitude;
  int ntout, ntoutStore, ntprog, nsteps, nbc, nstart, n, ntconserve, nonhydrostatic, cgsolver, maxiters, 
      qmaxiters, hprecond, qprecond, volcheck, masscheck, nonlinear, linearFS, newcells, wetdry, ramp_on, sponge_distance, 
    sponge_decay, thetaramptime, readSalinity, readTemperature, turbmodel, 
    TVD, horiTVD, vertTVD, TVDsalt, TVDtemp, TVDturb, laxWendroff, stairstep, AB, TVDmomentum, conserveMomentum,
    mergeArrays, computeSediments, subgrid, Intz0B, Intz0T, computeEnergy,C_zdependent, bottomPressure, simpleout, energyfluxCombine, nsolitons;
  int culvertmodel, marshmodel,wavemodel;
  FILE *FreeSurfaceFID, *HorizontalVelocityFID, *VerticalVelocityFID, *SalinityFID, *BGSalinityFID, 
       *InitSalinityFID, *InitTemperatureFID, *TemperatureFID, *PressureFID, *VerticalGridFID, *ConserveFID,    
       *StoreFID, *StartFID, *EddyViscosityFID, *ScalarDiffusivityFID, *EnergyFID; 
  interpolation interp; int prettyplot;
  int metmodel,  varmodel, outputNetcdf,  metncid, netcdfBdy, netcdfBdyFileID, readinitialnc, initialNCfileID, calcage, agemethod, calcaverage;
  int outputNetcdfFileID, averageNetcdfFileID, initialUNC, restartNC, restartAvgNC;
  REAL nctime, toffSet, gmtoffset;
  int nctimectr, avgtimectr, avgctr, avgfilectr, ntaverage, nstepsperncfile, ncfilectr;
  int sparsetimectr, sparsefilectr, sparseNetcdfFileID, outputNetcdfSparse, ntsparse, nksparse;
  REAL nugget, sill, range, Lsw, Cda, Ce, Ch;
  int wave_nesting, lowfreq_nudging, calcreynolds, skipAvgOutput, skipOutput;
  REAL TauL, TM2, ULm, ULz, UHtide, Uiw, VLm, VLz, Phitide, Phiiw, dW, Fhat_x, Fhat_y, 
       alpha1, alpha2, delta, drho, h2, D, lambda, xmid, ymid, interp_sponge, dtsoliton, bc_dt;
  char  starttime[15], basetime[15]; 
  char  INPUTZ0BFILE[BUFFERLENGTH], INPUTZ0TFILE[BUFFERLENGTH];
} propT;


/*
 * Public function declarations.
 *
 */
void Solve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void AllocatePhysicalVariables(gridT *grid, physT **phys, propT *prop);
void FreePhysicalVariables(gridT *grid, physT *phys, propT *prop);
void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void InitializeVerticalGrid(gridT **grid,int myproc);
void ReadProperties(propT **prop, gridT *grid, int myproc);
void SetDragCoefficients(gridT *grid, physT *phys, propT *prop);
REAL DepthFromDZ(gridT *grid, physT *phys, int i, int kind);
REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid);
void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc, interpolation interp, int subgridmodel) ;
void UpdateDZ(gridT *grid, physT *phys, propT *prop, int option);
void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void SetDensity(gridT *grid, physT *phys, propT *prop);

#endif
