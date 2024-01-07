/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include <time.h> // Pour une estimation des performances

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;


  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

/************************************************/

  double beta;
  double alpha;

  int incx;
  int incy;

/***********************************************/
  
  beta = 0.0;
  alpha = 1.0;
  incx = 1;
  incy = 1;

/***********************************************/

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, X, incx, beta, RHS, incy);

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  clock_t start_TRF;
  clock_t end_TRF;

  clock_t start_TRI;
  clock_t end_TRI;

  clock_t start_TRS;
  clock_t end_TRS;

  clock_t start_SV;
  clock_t end_SV;


  /* LU Factorization */
  if (IMPLEM == TRF) {
    clock_t start_TRF = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    clock_t end_TRF = clock();
  }


  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    clock_t start_TRI = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    clock_t end_TRI = clock();
  }
  

  if (IMPLEM == TRI || IMPLEM == TRF){
    clock_t start_TRS = clock();
    /* Solution (Triangular) */
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
      printf("\n INFO = %d\n",info);
    }
    clock_t end_TRS = clock();

  }

  /* It can also be solved with dgbsv */
  // TODO : use dgbsv
  if (IMPLEM == SV) {
    clock_t start_SV = clock();
    if (info==0) {
        dgbsv_("N", &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info!=0){printf("\n INFO DGBSV = %d\n",info);}

    } else {
        printf("\n INFO = %d\n",info);
    }
    clock_t end_SV = clock();
  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);

  /*Benchmarking of every routines used to compute a linear equation*/
  double cpu_time_used_TRF = ((double) (end_TRF - start_TRF)) / CLOCKS_PER_SEC;
  double cpu_time_used_TRI = ((double) (end_TRI - start_TRI)) / CLOCKS_PER_SEC;
  double cpu_time_used_TRS = ((double) (end_TRS - start_TRS)) / CLOCKS_PER_SEC;
  double cpu_time_used_SV = ((double) (end_SV - start_SV)) / CLOCKS_PER_SEC;
  
  printf("\nThe relative forward error is relres = %e\n",relres);
  printf("\nTemps d'exécution pour dgbtrf: %f secondes\n", cpu_time_used_TRF);
  printf("Temps d'exécution pour dgbtrftridiag: %f secondes\n", cpu_time_used_TRI);
  printf("Temps d'exécution pour dgbtrs: %f secondes\n", cpu_time_used_TRS);
  printf("Temps d'exécution pour dgbsv: %f secondes\n", cpu_time_used_SV);


  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}

