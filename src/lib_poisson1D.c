/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int i, j;
    for (i = 0; i < (*lab) * (*la); ++i) {
        AB[i] = 0.0;
    }

    for (int j = 0; j < *la; ++j) {
        AB[(*kv) + 1 + j * (*lab)] = 2.0;
        if (j < (*la) - 1) {
            AB[(*kv + 2) + j * (*lab)] = -1.0;
        }
        if (j > 0) {
            AB[(*kv) + j * (*lab)] = -1.0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

/* Correction du prof
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}
*/


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv) {
    int i, j;
    for (i = 0; i < (*lab) * (*la); i++) {
        AB[i] = 0.0;
    }
    for (j = 0; j < *la; j++) {
        AB[(*kv) + j * (*lab)] = 1.0;
    }
}

///////////////////////////////////////////////////////////////////////////////////

/* Correction du prof
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}
*/


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
      for (int i = 0; i < (*la); i++) {
        RHS[i] = 0.0;
        }
    RHS[0] = (*BC0);
    RHS[(*la) - 1] = (*BC1);
}

///////////////////////////////////////////////////////////////////////////

/* Correction du prof
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  
*/


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  //EX_SOL[(*la)] = EX_SOL[(*la) + 2];
  
  for (int i = 0; i < ((*la)); i++) {
    EX_SOL[i] = (*BC0) + X[i] * ((*BC1) - (*BC0));  
  } 
}

//////////////////////////////////////////////////////////////////////////////////////////////////

/* Correction du prof
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  
*/


void set_grid_points_1D(double* x, int* la){
    double dx = 1.0 / (*la + 1);
    int i;

    for (i = 0; i < *la; i++) {
        x[i] = (i + 1) * dx;
    }
}

////////////////////////////////////////////

/* Correction du prof 
void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}
*/


void LU_Facto(double* AB, int *lab, int *la, int *kv) {
    int n = *la;
    int i;
    if (*lab != 3) {
        printf("Erreur: La largeur de la bande n'est pas correcte pour une matrice tridiagonale.\n");
        return;
    }
    for (i = 0; i < n - 1; ++i) {
        if (AB[*kv + i * (*lab)] == 0) {
            printf("Erreur: Division par zéro lors de la factorisation LU.\n");
            return;
        }
        AB[(*kv + 1) + i * (*lab)] /= AB[*kv + i * (*lab)];
        AB[*kv + (i + 1) * (*lab)] -= AB[(*kv + 1) + i * (*lab)] * AB[(*kv - 1) + (i + 1) * (*lab)];
    }
}

double relative_forward_error(double* x, double* y, int* la){
  return 0;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
