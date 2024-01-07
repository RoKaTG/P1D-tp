/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include <math.h>

void eig_poisson1D(double* eigval, int *la) {
    int i;
    double h = 1.0 / (*la + 1);

    for (i = 0; i < *la; i++) {
        eigval[i] = 2.0 * (1 - cos((i + 1) * M_PI * h));
    }
}


double eigmax_poisson1D(int *la) {
    double h = 1.0 / (*la + 1);
    return 2.0 * (1 - cos(*la * M_PI * h));
}


double eigmin_poisson1D(int *la) {
    double h = 1.0 / (*la + 1);
    return 2.0 * (1 - cos(M_PI * h));
}


double richardson_alpha_opt(int *la) {
    double lambda_min = eigmin_poisson1D(la);
    double lambda_max = eigmax_poisson1D(la);
    return 2.0 / (lambda_min + lambda_max);
}


void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    double *AX = (double*) malloc(*la * sizeof(double)); // Pour stocker A*X
    double *Residual = (double*) malloc(*la * sizeof(double)); // Pour stocker le résidu
    int iter;
    double norm_res;

    for (int i = 0; i < *la; ++i) {
        X[i] = 0.0;
    }

    for (iter = 0; iter < *maxit; ++iter) {
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, AX, 1);

        for (int i = 0; i < *la; ++i) {
            Residual[i] = RHS[i] - AX[i];
        }

        norm_res = cblas_dnrm2(*la, Residual, 1);

        resvec[iter] = norm_res;

        if (norm_res < *tol) {
            break;
        }

        for (int i = 0; i < *la; ++i) {
            X[i] += *alpha_rich * Residual[i];
        }
    }

    *nbite = iter;

    free(AX);
    free(Residual);
}


void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
    for (int i = 0; i < *la; i++) {
        MB[i] = 1.0 / AB[*kv + i * (*lab)]; // Inverse des éléments de la diagonale
    }
}


void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv) {
    for (int j = 0; j < *la; j++) {
        for (int i = 0; i <= *kv; i++) {
            int index = i + j * (*lab);
            if (i <= j) {
                MB[index] = AB[index];
            } else {
                MB[index] = 0.0;
            }
        }
    }
}


void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    double *AX = (double*) malloc(*la * sizeof(double));
    double *Residual = (double*) malloc(*la * sizeof(double));
    double *Temp = (double*) malloc(*la * sizeof(double));
    int iter;
    double norm_res;

    for (int i = 0; i < *la; ++i) {
        X[i] = 0.0;
    }

    for (iter = 0; iter < *maxit; ++iter) {
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, AX, 1);

        for (int i = 0; i < *la; ++i) {
            Residual[i] = RHS[i] - AX[i];
        }

        for (int i = 0; i < *la; ++i) {
            Temp[i] = Residual[i] * MB[i];
        }

        norm_res = cblas_dnrm2(*la, Residual, 1);

        resvec[iter] = norm_res;

        if (norm_res < *tol) {
            break;
        }

        for (int i = 0; i < *la; ++i) {
            X[i] += Temp[i];
        }
    }

    *nbite = iter;

    free(AX);
    free(Residual);
    free(Temp);
}

void jacobi_GB(double *AB, double *RHS, double *X, int *lab, int *la, double *tol, int *maxit, double *resvec, int *nbite) {
    double *X_old = (double*) malloc(*la * sizeof(double));
    double *AX = (double*) malloc(*la * sizeof(double));
    int iter;
    double norm_res;

    for (int i = 0; i < *la; ++i) {
        X[i] = 0.0; 
    }

    for (iter = 0; iter < *maxit; ++iter) {
        for (int i = 0; i < *la; ++i) {
            X_old[i] = X[i];
        }

        // Calcul de AX
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *lab-1, *lab-1, 1.0, AB, *lab, X_old, 1, 0.0, AX, 1);

        for (int i = 0; i < *la; ++i) {
            X[i] = (RHS[i] - AX[i] + AB[*lab/2 + i * (*lab)] * X_old[i]) / AB[*lab/2 + i * (*lab)];
        }

        // Calcul du résidu et vérification de la convergence
        norm_res = 0;
        for (int i = 0; i < *la; ++i) {
            double temp = RHS[i] - AX[i];
            norm_res += temp * temp;
        }
        norm_res = sqrt(norm_res);

        resvec[iter] = norm_res;
        if (norm_res < *tol) {
            break;
        }
    }

    *nbite = iter;

    free(X_old);
    free(AX);
}

void gauss_seidel_GB(double *AB, double *RHS, double *X, int *lab, int *la, double *tol, int *maxit, double *resvec, int *nbite) {
    double *AX = (double*) malloc(*la * sizeof(double));
    int iter;
    double norm_res;

    for (int i = 0; i < *la; ++i) {
        X[i] = 0.0;
    }

    for (iter = 0; iter < *maxit; ++iter) {
        // Calcul de AX
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *lab-1, *lab-1, 1.0, AB, *lab, X, 1, 0.0, AX, 1);

        for (int i = 0; i < *la; ++i) {
            double sigma = 0.0;
            if (i > 0) {
                sigma += AB[*lab - 1 + i * (*lab)] * X[i - 1]; // Sous-diagonale
            }
            if (i < *la - 1) {
                sigma += AB[*lab + 1 + i * (*lab)] * X[i + 1]; // Sur-diagonale
            }
            X[i] = (RHS[i] - sigma) / AB[*lab/2 + i * (*lab)]; // Diagonale principale
        }

        // Calcul du résidu et vérification de la convergence
        norm_res = 0;
        for (int i = 0; i < *la; ++i) {
            double temp = RHS[i] - AX[i];
            norm_res += temp * temp;
        }
        norm_res = sqrt(norm_res);

        resvec[iter] = norm_res;
        if (norm_res < *tol) {
            break;
        }
    }

    *nbite = iter;

    free(AX);
}

void csr_poisson1D(double *AB, int *la, double *csr_values, int *csr_col_index, int *csr_row_ptr) {
    int n = *la; 
    int nnz = 3 * n - 2; 
    int row, col, idx = 0;

    csr_row_ptr[0] = 0;
    for (row = 0; row < n; ++row) {
        // Diagonale inférieure
        if (row > 0) {
            csr_values[idx] = -1; 
            csr_col_index[idx] = row - 1;
            idx++;
        }

        // Diagonale principale
        csr_values[idx] = 2; 
        csr_col_index[idx] = row;
        idx++;

        // Diagonale supérieure
        if (row < n - 1) {
            csr_values[idx] = -1; 
            csr_col_index[idx] = row + 1;
            idx++;
        }

        csr_row_ptr[row + 1] = idx;
    }
}

void csc_poisson1D(double *AB, int *la, double *csc_values, int *csc_row_index, int *csc_col_ptr) {
    int n = *la; 
    int nnz = 3 * n - 2;
    int row, col, idx = 0;

    csc_col_ptr[0] = 0;
    for (col = 0; col < n; ++col) {
        // Diagonale inférieure
        if (col > 0) {
            csc_values[idx] = -1;
            csc_row_index[idx] = col - 1;
            idx++;
        }

        // Diagonale principale
        csc_values[idx] = 2;
        csc_row_index[idx] = col;
        idx++;

        // Diagonale supérieure
        if (col < n - 1) {
            csc_values[idx] = -1;
            csc_row_index[idx] = col + 1;
            idx++;
        }

        csc_col_ptr[col + 1] = idx;
    }
}

void dcsrmv(int *la, double *csr_values, int *csr_col_index, int *csr_row_ptr, double *x, double *y) {
    int row;
    for (row = 0; row < *la; ++row) {
        y[row] = 0;
        for (int idx = csr_row_ptr[row]; idx < csr_row_ptr[row + 1]; ++idx) {
            y[row] += csr_values[idx] * x[csr_col_index[idx]];
        }
    }
}

void dcscmv(int *la, double *csc_values, int *csc_row_index, int *csc_col_ptr, double *x, double *y) {
    int col;
    for (col = 0; col < *la; ++col) {
        for (int idx = csc_col_ptr[col]; idx < csc_col_ptr[col + 1]; ++idx) {
            y[csc_row_index[idx]] += csc_values[idx] * x[col];
        }
    }
}
