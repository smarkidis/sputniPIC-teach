#ifndef GMRES_H
#define GMRES_H

#include <math.h>
#include "Basic.h"

typedef void (*GENERIC_IMAGE_GMRES) (FPfield *, FPfield *, EMfield *, interpDensSpecies *, grid *, parameters *);

void ApplyPlaneRotation(FPfield &dx, FPfield &dy, FPfield &cs, FPfield &sn);

void GMRes(GENERIC_IMAGE_GMRES FunctionImage, FPfield *xkrylov, int xkrylovlen, FPfield *b, int m, int max_iter, double tol, EMfield* field, interpDensSpecies* ids, grid* grd, parameters* param);



inline void GMRes(GENERIC_IMAGE_GMRES FunctionImage, FPfield *xkrylov, int xkrylovlen, FPfield *b, int m, int max_iter, double tol, EMfield* field, interpDensSpecies* ids, grid* grd, parameters* param) {
  
    // check size of m
    if (m > xkrylovlen) {
        std::cerr << "In GMRES the dimension of Krylov space(m) can't be > (length of krylov vector)/(# processors)" << std::endl;
        return;
    }
    
    bool GMRESVERBOSE = false;
    
    FPfield initial_error, normb, rho_tol, av, mu, htmp, tmp, delta = 0.001;
    
    // arrays for GMRes
    FPfield *r = new FPfield[xkrylovlen];
    FPfield *im = new FPfield[xkrylovlen];
    FPfield *v = new FPfield[xkrylovlen];
    FPfield *w = new FPfield[xkrylovlen];
    
    FPfield *s = new FPfield[m + 1];
    FPfield *cs = new FPfield[m + 1];
    FPfield *sn = new FPfield[m + 1];
    FPfield *y = new FPfield[m + 1];
    
    
    int k;
    eqValue(0.0, s, m + 1);
    eqValue(0.0, cs, m + 1);
    eqValue(0.0, sn, m + 1);
    eqValue(0.0, y, m + 1);
    
    
    // allocate H for storing the results from decomposition
    FPfield **H = newArr(FPfield, m + 1, m);
    for (int ii = 0; ii < m + 1; ii++)
        for (int jj = 0; jj < m; jj++)
            H[ii][jj] = 0;
   
    // allocate V
    FPfield **V = newArr(FPfield, xkrylovlen, m + 1);
    for (int ii = 0; ii < xkrylovlen; ii++)
        for (int jj = 0; jj < m + 1; jj++)
            V[ii][jj] = 0;
    
    
    
    if (GMRESVERBOSE) {
        std::cout << "------------------------------------" << std::endl;
        std::cout << "-             GMRES                -" << std::endl;
        std::cout << "------------------------------------" << std::endl;
        std::cout << std::endl;
    }
    
    for ( int itr = 0; itr < max_iter; itr++) {
        
        // r = b - A*x
        (*FunctionImage)(im, xkrylov, field, ids, grd, param);
        sub(r, b, im, xkrylovlen);
        
        // error as norm of th 1residual
        initial_error = normP(r, xkrylovlen);
        normb = normP(b, xkrylovlen);
        if (normb < 1E-16){
            //normb = 1.0;
            //initial_error = 0.0;
        }
        if (itr == 0) {
            
            std::cout << "Initial residual: " << initial_error << " norm b vector (source) = " << normb << std::endl;
            rho_tol = initial_error * tol;
            
            if ((initial_error / normb) <= tol) {
                
                std::cout << "GMRES converged without iterations: initial error < tolerance" << std::endl;
                delete[]r;
                delete[]im;
                delete[]s;
                delete[]v;
                delete[]cs;
                delete[]sn;
                delete[]w;
                delete[]y;
                delArr(H, m + 1);
                delArr(V, xkrylovlen);
                
                return;
            }
        }
        
        scale(v, r, (1.0 / initial_error), xkrylovlen);
        putColumn(V, v, 0, xkrylovlen);
        eqValue(0.0, s, m + 1);
        s[0] = initial_error;
        k = 0;
       
        
        // GMRes iteration
        while (rho_tol < initial_error && k < m) {
            
            // w= A*V(:,k)  );
            getColumn(v, V, k, xkrylovlen);
            (*FunctionImage) (w, v, field, ids, grd, param);
            putColumn(V, w, k + 1, xkrylovlen);
            av = normP(w, xkrylovlen);
            
            for ( int j = 0; j <= k; j++) {
                getColumn(v, V, j, xkrylovlen);
                H[j][k] = dotP(w, v, xkrylovlen);
                addscale(-H[j][k], w, v, xkrylovlen);
                
            }
            putColumn(V, w, k + 1, xkrylovlen);
            H[k + 1][k] = normP(w, xkrylovlen);
            
            if (av + delta * H[k + 1][k] == av) {
                
                for ( int j = 0; j <= k; j++) {
                    getColumn(v, V, j, xkrylovlen);
                    htmp = dotP(w, v, xkrylovlen);
                    H[j][k] = H[j][k] + htmp;
                    addscale(-htmp, w, v, xkrylovlen);
                }
                putColumn(V, w, k + 1, xkrylovlen);
                
                H[k + 1][k] = normP(w, xkrylovlen);
            }
            scale(w, (1.0 / H[k + 1][k]), xkrylovlen);
            
            putColumn(V, w, k + 1, xkrylovlen);
            
            if (0 < k) {
                
                for ( int j = 0; j < k; j++)
                    ApplyPlaneRotation(H[j + 1][k], H[j][k], cs[j], sn[j]);
                
                getColumn(y, H, k, m + 1);
            }
            
            mu = sqrt(H[k][k] * H[k][k] + H[k + 1][k] * H[k + 1][k]);
            cs[k] = H[k][k] / mu;
            sn[k] = -H[k + 1][k] / mu;
            H[k][k] = cs[k] * H[k][k] - sn[k] * H[k + 1][k];
            H[k + 1][k] = 0.0;
            
            ApplyPlaneRotation(s[k + 1], s[k], cs[k], sn[k]);
            initial_error = fabs(s[k]);
            k++;
        }
        
        k--;
        y[k] = s[k] / H[k][k];
        
        for ( int i = k - 1; i >= 0; i--) {
            tmp = 0.0;
            for ( int l = i + 1; l <= k; l++)
                tmp += H[i][l] * y[l];
            y[i] = (s[i] - tmp) / H[i][i];
            
        }
        
        
        for (int jj = 0; jj < xkrylovlen; jj++) {
            tmp = 0.0;
            for ( int l = 0; l < k; l++)
                tmp += y[l] * V[jj][l];
            xkrylov[jj] += tmp;
        }
        
        if (initial_error <= rho_tol) {
            
            std::cout << "GMRES converged at restart # " << itr << "; iteration #" << k << " with error: " << initial_error / rho_tol * tol << std::endl;
            
            // deallocate
            delete[]r;
            delete[]im;
            delete[]s;
            delete[]v;
            delete[]cs;
            delete[]sn;
            delete[]w;
            delete[]y;
            delArr(H, m + 1);
            delArr(V, xkrylovlen);
            
            return;
        }
        if (GMRESVERBOSE)
            std::cout << "Restart: " << itr << " error: " << initial_error / rho_tol * tol << std::endl;
        
    }
    
    std::cout << "GMRES not converged !! Final error: " << initial_error / rho_tol * tol << std::endl;
    
    delete[]r;
    delete[]im;
    delete[]s;
    delete[]v;
    delete[]cs;
    delete[]sn;
    delete[]w;
    delete[]y;
    delArr(H, m + 1);
    delArr(V, xkrylovlen);
    
    
    return;
}


inline void ApplyPlaneRotation(FPfield &dx, FPfield &dy, FPfield &cs, FPfield &sn) {
  FPfield temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}



#endif
