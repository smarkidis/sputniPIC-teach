#ifndef CG_H
#define CG_H

#include <math.h>

typedef void (*GENERIC_IMAGE) (FPfield *, FPfield  *, grid *);



inline bool CG(FPfield *xkrylov, int xkrylovlen, FPfield *b, int maxit, double tol, GENERIC_IMAGE FunctionImage, grid * grd) {
    
    // allocate residual, image, p, b, calculated on central points
    FPfield *r = new FPfield[xkrylovlen];
    FPfield *v = new FPfield[xkrylovlen];
    FPfield *z = new FPfield[xkrylovlen];
    FPfield *im = new FPfield[xkrylovlen];
    
    // aux variables
    FPfield c=0.0, t=0.0, d=0.0, initial_error=0.0, dotP=0.0;
    
    bool CONVERGED = false;
    bool CGVERBOSE = false;
    
    // i is the CG counter
    // don't use it for loops
    int i = 0;
    
    // initial guess for x: all the components are equal to 0
    for (int ii=0; ii < xkrylovlen; ii++)
        xkrylov[ii] = 0.0;
    
    // Compute r = b -Ax
    //sub(r, b, im, xkrylovlen);
    // v = r
    //eq(v, r, xkrylovlen);
    //c = dotP(r, r, xkrylovlen);
    
    (*FunctionImage) (im, xkrylov, grd);
    c = 0.0;
    for (int ii=0; ii < xkrylovlen; ii++){
        r[ii] = b[ii] - im[ii];
        v[ii] = r[ii];
        c +=  r[ii]*r[ii];
    }
    initial_error = sqrt(c);
    std::cout << "Initial error: " << initial_error << std::endl;
    
    if (CGVERBOSE) {
        std::cout << "------------------------------------" << std::endl;
        std::cout << "-               CG                 -" << std::endl;
        std::cout << "------------------------------------" << std::endl;
        std::cout << std::endl;
    }
    
    while (i < maxit) {
        
        // claculate image of Poisson equation
        (*FunctionImage) (z, v, grd);
        
        // t = c / dotP(v, z, xkrylovlen);
        dotP = 0.0;
        for (int ii=0; ii < xkrylovlen; ii++)
            dotP += v[ii]*z[ii];
        t = c / dotP;
        
        // x(i+1) = x + t*v - addscale(t, xkrylov, v, xkrylovlen);
        for (int ii=0; ii < xkrylovlen; ii++)
            xkrylov[ii] += t*v[ii];
        
        // r(i+1) = r - t*z - addscale(-t, r, z, xkrylovlen);
        for (int ii=0; ii < xkrylovlen; ii++)
            r[ii] -= t*z[ii];
        
        
        // d = dotP(r, r, xkrylovlen);
        d = 0.0;
        for (int ii=0; ii < xkrylovlen; ii++)
            d += r[ii]*r[ii];
        
        
        if (CGVERBOSE)
            std:: cout << "Iteration # " << i << " - norm of residual relative to initial error " << sqrt(d) / initial_error << std::endl;
        if (sqrt(d) < tol * initial_error) {
            std::cout << "CG converged at iteration # " << i << std::endl;
            CONVERGED = true;
            break;
        }
        else if (sqrt(d) > 10E8 * initial_error) {
            // HERE
            std::cerr << "CG not converging" << std::endl;
            std::cerr << "CG stopped" << std::endl;
           
            CONVERGED = false;
            break;
        }
        
        // addscale(1, d/c, v, r, xkrylovlen); DOUBLE CHECK THIS
        for ( int ii = 0; ii < xkrylovlen; ii++)
            v[ii] = v[ii] * d/c + r[ii];
        
        c = d;
        // update the CG iteration counter
        i++;
        
    }
    
    
    // deallocate
    delete[]r;
    delete[]im;
    delete[]v;
    delete[]z;
    
    // return true if converged
    return (CONVERGED);
}
#endif
