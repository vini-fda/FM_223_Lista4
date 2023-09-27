/* SEGPLOT — This program produces a nonrigorous but nevertheless very
accurate picture of the stable manifold of the fixed point in the first
* quadrant of the Henon map with the usual parameters. It implements
* most of the basic ideas described in Sections 2 and 3 of the paper,
* but it does not require a knowledge of Lipschitz constants of the
* map or its first derivatives. This program is intended to be
* completely portable. it can be executed with the command line
* “segplot 6 0.2" to produce about 61,000 points on the stable
* manifold. See comments at the end of the program and additional
* details provided in the appendix.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define inbox(y) (y[0] >= boxmin[0] && y[0] <= boxmax[0] && \
                  y[1] >= boxmin[1] && y[1] <= boxmax[1])

#define SMALL 1.0e-06 /* minimum relative change in Y */
#define JUMP 1.0e-03 /* length of initial curve */
#define TINY 1.0e-16 /* minimum step size at any point */

double modval[2] = {0.0, 0.0}; /* nonzero if X or Y is
                                * taken modulo something */

/* Gamma(0) is the fixed point, here determined
* analytically ahead of time. */

// double endpoint[2] = {0.88389626793, 0.88389626793}; /* gamma(0) */

/* Gamma(1) is the fixed point together with a step of length JUMP in the
* direction of the unstable eigenvector of the inverse map, here
* determined analytically ahead of time. */
// double dgamma[2] = {JUMP, -JUMP*6.41246286050}; /* gamma(1)-gamma(0) */
/* BOXMIN/MAX is the lockout region. It should produce an
* accurate picture of the unstable manifold (of the fixed point for
the inverse map) in the smaller box [-3, -3] x [3, 3].
*/

double boxmin[2]={-3.4, -3.4}; /* min x, y coords of box */
double boxmax[2] = {3.4, 3.4}; /* max x, y coords of box */
double dbox[2] = {6.8, 6.8}; /* boxmax — boxmin */

/*---------------------------------------*/
/* MAP — return the Mth iterate of Y. This is the 2nd iterate of
* the inverse Henon map with the usual parameters, Eq. (16).
*/

void map(double *y, int m, double alpha, double beta) {
    int j;
    double xnext, ynext;

    for (j=0; j<2*m;j++) {
        xnext=y[1];
        ynext=(y[1]*y[1] + y[0] - alpha) / beta;
        y[0]=xnext;
        y[1]=ynext;
    }
    return;
}

// void map(double *y, int m, double alpha, double beta) {
//     int j;
//     double xnext, ynext;

//     for (j=0; j<2*m;j++) {
//         xnext=alpha-y[0]*y[0]+beta*y[1];
//         ynext=y[0];
//         y[0]=xnext;
//         y[1]=ynext;
//     }
//     return;
// }

/*---------------------------------------*/
/* YCHANGE — Determine how far T^m(s) is from Y in units of pixels.
* On return, YC holds T^m(s + ds).

* The return value is the max norm of Y-YC in pixels. A minimum value of
* SMALL is returned if Y-YC is too small.

* This handles modulo variables provided MAP returns them in the half-open
* interval [0, modval)
*/

double ychange (double s, double *y, double *yc, int m, double endpoint[2], double dgamma[2], double alpha, double beta) {
    double dx, dy;

    yc[0] = endpoint[0] + s * dgamma[0];
    yc[1] = endpoint[1] + s * dgamma[1];
    map(yc, m, alpha, beta);
    dx = fabs(y[0] - yc[0] - modval[0]) / dbox[0];
    dy = fabs(y[1] - yc[1] - modval[1]) / dbox[1];
    dx = max(dx, dy);

    return max(dx, SMALL);
}
/*---------------------------------------*/
/* REVISE — return an estimate of DS so that a small step from Y
* (= gamma(s)) produces an image that does not move more than about
* EPSILON from the image of Y.
*/

double revise(double s, double ds, double *y, double *yc, int m, double endpoint[2], double dgamma[2], double alpha, double beta) {
    double previous = ds, delta;

    while(((delta = ychange(s + ds, y, yc, m, endpoint, dgamma, alpha, beta)) > JUMP || 
          (ds *= 0.2 + 0.8 / delta) <= previous * 0.5) && ds > TINY)
        ds = previous = previous * 0.5;
    ds = min(ds, previous * 2.0);

    if (ds < TINY)
        fprintf(stderr, "stepsize reset from %lg to %lg at s=%lg, m=%d\n",
                ds, TINY, s, m);
    return max(ds, TINY);
}
/*---------------------------------------*/
/* ITERATE — iterate the map as follows.
* On entry, if M==0 or Y is not in the box, then Y is set to gamma(s).
* Otherwise, Y is assumed to contain the Mth image of some starting
* point and S is ignored. We iterate until K = =N or the (K + 1)st image
* of Y falls outside the box. In ether case, Y is set to the Kth iterate
* and we retum K.
*/
int iterate(double s, double *y, int m, int n, double endpoint[2], double dgamma[2], double alpha, double beta) {
    double prev[2];
    int ok; /*nonzero as long as we stay in the box */
    if (m == 0 || !inbox(y)) {
        y[0] = endpoint[0] + s * dgamma[0];
        y[1] = endpoint[1] + s * dgamma[1];
        m=0;
    }
    prev[0] = y[0]; prev[1] = y[1];
    if (m < n) do {
        map(y, 1, alpha, beta);
        if ((ok = inbox(y))) {
            prev[0] = y[0]; prev[1] = y[1];
        } else {
            y[0] = prev[0]; y[1] = prev[1];
        } 
    } while (ok && ++m < n);
    return m;
}
/*---------------------------------------*/
/* calc_manifold.
* n is the number of iterations to proceed. 
* smax is a number between 0 and 1 stating
*   the fraction of the inital curve that we want to plot. Compile the
* program with the Unix command line ‘make segplot”.
* Then the command “segplot 6 0.2" plots the 6th iterate of the initial
* line segment gamma under the inverse Henon map. Gamma is parametrized
* from 0 to 1, and we plot the 6th image of the fist 20% of the curve.
* The latter command produces about 61,000 dots.
*/

void calc_manifold(int n, double smax, double alpha, double beta, double endpoint[2], double eigvec[2], double* restrict xv, double* restrict yv, size_t capacity) {
    int ok = 0, m = 0;
    double s, y[2], yc[2];
    double ds = 1.0e-05; /* a reasonable first step size */
    size_t i = 0;
    double dgamma[2] = {JUMP*eigvec[0], JUMP*eigvec[1]};

    for (s = 0.0; s < smax; s = min(s + ds, smax)) {
        if(!ok || m < n)
            m = iterate(s, y, m, n, endpoint, dgamma, alpha, beta);
        ds = revise(s, ds, y, yc, m, endpoint, dgamma, alpha, beta);
        if (m < n) { /* try going one more time */
            map(y, 1, alpha, beta); map(yc, 1, alpha, beta);
            if (inbox(y) || inbox(yc))
                ds = revise(s, ds, y, yc, ++m, endpoint, dgamma, alpha, beta);
        }

        y[0] = yc[0]; y[1] = yc[1];

        if ((ok = inbox(y)) && m == n) {
            xv[i] = y[0];
            yv[i] = y[1];
            i++;
            if (i >= capacity) {
                fprintf(stderr, "error: capacity exceeded\n");
                return;
            }
        }
    }
}
