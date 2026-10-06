/*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Patricia Azike
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TRANSPORT2D_FORT_H
#define TRANSPORT2D_FORT_H

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

/* ---------------------- SET METRIC TERMS ---------------------------------------- */

#define TRANSPORT2D_SETAUX_METRIC FCLAW_F77_FUNC(transport2d_setaux_metric, \
                                                   TRANSPORT2D_SETAUX_METRIC)

void TRANSPORT2D_SETAUX_METRIC(const int* blockno, const int* mx, const int* my,
                                 const int* mbc, 
                                 const double* xlower, const double* ylower,
                                 const double* dx, const double* dy, 
                                 double area[],double edgelengths[],
                                 double aux[],const int* maux);


/* ------------------------------- SET_VELOCITY --------------------------------------- */

#define TRANSPORT2D_SETAUX_VELOCITY FCLAW_F77_FUNC(transport2d_setaux_velocity, \
                                             TRANSPORT2D_SETAUX_VELOCITY)

void TRANSPORT2D_SETAUX_VELOCITY(const int* blockno, const int* mx, const int* my,
                              const int* mbc, const double* dx, const double* dy,
                              const double* xlower, const double* ylower,
                              const double *t, double xp[], double yp[], 
                              double zp[], double xnormals[],double ynormals[],
                              double surfnormals[], double aux[],const int* maux);


/* ------------------------------- Riemann solvers ------------------------------------ */

#define RPN2_CONSERVATIVE_UPDATE FCLAW_F77_FUNC(rpn2_conservative_update, \
                                                RPN2_CONSERVATIVE_UPDATE)

void RPN2_CONSERVATIVE_UPDATE(const int* meqn, const int* maux, 
                              const int* idir, const int* iface,
                              double q[], double aux_center[], 

                              double aux_edge[], double flux[]);

#define RPN2_FORWARD FCLAW_F77_FUNC(rpn2_forward, \
                          RPN2_FORWARD)
void RPN2_FORWARD(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves, 
                          const int* mbc, const int* mx, 
                          double ql[], double qr[],
                          double auxl[], double auxr[], 
                          double fwave[], double s[], 
                          double amdq[], double apdq[]);


#define RPT2_FORWARD FCLAW_F77_FUNC(rpt2_forward, \
                       RPT2_FORWARD)

void RPT2_FORWARD(const int* ixy, const int* maxm, 
                       const int* meqn, const int* mwaves,
                       const int* mbc, const int* mx, 
                       double ql[], double qr[],
                       double aux1[], double aux2[], 
                       double aux3[], const int* imp,
                       double dsdq[], double bmasdq[], double bpasdq[]);


#define RPN2_ADJOINT FCLAW_F77_FUNC(rpn2_adjoint, \
                          RPN2_ADJOINT)
void RPN2_ADJOINT(const int* ixy, const int* maxm, const int* meqn, 
                          const int* mwaves, 
                          const int* mbc, const int* mx, 
                          double ql[], double qr[],
                          double auxl[], double auxr[], 
                          double wave[],double s[], 
                          double amdq[], double apdq[]);


#define RPT2_ADJOINT FCLAW_F77_FUNC(rpt2_adjoint, \
                       RPT2_ADJOINT)

void RPT2_ADJOINT(const int* ixy, const int* maxm, 
                       const int* meqn, const int* mwaves,
                       const int* mbc, const int* mx, 
                       double ql[], double qr[],
                       double aux1[], double aux2[], 
                       double aux3[], const int* imp,
                       double asdq[], double bmasdq[], double bpasdq[]);


/* -------------------------- Non-mapped Riemann solvers ----------------------------- */

/* These are used in the filament example */

#define CLAWPACK46_RPN2ADV FCLAW_F77_FUNC(clawpack46_rpn2adv,CLAWPACK46_RPN2ADV)
void CLAWPACK46_RPN2ADV(const int* ixy,const int* maxm, const int* meqn, 
                        const int* mwaves,
                        const int* mbc,const int* mx, 
                        double ql[], double qr[],
                        double auxl[], double auxr[], 
                        double wave[],
                        double s[], double amdq[], double apdq[]);

#define CLAWPACK46_RPT2ADV FCLAW_F77_FUNC(clawpack46_rpt2adv, CLAWPACK46_RPT2ADV)
void CLAWPACK46_RPT2ADV(const int* ixy, const int* maxm, const int* meqn, 
                        const int* mwaves,
                        const int* mbc, const int* mx, 
                        double ql[], double qr[],
                        double aux1[], double aux2[], 
                        double aux3[], const int* imp,
                        double dsdq[], double bmasdq[], double bpasdq[]);

/* ------------------------------------ Source terms  -------------------------------------- */

#define SRC2_FORWARD    FCLAW_F77_FUNC(src2_forward,   SRC2_FORWARD)
void SRC2_FORWARD(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx,const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, double q[],
                     const int* maux, double aux[], const double* t,
                     const double* dt);


#define SRC2_ADJOINT    FCLAW_F77_FUNC(src2_adjoint,   SRC2_ADJOINT)
void SRC2_ADJOINT(const int* maxmx, const int* maxmy, const int* meqn,
                     const int* mbc, const int* mx,const int* my,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy, double q[],
                     const int* maux, double aux[], const double* t,
                     const double* dt);

/* ------------------------------------ Tagging  -------------------------------------- */

#define USER_EXCEEDS_THRESHOLD FCLAW_F77_FUNC(user_exceeds_threshold, \
                                              USER_EXCEEDS_THRESHOLD)

int USER_EXCEEDS_THRESHOLD(int* blockno,
                           double qval[], 
                           double* qmin, double *qmax,
                           double quad[], 
                           double *dx, double *dy, 
                           double *xc, double *yc, 
                           int* tag_threshold, 
                           int* init_flag,
                           int* is_ghost);


/* ------------------------------------ Tagging  -------------------------------------- */
#define TRANSPORT_MAPC2M_LATLONG FCLAW_F77_FUNC(transport_mapc2m_latlong, \
                                                 TRANSPORT_MAPC2M_LATLONG)

void TRANSPORT_MAPC2M_LATLONG(double* xc, double *yc, 
                              double* xp,double *yp, double*zp);

/* ------------------------------------ GAUGES -------------------------------------- */


#define MAP_LATLONG2UNITREGION FCLAW_F77_FUNC(map_latlong2unitregion, \
                                              MAP_LATLONG2UNITREGION)

void MAP_LATLONG2UNITREGION(double *x_long, double *y_lat, double *xc, double* yc);


#define MAP_BRICK2LATLONG FCLAW_F77_FUNC(map_brick2latlong,MAP_BRICK2LATLONG)
void MAP_BRICK2LATLONG(const int* blockno,const double* xc, double* yc,
                       double *x_long,double *y_lat);




#define SMOKE2D_UPDATE_GAUGE FCLAW_F77_FUNC(smoke2d_update_gauge, \
                                            SMOKE2D_UPDATE_GAUGE)
void SMOKE2D_UPDATE_GAUGE (int* blockno, int* mx,int* my, 
                           int* mbc, int* meqn, 
                           double* xlower,double* ylower, 
                           double* dx, double* dy, 
                           double *q,
                           double* xc,double* yc,
                           double *qvar);

#ifdef __cplusplus
}
#endif

#endif
