/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef ADJOINT_USER_H
#define ADJOINT_USER_H

// #include "../all/advection_user.h"

#include <transport2d.h>
#include <transport2d_fort.h>


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

typedef struct adjoint_options
{
    int example;
    int pseudo_1d;
    int mdata;
    int initial_condition;
    double eps_1d;
    double eps_2d;
    double beta;
    double x0;
    double y0;

    double *dm; /* data observation*/
    const char *dm_string;

    double *W_eps; /* data weight*/
    const char *W_eps_string;

    double *xm;   /* x position of mdata */
    const char *xm_string;

    double *ym;   /* y position of m data */
    const char *ym_string;
    
    double *tm;   /* temporal position of m data */
    const char *tm_string;
    
    int is_registered;

} adjoint_options_t;


adjoint_options_t* adjoint_options_register (fclaw_app_t * app,
                                           const char *configfile);

void adjoint_options_store (fclaw_global_t* glob, adjoint_options_t* user);

const adjoint_options_t* adjoint_get_options(fclaw_global_t* glob);

void adjoint_link_solvers(fclaw_global_t *glob);

/* adjoint */
void adjoint_create_domain(fclaw_global_t *glob);
void adjoint_initialize(fclaw_global_t* glob);
void adjoint_finalize(fclaw_global_t* glob);


#ifdef __cplusplus
}
#endif

#endif
