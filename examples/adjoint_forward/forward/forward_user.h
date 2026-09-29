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

#ifndef FORWARD_USER_H
#define FORWARD_USER_H

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

typedef struct forward_options
{
    int example;
    int initial_condition;
    double W_f;
    double W_i;
    double beta;
    double x0;
    double y0;

    int is_registered;

} forward_options_t;


forward_options_t* forward_options_register (fclaw_app_t * app,
                                           const char *configfile);

void forward_options_store (fclaw_global_t* glob, forward_options_t* user);

const forward_options_t* forward_get_options(fclaw_global_t* glob);

void forward_link_solvers(fclaw_global_t *glob);

/* forward */
void forward_create_domain(fclaw_global_t *glob);
void forward_initialize(fclaw_global_t* glob);
void forward_finalize(fclaw_global_t* glob);


#ifdef __cplusplus
}
#endif

#endif
