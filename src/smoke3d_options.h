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

#ifndef SMOKE3D_OPTIONS_H
#define SMOKE3D_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw_global;

typedef struct smoke3d_options
{

    int time_dependent_velocity;
    int forward_model;
    int adjoint_model;

    int is_registered;
}
smoke3d_options_t;



smoke3d_options_t*  smoke3d_options_register (fclaw_app_t * app,
                                              const char *section,
                                              const char *configfile);

void smoke3d_options_store (struct fclaw_global* glob, 
                            smoke3d_options_t* smoke3d_opt);

const smoke3d_options_t* smoke3d_get_options(struct fclaw_global* glob);



#ifdef __cplusplus
}
#endif

#endif /* SMOKE3D_OPTIONS_H */
