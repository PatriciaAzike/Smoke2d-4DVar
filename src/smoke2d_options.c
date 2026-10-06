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

#include "smoke2d_options.h"

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_pointer_map.h>


static void *
smoke2d_register (smoke2d_options_t *smoke_opt, sc_options_t * opt)
{
    sc_options_add_bool (opt, 0, "time-dependent-velocity", 
                         &smoke_opt->time_dependent_velocity, 0,
                         "Time dependent velocity (T/F) [F]");


   sc_options_add_bool (opt, 0, "forward-model", 
                         &smoke_opt->forward_model, 0,
                         "Forward model (T/F) [T]");


   sc_options_add_bool (opt, 0, "adjoint-model", 
                         &smoke_opt->adjoint_model, 0,
                         "Adjoint model (T/F) [F]");

    smoke_opt->is_registered = 1;
    return NULL;
}

static
fclaw_exit_type_t
smoke2d_postprocess (smoke2d_options_t *smoke_opt)
{
    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
smoke2d_check (smoke2d_options_t *smoke_opt)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}


static void
smoke2d_destroy (smoke2d_options_t *smoke_opt)
{
}


/* ------- Generic option handling routines that call above routines ----- */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    smoke2d_options_t *smoke2d_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    smoke2d_opt = (smoke2d_options_t*) package;

    return smoke2d_register(smoke2d_opt,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    smoke2d_options_t *smoke_opt = (smoke2d_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(smoke_opt->is_registered);

    /* Convert strings to arrays */
    return smoke2d_postprocess (smoke_opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    smoke2d_options_t           *smoke_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    smoke_opt = (smoke2d_options_t*) package;

    return smoke2d_check(smoke_opt);
}


static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    smoke2d_options_t *smoke_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    smoke_opt = (smoke2d_options_t*) package;
    FCLAW_ASSERT (smoke_opt->is_registered);

    smoke2d_destroy (smoke_opt);

    FCLAW_FREE (smoke_opt);
}


static const
fclaw_app_options_vtable_t smoke2d_options_vtable =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};


/* ------------- smoke_opt options access functions --------------------- */

#if 0
fc2d_clawpack46_options_t*  fc2d_clawpack46_options_register (fclaw_app_t * app,
                                                              const char *section,
                                                              const char *configfile)
{
    fc2d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);

    clawopt = FCLAW_ALLOC (fc2d_clawpack46_options_t, 1);
    fclaw_app_options_register (app, section, configfile,
                                &clawpack46_options_vtable, clawopt);
    
    fclaw_app_set_attribute(app, section, clawopt);
    return clawopt;
}
#endif


smoke2d_options_t*  smoke2d_options_register (fclaw_app_t * app,
                                              const char *section,
                                              const char *configfile)
{
    FCLAW_ASSERT (app != NULL);

    smoke2d_options_t *smoke_opt = FCLAW_ALLOC (smoke2d_options_t, 1);
    fclaw_app_options_register (app, section, configfile, &smoke2d_options_vtable,
                                smoke_opt);

    fclaw_app_set_attribute(app,section,smoke_opt);
    return smoke_opt;
}


const smoke2d_options_t* smoke2d_get_options(fclaw_global_t *glob)
{
    smoke2d_options_t* smoke_opt = (smoke2d_options_t*) 
    fclaw_pointer_map_get(glob->options, "smoke2d");
    FCLAW_ASSERT(smoke_opt != NULL);
    return smoke_opt;
}

void smoke2d_options_store (fclaw_global_t* glob, smoke2d_options_t* smoke_opt)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"smoke2d") == NULL);
    fclaw_pointer_map_insert(glob->options, "smoke2d", smoke_opt, NULL);
}




