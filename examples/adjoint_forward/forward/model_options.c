/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "model_user.h"

static void *
model_register (model_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0 = cart (brick); 1 = 5-patch; 2 = 5-patch [0]");

    sc_options_add_int (opt, 0, "initial-condition", &user->initial_condition, 0,
                           "Initial conditions (zero by default) [0]");

    sc_options_add_double (opt, 0, "W_f", &user->W_f, 0.25,
                        "weight assoc. with error in forcing [0.25]");


    sc_options_add_double (opt, 0, "W_i", &user->W_i, 5e-2,
                            "weight assoc. with error in initial cond. [5e-2]");

    sc_options_add_double (opt, 0, "beta", &user->beta, 10,
                            "width of Gaussian initial cond. [10]");

    sc_options_add_double (opt, 0, "x0", &user->x0, 0.3,
                            "x-coordinate of the center [0.3]");

    sc_options_add_double (opt, 0, "y0", &user->y0, 0.58,
                            "y-coordinate of the center [0.58]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
model_postprocess(model_options_t *user)
{
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
model_check (model_options_t *user)
{
    return FCLAW_NOEXIT;
}


static void
model_destroy (model_options_t *user)
{
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    model_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (model_options_t*) package;

    return model_register(user,opt);
}


static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    model_options_t *user = (model_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return model_postprocess (user);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    model_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (model_options_t*) package;
    return model_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    model_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (model_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    model_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};


/* ------------- User options access functions --------------------- */

model_options_t* model_options_register (fclaw_app_t * app,
                                           const char *configfile)
{
    model_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (model_options_t, 1);
    fclaw_app_options_register (app,"model-user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"model-user",user);
    return user;
}

void model_options_store (fclaw_global_t* glob, model_options_t* user)
{
    fclaw_global_options_store(glob, "model-user", user);
}

const model_options_t* model_get_options(fclaw_global_t* glob)
{
    return (model_options_t*) fclaw_global_get_options(glob, "model-user");
}
