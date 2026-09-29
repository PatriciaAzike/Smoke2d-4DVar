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

#include "adjoint_user.h"

static void *
adjoint_register (adjoint_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0 = cart (brick); 1 = 5-patch; 2 = 5-patch [0]");

    sc_options_add_int (opt, 0, "mdata", &user->mdata, 1,
                        "number of data points [1]");

    sc_options_add_int (opt, 0, "initial-condition", &user->initial_condition, 0,
                           "Initial conditions (zero by default) [0]");

    sc_options_add_double (opt, 0, "epsilon-1d", &user->eps_1d, 1e-3,
                            "width of heat kernel for pseudo-1d case [1e-3]");

    sc_options_add_double (opt, 0, "epsilon-2d", &user->eps_2d, 1e-2,
                            "width of heat kernel for full-2d case [1e-2]");

    sc_options_add_double (opt, 0, "beta", &user->beta, 10,
                            "width of Gaussian initial cond. [10]");

    sc_options_add_double (opt, 0, "x0", &user->x0, 0.3,
                            "x-coordinate of the center [0.3]");

    sc_options_add_double (opt, 0, "y0", &user->y0, 0.58,
                            "y-coordinate of the center [0.58]");

    fclaw_options_add_double_array(opt, 0, "xm", &user->xm_string,
                                   NULL, &user->xm, user->mdata,
                                   "[user] x position of m data");

    fclaw_options_add_double_array(opt, 0, "ym", &user->ym_string,
                                   NULL, &user->ym, user->mdata,
                                   "[user] y positin of m data");

    fclaw_options_add_double_array(opt, 0, "tm", &user->tm_string,
                                   NULL, &user->tm, user->mdata,
                                   "[user] temporal position of m data");

    fclaw_options_add_double_array(opt, 0, "dm", &user->dm_string,
                                   NULL, &user->dm, user->mdata,
                                   "[user] data observation");

    fclaw_options_add_double_array(opt, 0, "W_eps", &user->W_eps_string,
                                   NULL, &user->W_eps, user->mdata,
                                   "[user] weight assoc. with error at m data.");

    sc_options_add_int (opt, 0, "pseudo-1d", &user->pseudo_1d, 1,
                        "dimension of code: 1 = psuedo-1d; 2 = full-2d; [1]");


    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
adjoint_postprocess(adjoint_options_t *user)
{
    fclaw_options_convert_double_array(user->xm_string, &user->xm, user->mdata);
    fclaw_options_convert_double_array(user->ym_string, &user->ym, user->mdata);
    fclaw_options_convert_double_array(user->tm_string, &user->tm, user->mdata);
    fclaw_options_convert_double_array(user->dm_string, &user->dm, user->mdata);
    fclaw_options_convert_double_array(user->W_eps_string, &user->W_eps, user->mdata);
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
adjoint_check (adjoint_options_t *user)
{
    return FCLAW_NOEXIT;
}


static void
adjoint_destroy (adjoint_options_t *user)
{
    fclaw_options_destroy_array (user->xm);
    fclaw_options_destroy_array (user->ym);
    fclaw_options_destroy_array (user->tm);
    fclaw_options_destroy_array (user->dm);
    fclaw_options_destroy_array (user->W_eps);
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    adjoint_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (adjoint_options_t*) package;

    return adjoint_register(user,opt);
}


static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    adjoint_options_t *user = (adjoint_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return adjoint_postprocess (user);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    adjoint_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (adjoint_options_t*) package;
    return adjoint_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    adjoint_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (adjoint_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    adjoint_destroy (user);

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

adjoint_options_t* adjoint_options_register (fclaw_app_t * app,
                                           const char *configfile)
{
    adjoint_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (adjoint_options_t, 1);
    fclaw_app_options_register (app,"adjoint-user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"adjoint-user",user);
    return user;
}

void adjoint_options_store (fclaw_global_t* glob, adjoint_options_t* user)
{
    fclaw_global_options_store(glob, "adjoint-user", user);
}

const adjoint_options_t* adjoint_get_options(fclaw_global_t* glob)
{
    return (adjoint_options_t*) fclaw_global_get_options(glob, "adjoint-user");
}
