/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton, Patricia Azike
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

#include "forward_user.h"
#include "forward_fort.h"
#include "model_user.h"
#include "user_run_fort.h"

#include "fc2d_clawpack46_options.h"
#include <fclaw_clawpatch.h>

void forward_problem_setup(fclaw_global_t *glob)
{
    const forward_options_t* user = forward_get_options(glob);
    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.6f   %s",user->W_f,"\% W_f\n");
        fprintf(f,  "%-24.6f   %s",user->W_i,"\% W_i\n");
        fprintf(f,  "%-24.6f   %s",user->beta,"\% beta\n");
        fprintf(f,  "%-24.6f   %s",user->x0,"\% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,"\% y0\n");
        fprintf(f,"%-24d %s\n",fclaw_opt->moving_gauges,"\% moving_gauges");

        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    FORWARD_SETPROB();
}

void forward_link_solvers(fclaw_global_t *glob)
{
    //const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    fclaw_vtable_t *claw_vt = fclaw_vt(glob);
    claw_vt->problem_setup = forward_problem_setup;

    fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);

    clawpack46_vt->fort_setprob   = &FORWARD_SETPROB;
    clawpack46_vt->fort_qinit     = &FORWARD_QINIT;
    clawpack46_vt->fort_rpn2      = &RPN2_FORWARD;
    clawpack46_vt->fort_rpt2      = &RPT2_FORWARD;
    //clawpack46_vt->fort_rpn2_cons = &RPN2_CONSERVATIVE_UPDATE; 

    clawpack46_vt->fort_src2      = &FORWARD_SRC2;
    clawpack46_vt->fort_setaux    = &FORWARD_SETAUX;

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d2->fort_copy_face          = FCLAW2D_CLAWPATCH46_FORT_COPY_FACE;
    clawpatch_vt->d2->fort_copy_corner        = FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER;
    clawpatch_vt->d2->fort_average_face       = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE;
    clawpatch_vt->d2->fort_average_corner     = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER;
    clawpatch_vt->d2->fort_average2coarse     = FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE;
    clawpatch_vt->d2->fort_interpolate_face   = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE;
    clawpatch_vt->d2->fort_interpolate_corner = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER;
    clawpatch_vt->d2->fort_interpolate2fine   = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE;
    clawpatch_vt->d2->fort_local_ghost_pack   = FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;
    clawpatch_vt->d2->fort_timeinterp         = FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP;
    clawpatch_vt->d2->fort_tag4refinement     = CLAWPATCH46_FORT_TAG4REFINEMENT;
    clawpatch_vt->d2->fort_tag4coarsening     = CLAWPATCH46_FORT_TAG4COARSENING;
}

void model_problem_setup(fclaw_global_t *glob)
{
    /* Write setprob.data from [model-user] in model_options.ini.
       Do not copy ../forward/setprob.data: the model run starts
       before the forward run has written that file. */
    const model_options_t* user = model_get_options(glob);
    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24d   %s",user->example,"\% example\n");
        fprintf(f,"%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,"%-24.6f   %s",user->W_f,"\% W_f\n");
        fprintf(f,"%-24.6f   %s",user->W_i,"\% W_i\n");
        fprintf(f,"%-24.6f   %s",user->beta,"\% beta\n");
        fprintf(f,"%-24.6f   %s",user->x0,"\% x0\n");
        fprintf(f,"%-24.6f   %s",user->y0,"\% y0\n");
        fprintf(f,"%-24d %s\n",fclaw_opt->moving_gauges,"\% moving_gauges");
        fclose(f);
    }
    fclaw_domain_barrier(glob->domain);
    FORWARD_SETPROB();
}

void model_link_solvers(fclaw_global_t *glob)
{
    //const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    fclaw_vtable_t *claw_vt = fclaw_vt(glob);
    claw_vt->problem_setup = model_problem_setup;

    fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);

    clawpack46_vt->fort_setprob   = &FORWARD_SETPROB;
    clawpack46_vt->fort_qinit     = &FORWARD_QINIT;
    clawpack46_vt->fort_rpn2      = &RPN2_FORWARD;
    clawpack46_vt->fort_rpt2      = &RPT2_FORWARD;
    //clawpack46_vt->fort_rpn2_cons = &RPN2_CONSERVATIVE_UPDATE; 

    clawpack46_vt->fort_src2      = &FORWARD_SRC2;
    clawpack46_vt->fort_setaux    = &FORWARD_SETAUX;

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d2->fort_copy_face          = FCLAW2D_CLAWPATCH46_FORT_COPY_FACE;
    clawpatch_vt->d2->fort_copy_corner        = FCLAW2D_CLAWPATCH46_FORT_COPY_CORNER;
    clawpatch_vt->d2->fort_average_face       = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_FACE;
    clawpatch_vt->d2->fort_average_corner     = FCLAW2D_CLAWPATCH46_FORT_AVERAGE_CORNER;
    clawpatch_vt->d2->fort_average2coarse     = FCLAW2D_CLAWPATCH46_FORT_AVERAGE2COARSE;
    clawpatch_vt->d2->fort_interpolate_face   = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_FACE;
    clawpatch_vt->d2->fort_interpolate_corner = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE_CORNER;
    clawpatch_vt->d2->fort_interpolate2fine   = FCLAW2D_CLAWPATCH46_FORT_INTERPOLATE2FINE;
    clawpatch_vt->d2->fort_local_ghost_pack   = FCLAW2D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK;
    clawpatch_vt->d2->fort_timeinterp         = FCLAW2D_CLAWPATCH46_FORT_TIMEINTERP;
    clawpatch_vt->d2->fort_tag4refinement     = CLAWPATCH46_FORT_TAG4REFINEMENT;
    clawpatch_vt->d2->fort_tag4coarsening     = CLAWPATCH46_FORT_TAG4COARSENING;
}
