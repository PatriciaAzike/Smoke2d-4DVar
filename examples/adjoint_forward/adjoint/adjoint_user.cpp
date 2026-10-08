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
#include "adjoint_fort.h"
#include "user_run_fort.h"

#include <fclaw_clawpatch.h>

void adjoint_problem_setup(fclaw_global_t *glob)
{
    const adjoint_options_t* user = adjoint_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24d   %s",user->mdata,"\% mdata\n");
        fprintf(f,  "%-24d   %s",user->initial_condition,"\% initial_condition\n");
        fprintf(f,  "%-24.6f   %s",user->eps_1d,"\% epsilon-1d\n");
        fprintf(f,  "%-24.6f   %s",user->eps_2d,"\% epsilon-2d\n");
        fprintf(f,  "%-24.6f   %s",user->beta,"\% beta\n");
        fprintf(f,  "%-24.6f   %s",user->x0,"\% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,"\% y0\n");
        fprintf(f,  "%-24.6f   %s",fclaw_opt->tfinal,"\% tfinal\n");


        for (int i = 0; i < user->mdata; i++)
        {
            fprintf(f,  "%-24.6f   %s",user->xm[i],"\% xm\n");
            fprintf(f,  "%-24.6f   %s",user->ym[i],"\% ym\n");
            fprintf(f,  "%-24.6f   %s",user->tm[i],"\% tm\n");
            fprintf(f,  "%-24.6f   %s",user->dm[i],"\% dm\n");
            fprintf(f,  "%-24.6f   %s",user->W_eps[i],"\% W_eps\n");
        }

        fprintf(f,  "%-24s   %s",user->pseudo_1d_experiment ? "T" : "F",
                "\% pseudo-1d-experiment\n");
        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    ADJOINT_SETPROB();
}

void adjoint_link_solvers(fclaw_global_t *glob)
{
    //const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    fclaw_vtable_t *claw_vt = fclaw_vt(glob);
    claw_vt->problem_setup = adjoint_problem_setup;

    fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);

    clawpack46_vt->fort_setprob   = &ADJOINT_SETPROB;
    clawpack46_vt->fort_qinit     = &ADJOINT_QINIT;
    
    clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
    clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;
    clawpack46_vt->fort_src2      = &ADJOINT_SRC2;

    clawpack46_vt->fort_setaux    = &ADJOINT_SETAUX;

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
