/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Patricia Azike
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

#include "../smoke3d_options.h"

#include "transport2d.h"


void transport2d_setaux(fclaw_global_t *glob,
                        fclaw_patch_t *patch,
                        int blockno,
                        int patchno)
{
    int maux;
    double *aux;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);


    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *area, *edgelengths,*curvature;
    fclaw_clawpatch_2d_metric_scalar(glob, patch,&area,&edgelengths,
                                    &curvature);

    /* Handles both non-conservative (ex 1-2) and conservative (ex 3-4) forms */
    FCLAW_ASSERT(maux == 7);
    TRANSPORT2D_SETAUX_METRIC(&blockno, &mx,&my,&mbc, &xlower,&ylower,
                                &dx,&dy, area, edgelengths, aux, &maux);

    /* If velocity field is time dependent, it will also be set in 
       b4step2 */
    double *xp, *yp, *zp, *xd, *yd, *zd, *areas;
    fclaw_clawpatch_2d_metric_data(glob,patch,&xp,&yp,&zp,&xd,&yd,&zd,&areas);

    double *xnormals,*ynormals,*xtangents,*ytangents,*surfnormals;
    fclaw_clawpatch_2d_metric_vector(glob,patch,
                                    &xnormals, &ynormals,
                                    &xtangents, &ytangents,
                                    &surfnormals);

    /* If the velocity field is time dependent, b4step2 will set the velocity.  
       No need to do it here. */
    const smoke3d_options_t *smoke_opt = smoke3d_get_options(glob);
    double t = 0;
    if (smoke_opt->time_dependent_velocity == 0)
        TRANSPORT2D_SETAUX_VELOCITY(&blockno, &mx, &my, &mbc,
                                    &dx, &dy, &xlower, &ylower,
                                    &t, xp, yp, zp, 
                                    xnormals,ynormals, surfnormals,
                                    aux,&maux);

}
