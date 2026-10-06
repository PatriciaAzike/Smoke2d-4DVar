 /*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Patricia Azike, 
Sandra Babyale
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

#ifndef TRANSPORT2D_H
#define TRANSPORT2D_H

#include <fclaw_include_all.h>

#include <fclaw_clawpatch_pillow.h>

/* Headers for both Clawpack 4.6 and  Clawpack 5.0 */
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>

/* Clawpack 4.6 headers */
#include <fc2d_clawpack46.h>  
#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46_fort.h>  
#include <clawpack46_user_fort.h>  
#include <fclaw2d_clawpatch46_fort.h>


#include "transport2d_fort.h"

#include "../smoke2d_options.h"
#include "../smoke2d_gauges.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

//struct fclaw_options;
//struct user_options;
struct fclaw_patch;
struct fclaw_domain;


void transport2d_setaux(struct fclaw_global *glob,
                        struct fclaw_patch *patch,
                        int blockno,
                        int patchno);

void transport2d_b4step2(struct fclaw_global *glob,
                         struct fclaw_patch *patch,
                         int blockno,
                         int patchno,
                         double t,
                         double dt);

/* --------------------------------- Sphere mappings ---------------------------------- */

#if 0
fclaw_map_context_t * fclaw_map_new_cubedsphere (const double scale[],
                                                     const double rotate[]);

fclaw_map_context_t * fclaw_map_new_pillowsphere (const double scale[],
                                                      const double rotate[]);
#endif                                                      


#ifdef __cplusplus
}
#endif

#endif
