/*
  Copyright (c) 2012-2025 Carsten Burstedde, Donna Calhoun, Patricia Azike
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

#ifndef USER_RUN_FORT_H
#define USER_RUN_FORT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <fclaw2d_clawpatch46_fort.h>
  
#define CLAWPATCH46_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(clawpatch46_fort_tag4refinement, \
                                                                CLAWPATCH46_FORT_TAG4REFINEMENT)
void CLAWPATCH46_FORT_TAG4REFINEMENT (const int* mx,const int* my,
                        const int* mbc,const int* meqn,
                        const double* xlower, const double* ylower,
                        const double* dx, const double* dy,
                        const int* blockno,
                        double q[],
                        const double* refine_threshold,
                        const int* init_flag,
                        int* tag_patch); 

#define CLAWPATCH46_FORT_TAG4COARSENING FCLAW_F77_FUNC(clawpatch46_fort_tag4coarsening, \
                                                        CLAWPATCH46_FORT_TAG4COARSENING)
void CLAWPATCH46_FORT_TAG4COARSENING (const int* mx, const int* my,
                        const int* mbc, const int* meqn,
                        double xlower[], 
                                                double ylower[],
                        const double* dx, const double* dy,
                        const int* blockno,
                        double q0[],double q1[],
                        double q2[],double q3[],
                        const double* coarsen_threshold,
                        const int* init_flag,
                        int* tag_patch);
    
#ifdef __cplusplus
}
#endif


#endif
