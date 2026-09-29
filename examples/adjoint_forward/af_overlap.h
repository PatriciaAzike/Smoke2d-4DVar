/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton,
Hannes Brandt, Patricia Azike
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

#ifndef AF_OVERLAP_H
#define AF_OVERLAP_H

#include <fclaw_options.h>

#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw_global.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting below */    
#endif

typedef struct overlap_adjoint_metadata
{
  fclaw_global_t   *glob;
  fclaw_domain_t   *adjoint_domain;

  int meqn;
  int isset;
} overlap_adjoint_metadata_t;

typedef struct overlap_forward_metadata
{
  fclaw_global_t   *glob;
  fclaw_domain_t   *domain;

  sc_array_t         *query_points;
  size_t              cell_idx;
  /* Number of cell-centered overlap samples stored for each target patch.
   * With AMR this is mx*my, not one point per patch. */
  int                 num_cells_in_patch;
} overlap_forward_metadata_t;



/* One overlap point is one target cell center.  The point stores its physical
 * coordinates for the ForestClaw overlap search, and also the target cell
 * location so the interpolated value can be written back into the right q/aux
 * entry after the source patch is found. */
typedef struct overlap_point 
{
  size_t  lnum;
  double  xy[2];  // In [ax,bx]x[ay,by] forward domain coordinates

  /* Target cell indices and patch dimensions.  ForestClaw gives us the patch
   * pointer during query construction, but only this overlap_point_t survives
   * through fclaw_overlap_exchange. */
  int i;
  int j;
  int mx;
  int my;
  int mbc;

  int maux;
  double* aux_forward;  /* Store adjoint solution here */

  int meqn;
  double* q_forward;  /* Store representer/model solution updates here */

  /* This is needed so we can interpolate from the adjoint patch. */
  overlap_adjoint_metadata_t   adjoint_metadata;
}
overlap_point_t;

// Geometry for the adjoint problem
typedef struct adjoint_geometry
{
    fclaw_options_t *fclaw_opt;
    fclaw_block_t *blocks;
}
adjoint_geometry_t;

typedef struct overlap_patch_user
{
    int init_flag;
    void *meta;
    adjoint_geometry_t *geo;
    int rep_idx;
    double *beta;
} overlap_patch_user_t;


void create_query_points (overlap_forward_metadata_t * c);

/*int apply_inverse_adjoint_mapping (overlap_point_t * op, double xy[2],
                                    int blockno, overlap_geometry_t * geo);*/

int overlap_interpolate (fclaw_domain_t * domain, 
                         fclaw_patch_t * patch,
                         int blockno, int patchno, 
                         void *point, void *user);

int overlap_interpolate_patch (fclaw_domain_t * domain, 
                               fclaw_patch_t * patch,
                               int blockno, int patchno, 
                               void *point, void *user);

int overlap_interpolate_rep_to_initial_estimate (fclaw_domain_t *representer_domain, 
                               fclaw_patch_t  *representer_patch,
                               int blockno, int patchno, 
                               void *point, void *user);
                         
void output_query_points (overlap_forward_metadata_t * c);


#if 0
void add_cell_centers (fclaw_domain_t * domain, fclaw_patch_t * patch,
                       int blockno, int patchno, void *user);
#endif


/* ------------------------------ FORTRAN SOLUTION ------------------------------ */
#define SETAUX_FROM_ADJOINT FCLAW_F77_FUNC(setaux_from_adjoint, \
                                            SETAUX_FROM_ADJOINT)
void SETAUX_FROM_ADJOINT(const int* mx,
                         const int* my,
                         const int* mbc,
                         const int* meqn,
                         const int* maux,
                         double adjoint[],
                         double aux_representer[]);


#define SETAUX_FROM_ADJOINT_POINT FCLAW_F77_FUNC(setaux_from_adjoint_point, \
                                                 SETAUX_FROM_ADJOINT_POINT)
void SETAUX_FROM_ADJOINT_POINT(const int* i,
                               const int* j,
                               const int* mx,
                               const int* my,
                               const int* mbc,
                               const int* meqn,
                               const int* maux,
                               double qvar[],
                               double aux_representer[]);



#define SETQ_FROM_ADJOINT FCLAW_F77_FUNC(setq_from_adjoint, \
                                            SETQ_FROM_ADJOINT)
void SETQ_FROM_ADJOINT(const int* mx,
                         const int* my,
                         const int* mbc,
                         const int* meqn,
                         double adjoint[],
                         double q_representer[]);


#define SETQ_FROM_ADJOINT_POINT FCLAW_F77_FUNC(setq_from_adjoint_point, \
                                               SETQ_FROM_ADJOINT_POINT)
void SETQ_FROM_ADJOINT_POINT(const int* i,
                             const int* j,
                             const int* mx,
                             const int* my,
                             const int* mbc,
                             const int* meqn,
                             double qvar[],
                             double q_representer[]);


#define UPDATEQ_FROM_REPRESENTER   FCLAW_F77_FUNC(updateq_from_representer, \
                                                  UPDATEQ_FROM_REPRESENTER)
void UPDATEQ_FROM_REPRESENTER(const int* mx, 
                              const int* my,
                              const int* mbc, 
                              const int* meqn,
                              double q_rep[],
                              double q_model[],
                              double *beta_j);


#define UPDATEQ_FROM_REPRESENTER_POINT FCLAW_F77_FUNC(updateq_from_representer_point, \
                                                      UPDATEQ_FROM_REPRESENTER_POINT)
void UPDATEQ_FROM_REPRESENTER_POINT(const int* i,
                                    const int* j,
                                    const int* mx,
                                    const int* my,
                                    const int* mbc,
                                    const int* meqn,
                                    double qvar[],
                                    double q_model[],
                                    double *beta_j);

void dgesv_(int *n,
                int *nrhs,
                double *A,
                int *lda,
                int *ipiv,
                double *B,
                int *ldb,
                int *info);

#ifdef __cplusplus
}
#endif

#endif
