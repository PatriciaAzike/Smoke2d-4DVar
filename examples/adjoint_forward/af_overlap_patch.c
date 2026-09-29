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

#include "af_overlap.h"


#include <fclaw_forestclaw.h>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>


static
void interpolate_point(fclaw_global_t *glob,
                       fclaw_patch_t *patch,
                       double xref,
                       double yref,
                       double qvar[])
{
    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob, patch, &q, &meqn);

    double *aux;
    int maux;
    fclaw_clawpatch_aux_data(glob, patch, &aux, &maux);

    int mx, my, mbc;
    double xlower, ylower, dx, dy;
    fclaw_clawpatch_2d_grid_data(glob, patch, &mx, &my, &mbc,
                                 &xlower, &ylower, &dx, &dy);

    xlower = patch->xlower;
    ylower = patch->ylower;
    dx = (patch->xupper - patch->xlower)/mx;
    dy = (patch->yupper - patch->ylower)/my;

    double avar[maux];
    int gauge_num = -1;
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d2->fort_gauge_update(&gauge_num, &mx, &my, &mbc, &meqn,
                                        &xlower, &ylower, &dx, &dy,
                                        q, &maux, aux, &xref, &yref,
                                        qvar, avar);
}


/* This is a call to each forward (consumer) patch.
    In this call, we determine which points on the forward 
    patch need to be stored as query points. 
    This routine doesn't know anything about the adjoint problem 
*/

static
void add_patch(fclaw_domain_t * domain, 
               fclaw_patch_t * forward_patch,
               int blockno, int patchno, void *user)
{
    /* assert that user is a valid overlap_forward_metadata_t */
    overlap_forward_metadata_t *c = (overlap_forward_metadata_t *) user;
    FCLAW_ASSERT (c != NULL);
    FCLAW_ASSERT (c->domain != NULL);
    FCLAW_ASSERT (c->query_points != NULL);

    int mx, my, mbc;
    double xlower, ylower, dx, dy;
    fclaw_clawpatch_2d_grid_data (c->glob, forward_patch, &mx, &my, &mbc,
                                 &xlower, &ylower, &dx, &dy);

    double *aux;
    int maux;
    fclaw_clawpatch_aux_data(c->glob,forward_patch,&aux,&maux);

    double* q;
    int meqn;
    fclaw_clawpatch_soln_data(c->glob,forward_patch,&q,&meqn);    

    for (int j = 1; j <= my; j++)
    {
        for (int i = 1; i <= mx; i++)
        {
            /*  Here, we build overlap_point_t for every cell in a patch.
                AMR patches can have different physical sizes, but each patch
                still has mx*my cells. We create one query per local cell so
                overlap interpolation still works for nonmatching patch layouts. */
            overlap_point_t *op;
            op = (overlap_point_t *) sc_array_index(c->query_points,
                                                    c->cell_idx);

            memset(op, -1, sizeof(overlap_point_t));

            op->xy[0]= xlower + (i - 0.5)*dx;
            op->xy[1]= ylower + (j - 0.5)*dy;

            op->i = i;
            op->j = j;
            op->mx = mx;
            op->my = my;
            op->mbc = mbc;
            op->aux_forward = aux;
            op->maux = maux;
            op->q_forward = q;
            op->meqn = meqn;

            op->lnum = c->cell_idx++;
            op->adjoint_metadata.isset = 0;
        }
    }
}

void create_query_points (overlap_forward_metadata_t * c)
{
    /* We create a process-local set of query points, for which we want to
     * obtain interpolation data from the adjoint side. We query the
     * center-point of every local cell; this is what makes the overlap usable
     * after ForestClaw has refined only part of the domain. */
    const fclaw_clawpatch_options_t *clawpatch_opt =
        fclaw_clawpatch_get_options(c->glob);
    int points_per_patch = clawpatch_opt->mx*clawpatch_opt->my;

    c->query_points = sc_array_new_count (sizeof (overlap_point_t),
                                          points_per_patch*
                                          c->domain->local_num_patches);
    c->cell_idx = 0;
    
    fclaw_domain_iterate_patches (c->domain, add_patch, c);

    /* verify that we created as many query_points as expected */
    FCLAW_ASSERT (c->cell_idx ==
                  (size_t) points_per_patch*c->domain->local_num_patches);
}

/* 
    This maps a query point (assumed here to be given in the coordinate 
   system of the adjoint problem) into the native p4est coordinate system.  In this 
   coordinate system, patch coordinates are relative to a [0,1]x[0,1] block 
   coordinate system.

   This is necessary since we may need to determine if a query point is in a non-leaf
   quadrant.
*/   

static
int apply_inverse_adjoint_mapping (overlap_point_t * forward_op, int blockno, 
                                   adjoint_geometry_t * adjoint_geo, double xy[2])
{

    /* Adjoint domain */
    double adjoint_ax = adjoint_geo->fclaw_opt->ax;
    double adjoint_bx = adjoint_geo->fclaw_opt->bx;
    double adjoint_ay = adjoint_geo->fclaw_opt->ay;
    double adjoint_by = adjoint_geo->fclaw_opt->by;


    /* check, if the point lies in the adjoint domain: do this properly */
    if (forward_op->xy[0] < adjoint_ax || forward_op->xy[0] > adjoint_bx ||
        forward_op->xy[1] < adjoint_ay || forward_op->xy[1] > adjoint_by)
    {
        return 0;
    }
    

    /* Assign xy[] point in [0,1]x[0,1] relative to adjoint domain. */
    xy[0] = forward_op->xy[0] - adjoint_ax;
    xy[1] = forward_op->xy[1] - adjoint_ay;
    

    /* We scale from the physical extent in each dimension to the brick extent. */
#if 0    
    xy[0] = xy[0] * (adjoint_fclaw_opt->mi / (adjoint_bx - adjoint_ax));
    xy[1] = xy[1] * (adjoint_fclaw_opt->mj / (adjoint_by - adjoint_ay));
#endif

    /* Adjoint options */
    fclaw_options_t *adjoint_opt = adjoint_geo->fclaw_opt;    

    /* This is now scaled in [0,mi]x[0,mj] */
    xy[0] = (forward_op->xy[0] - adjoint_ax)/(adjoint_bx - adjoint_ax)*adjoint_opt->mi;
    xy[1] = (forward_op->xy[1] - adjoint_ay)/(adjoint_by - adjoint_ay)*adjoint_opt->mj;
   
    // this would be useful when dealing with latlong 3d
    /* The coordinates are now in the [0,mi]x[0,mj] reference coordinate
     * system of the whole brick. Next, we shift xy back to the
     * [0,1]x[0,1] reference system of the block with index blockno on
     *  which we are operating right now. 
     *
     * This is needed so that we can compare xy[] to patch coordinates in the [0,1]x[0,1]
     * p4est patch coordinate system. 
     */
    xy[0] = xy[0] - adjoint_geo->blocks[blockno].vertices[0];
    xy[1] = xy[1] - adjoint_geo->blocks[blockno].vertices[1];

    return 1;                   /* the point lies in the domain */
}


/* Interpolate adjoint data into the representer/forward mesh at the consumer (target)
 * cell centers.  fclaw_overlap_exchange may call this on meta-patches as it
 * searches, but patchno >= 0 means ForestClaw has found the producer (source) leaf patch
 * that owns this query point. */
int overlap_interpolate_patch (fclaw_domain_t *adjoint_domain, 
                               fclaw_patch_t  *adjoint_patch,
                               int blockno, int patchno, 
                               void *point, void *user)
{   
    
    /* assert that we got passed a valid overlap_point_t */
    FCLAW_ASSERT (point != NULL);

    //int init_flag = *(int*) user;


    /* Overlap point in Forward Coordinates */
    overlap_point_t *op = (overlap_point_t *) point;


    /* Assert that we got passed a valid adjoint_geometry_t.
     * We have to pass the fclaw_blocks_t array via the user pointer, because
     * the input domain to this callback is not equal to the forward_domain
     * passed to fclaw_exchange. Whenever the input domain is artificial
     * (domain_is_meta(domain) evaluates to true), domain->blocks is NULL. */
    FCLAW_ASSERT (user != NULL);
    overlap_patch_user_t *u = (overlap_patch_user_t *) user;
    FCLAW_ASSERT(u->geo != NULL);
    FCLAW_ASSERT(u->geo->blocks != NULL);

    /* Apply the inverse mapping of the adjoint side to the point. The result
     * lies in the same reference coordinate system as the patch-boundaries.
     * The inversely mapped point is stored in xy, which we will use for further
     * geometrical operations.
     * If the point lies outside of the domain, we immediately return 0. */
 
    /* Determine if adjoint domain overlaps the current domain.  If it does, the point
       xy[] will contain the location in [0,1]x[0,1] coordinates of  */
    double xy[2];
    int in_adjoint_domain = apply_inverse_adjoint_mapping (op, blockno, u->geo, xy);
    if (!in_adjoint_domain)
    {
        return 0;
    }
    double tol = SC_1000_EPS;

    /* we check if the query point intersects the patch.  Keep in mind that the patch 
      patch coordinates [p->xlower,p->xupper]x[p->ylower,p->yupper] are relative to a
      block [0,1]x[0,1] coordinates.  */
    if ((   xy[0] < adjoint_patch->xlower - tol || xy[0] > adjoint_patch->xupper + tol)
        || (xy[1] < adjoint_patch->ylower - tol || xy[1] > adjoint_patch->yupper + tol))
    {
        /* this IS the actual check for overlapping a point with a patch. */
        return 0;
    }

    /* Check to see if quadrants is a leaf patch */
    if (patchno >= 0)
    {
        /* We have an adjoint leaf patch and so can access all of the 
            data associated with a leaf 
        */
        overlap_adjoint_metadata_t *c = (overlap_adjoint_metadata_t *) u->meta;
        FCLAW_ASSERT (c != NULL);
        FCLAW_ASSERT (c->adjoint_domain != NULL);

        double qvar[op->meqn];
        interpolate_point(c->glob, adjoint_patch,
                          xy[0], xy[1],
                          qvar);
    
        /* For all t, we store the adjoint solution in the aux array for the forward problem.
         This aux array data is used to update the source term. The source term 
         is updated at every time step t>=0  */
        SETAUX_FROM_ADJOINT_POINT(&op->i, &op->j,
                                  &op->mx, &op->my, &op->mbc,
                                  &op->meqn, &op->maux,
                                  qvar, op->aux_forward);
        
        if (u->init_flag != 0)
        {
           //For t=0, we also initialize the representer solution using the adjoint solution
           SETQ_FROM_ADJOINT_POINT(&op->i, &op->j,
                                   &op->mx, &op->my, &op->mbc,
                                   &op->meqn,
                                   qvar, op->q_forward);

        }
                   
        op->adjoint_metadata.isset++;    
    }
     
    return 1;
}

int overlap_interpolate_rep_to_initial_estimate (fclaw_domain_t *representer_domain, 
                               fclaw_patch_t  *representer_patch,
                               int blockno, int patchno, 
                               void *point, void *user)
{   
    /* Same overlap mechanism as overlap_interpolate_patch, but the producer (source) is a
     * representer mesh and the consumer (target) is the model q array.  The consumer (target) q is
     * accumulated with beta_j so the final model update forms sum_j beta_j*r_j. */
    
    /* assert that we got passed a valid overlap_point_t */
    FCLAW_ASSERT (point != NULL);

    //double beta = *((double*) user);


    /* Overlap point in Forward (representer) Coordinates */
    overlap_point_t *op = (overlap_point_t *) point;


    /* Assert that we got passed a valid adjoint_geometry_t.
     * We have to pass the fclaw_blocks_t array via the user pointer, because
     * the input domain to this callback is not equal to the forward_domain
     * passed to fclaw_exchange. Whenever the input domain is artificial
     * (domain_is_meta(domain) evaluates to true), domain->blocks is NULL. */
    FCLAW_ASSERT (user != NULL);
    overlap_patch_user_t *u = (overlap_patch_user_t *) user;
    FCLAW_ASSERT(u->geo != NULL);
    FCLAW_ASSERT(u->geo->blocks != NULL);

    /* Apply the inverse mapping of the adjoint side to the point. The result
     * lies in the same reference coordinate system as the patch-boundaries.
     * The inversely mapped point is stored in xy, which we will use for further
     * geometrical operations.
     * If the point lies outside of the domain, we immediately return 0. */
 
    /* Determine if adjoint domain overlaps the current domain.  If it does, the point
       xy[] will contain the location in [0,1]x[0,1] coordinates of  */
    double xy[2];
    int in_representer_domain = apply_inverse_adjoint_mapping (op, blockno, u->geo, xy);
    if (!in_representer_domain)
    {
        return 0;
    }
    double tol = SC_1000_EPS;

    /* we check if the query point intersects the patch.  Keep in mind that the patch 
      patch coordinates [p->xlower,p->xupper]x[p->ylower,p->yupper] are relative to a
      block [0,1]x[0,1] coordinates.  */
    if ((   xy[0] < representer_patch->xlower - tol || xy[0] > representer_patch->xupper + tol)
        || (xy[1] < representer_patch->ylower - tol || xy[1] > representer_patch->yupper + tol))
    {
        /* this IS the actual check for overlapping a point with a patch. */
        return 0;
    }

    /* Check to see if quadrants is a leaf patch */
    if (patchno >= 0)
    {
        /* We have an adjoint leaf patch and so can access all of the 
            data associated with a leaf 
        */
        overlap_forward_metadata_t *c = (overlap_forward_metadata_t *) u->meta;
        FCLAW_ASSERT (c != NULL);
        FCLAW_ASSERT (c->domain != NULL);

        double qvar[op->meqn];
        interpolate_point(c->glob, representer_patch,
                          xy[0], xy[1],
                          qvar);

        int j = u->rep_idx;
        double b_j = u->beta[j];
    
        UPDATEQ_FROM_REPRESENTER_POINT(&op->i, &op->j,
                                       &op->mx, &op->my, &op->mbc,
                                       &op->meqn,
                                       qvar, op->q_forward, &b_j);
        
                   
        op->adjoint_metadata.isset++;    
    }
     
    return 1;
}

void output_query_points (overlap_forward_metadata_t * c)
{
    /* Printing every query point is overwhelming once AMR uses one query per
     * cell.  Keep a compact diagnostic that confirms the local query count
     * and how many points have already received producer (source)-side data. */
    int npz = c->query_points->elem_count;
    int isset = 0;
    for (int iz = 0; iz < npz; iz++)
    {
        overlap_point_t *op = (overlap_point_t *) sc_array_index (c->query_points, iz);
        if (op->adjoint_metadata.isset)
        {
            isset++;
        }
    }
    fclaw_global_infof("Created %d overlap query points on process %d (%d set).\n",
                       npz, c->domain->mpirank, isset);
}
