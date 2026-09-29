/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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
#include <math.h>

#include "user.h"
#include "af_overlap.h"

#include <fclaw_patch.h>
#include <fclaw_include_all.h>


#include <fclaw_forestclaw.h>
#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_advance.h>
#include <fclaw_regrid.h>
#include <fclaw_output.h>
#include <fclaw_diagnostics.h>
#include <fclaw_vtable.h>
#include <fclaw_context.h>
#include <fclaw_packing.h>

#include "fclaw_math.h"



/*  -----------------------------------------------------------------
    Time stepping
    -- saving time steps
    -- restoring time steps
    -- Time stepping, based on when output files should be created.
    ----------------------------------------------------------------- */

static
void cb_restore_time_step(fclaw_domain_t *domain,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_patch_restore_step(s->glob,this_patch);
}

static
void restore_time_step(fclaw_global_t *glob)
{
    fclaw_global_iterate_patches(glob,cb_restore_time_step,(void *) NULL);

    //fclaw_options_t *fopt = fclaw2d_get_options(glob);
    //fclaw2d_time_sync_reset(glob,fopt->minlevel,fopt->maxlevel,0);
}

static
void cb_save_time_step(fclaw_domain_t *domain,
                       fclaw_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_patch_save_step(s->glob,this_patch);
}

static
void save_time_step(fclaw_global_t *glob)
{
    fclaw_global_iterate_patches(glob,cb_save_time_step,(void *) NULL);
}


/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */


void destroy_first_call(void* user)
{
    int* first_call = (int*) user;
    FCLAW_FREE(first_call);
}

int is_first_call(fclaw_global_t *glob)
{
    int* first_call = (int*) fclaw_global_get_attribute(glob, "fclaw_run_first_call");
    if (first_call == NULL)
    {
        first_call = FCLAW_ALLOC(int, 1);
        *first_call = 0;
        fclaw_global_attribute_store(glob, "fclaw_run_first_call", first_call, NULL, destroy_first_call);
        return 1;
    }
    return *first_call;
}

/* -------------------------------------------------------------------------------
   Output style 1
   Output times are at times [0,dT, 2*dT, 3*dT,...,Tfinal], where dT = tfinal/nout
   -------------------------------------------------------------------------------- */
static
int outstyle_1(fclaw_global_t *glob)
{
    fclaw_domain_t** domain = &glob->domain;
    /* get context object from glob, if there is none, a new one will be created */
    fclaw_context_t *ctx = fclaw_context_get(glob, "fclaw_run_outstyle1_ctx");

    /* true if this is the first call to the run routine for this glob */
    int first_call = is_first_call(glob);

    /* Set error to 0 */
    int init_flag = first_call;  /* Store anything that needs to be stored */
    fclaw_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    int iframe = 0;
    fclaw_context_get_int(ctx, "iframe", &iframe);
    printf("iframe = %d\n", iframe);
    
    double t0 = 0;
    /*fclaw_context_get_double(ctx, "t0", &t0);
    printf("t0 = %f\n", t0);*/
    
    if(iframe == 0 && first_call)
    {
        fclaw_output_frame(glob,iframe);
        /* save context values */
        fclaw_context_save(ctx);
        /* output checkpoint */
        fclaw_output_checkpoint(glob, iframe);
    }

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    double final_time = fclaw_opt->tfinal;
    int nout = fclaw_opt->nout;
    double initial_dt = fclaw_opt->initial_dt;
    int level_factor = pow_int(2,fclaw_opt->maxlevel - fclaw_opt->minlevel);
    double dt_minlevel = initial_dt;
    fclaw_context_get_double(ctx, "dt_minlevel", &dt_minlevel);


    double dt_outer = (final_time-t0)/((double) nout);

    double t_curr = 0;
    fclaw_context_get_double(ctx, "t_curr", &t_curr);
    printf("t_curr = %f\n", t_curr);
    

    int n_inner = 0;
    fclaw_context_get_int(ctx, "n_inner", &n_inner);

    int n = 0;
    fclaw_context_get_int(ctx, "n", &n);
    printf("n = %d\n", n);


    double tend = 0;
    fclaw_context_get_double(ctx, "tend", &tend);

    double dt_step=0, dt_step_desired=0, maxcfl_step=0;
    int took_small_step=0, took_big_step=0;

    fclaw_context_get_double(ctx, "dt_step", &dt_step);
    fclaw_context_get_double(ctx, "dt_step_desired", &dt_step_desired);
    fclaw_context_get_double(ctx, "maxcfl_step", &maxcfl_step);
    fclaw_context_get_int(ctx, "took_small_step", &took_small_step);
    fclaw_context_get_int(ctx, "took_big_step", &took_big_step);

    
    if (!first_call)
    {  
        //fclaw_context_save(ctx);
        goto AFTER_TIMESTEP;
    }

    while(n < nout)
    {
        double tstart = t_curr;

        glob->curr_time = t_curr;
        tend = tstart + dt_outer;

       
        while (t_curr < tend)
        {
            /* In case we have to reject this step */
            if (!fclaw_opt->use_fixed_dt)
            {
                save_time_step(glob);
            }

            /* Use the tolerance to make sure we don't take a tiny time
               step just to hit 'tend'.   We will take a slightly larger
               time step now (dt_cfl + tol) rather than taking a time step
               of 'dt_minlevel' now, followed a time step of only 'tol' in
               the next step.  Of course if 'tend - t_curr > dt_minlevel',
               then dt_minlevel doesn't change. */

            dt_step = dt_minlevel;
            if (fclaw_opt->advance_one_step)
            {
                dt_step /= level_factor;
            }

            double tol = 1e-2*dt_step;
            took_small_step = 0;
            took_big_step = 0;
            dt_step_desired = dt_step;
            if (!fclaw_opt->use_fixed_dt)
            {
                double small_step = tend-(t_curr+dt_step);
                if (small_step  < tol)
                {
                    dt_step = tend - t_curr;  // <= 'dt_minlevel + tol'
                    if (small_step < 0)
                    {
                        /* We have (tend-t_curr) < dt_minlevel, and
                           we have to take a small step to hit tend */
                        took_small_step = 1;
                    }
                    else
                    {
                        /* Take a bigger step now to avoid small step
                           in next time step. */
                        took_big_step = 1;
                    }
                }
            }

            glob->curr_dt = dt_step;  
            maxcfl_step = fclaw_advance_all_levels(glob, t_curr,dt_step);

            if (fclaw_opt->reduce_cfl)
            {
                /* If we are taking a variable time step, we have to reduce the 
                   maxcfl so that every processor takes the same size dt */
                fclaw_timer_start (&glob->timers[FCLAW_TIMER_CFL_COMM]);
                maxcfl_step = fclaw_domain_global_maximum (*domain, maxcfl_step);
                fclaw_timer_stop (&glob->timers[FCLAW_TIMER_CFL_COMM]);                
            }


            double tc = t_curr + dt_step;
            fclaw_global_productionf("Level %d (%d-%d) step %5d : dt = %12.3e; maxcfl (step) = " \
                                     "%16.8f; Final time = %12.4f\n",
                                     fclaw_opt->minlevel,
                                     (*domain)->global_minlevel,
                                     (*domain)->global_maxlevel,
                                     n_inner+1,dt_step,
                                     maxcfl_step, tc);

            if ((maxcfl_step > fclaw_opt->max_cfl) & fclaw_opt->reduce_cfl)
            {
                fclaw_global_essentialf("   WARNING : Maximum CFL exceeded; "    \
                                        "retaking time step\n");

                if (!fclaw_opt->use_fixed_dt)
                {
                    restore_time_step(glob);
                
                    /* Modify dt_level0 from step used. */
                    dt_minlevel = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;

                    /* Got back to start of loop, without incrementing
                       step counter or time level */
                    continue;
                }
            }

            
            /* We are happy with this step */
            n_inner++;
            t_curr += dt_step;
            glob->curr_time = t_curr;

            /* This custom runner returns after one successful step so the
             * caller can exchange overlap data before the next step.  Regrid
             * before returning, otherwise the caller would rebuild overlap
             * on the old mesh and then advance on a newly adapted mesh. */
            if (fclaw_opt->regrid_interval > 0)
            {
                if (n_inner % fclaw_opt->regrid_interval == 0)
                {
                    fclaw_global_infof("regridding at step %d\n",n);
                    fclaw_regrid(glob);
                }
            }

            // Jump here after an overlap exchange.
            fclaw_context_save(ctx);
            return n;  

            AFTER_TIMESTEP:

            /* Update this step, if necessary */
            if (!fclaw_opt->use_fixed_dt)
            {
                double step_fraction = 100.0*dt_step/dt_step_desired;
                if (took_small_step)
                {
                    fclaw_global_infof("   WARNING : Took small time step which was " \
                                       "%6.1f%% of desired dt.\n",
                                       step_fraction);
                }
                if (took_big_step)
                {
                    fclaw_global_infof("   WARNING : Took big time step which was " \
                                       "%6.1f%% of desired dt.\n",step_fraction);

                }


                /* New time step, which should give a cfl close to the
                   desired cfl. */
                double dt_new = dt_minlevel*fclaw_opt->desired_cfl/maxcfl_step;
                if (!took_small_step)
                {
                    dt_minlevel = dt_new;
                }
                else
                {
                    /* use time step that would have been used had we
                       not taken a small step */
                }
            }

            if (fclaw_opt->advance_one_step)
            {
                fclaw_diagnostics_gather(glob, init_flag);                
            }

        }

        /* Output file at every outer loop iteration */
        fclaw_diagnostics_gather(glob, init_flag);
        //glob->curr_time = t_curr;
        iframe++;
        fclaw_output_frame(glob,iframe);

        /* increment n before checkpoint since we want to start at next
           iteration when restarting */
        n++;

        /* save context values */
        fclaw_context_save(ctx);
        /* output checkpoint */
        fclaw_output_checkpoint(glob, iframe);
    }
    //t0 = t_curr;
    
    fclaw_context_save(ctx);

    return n;
}



/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

/*void user_run(fclaw_global_t *forward_glob,\
              fclaw_global_t *adjoint_glob,
              overlap_forward_metadata_t *c,
              adjoint_geometry_t *adjoint_geo)*/
int user_run(fclaw_global_t *glob)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    if (fclaw_opt->outstyle != 1)
    {
        fclaw_global_essentialf("Outstyle %d not implemented yet: Forcing outstyle_1\n",
                                 fclaw_opt->outstyle);
    }

    // run_type : Forward/backward (use enum?) 
    // t_curr_forward : Current forward time - adjoint solver needs to run
    // to this time.
    int n_new = outstyle_1(glob);

    return n_new;
}

#if 0
{
    /* Adjoint setup */
    fclaw_options_t *adjoint_opt = fclaw_get_options(adjoint_glob);

    /* This calls setprob.f90 for the adjoint problem */
    fclaw_problem_setup(adjoint_glob);

    // Run adjoint and produce N+1 restart files
    outstyle_1(adjoint_glob);

    
}
#endif
