/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton, Patricia Azike
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

/* This example demonstrates the use of fclaw_overlap_exchange to exchange
 * interpolation data between two meshes.
 * Both meshes have a clearly assigned role:
 *  forward(renamed as forward) (e.g. Gemini) - queries data for points - represented by forward
 *  producer(renamed as adjoint) (e.g. MAGIC) - provides data - represented by adjoint. */

#include <sys/stat.h>
#include <sys/types.h>


#include "adjoint/adjoint_user.h"
#include "forward/forward_user.h"
#include "forward/model_user.h"
#include "forward/forward_fort.h"
#include "adjoint/adjoint_fort.h"

#include "user.h"

#include "af_overlap.h"

#include <fclaw_global.h>
#include <fclaw_restart.h>
#include <fclaw_output.h>
#include <fclaw_context.h>

#include <fclaw_gauges.h>
#include <fclaw_clawpatch_gauges.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw_diagnostics.h>

#include "fc2d_clawpack46_options.h"

#include <forestclaw2d.h>


typedef struct fclaw_gauge_acc
{
    int dim;
    int num_gauges;
    int num_gauges_set;  /* In case some gauges are not in the domain */    
    fclaw_gauge_t *gauges;
} fclaw_gauge_acc_t;


int m = 1;

static
void setup_overlap(fclaw_global_t *forward_glob,\
                   fclaw_global_t *adjoint_glob,
                   overlap_forward_metadata_t *c,
                   adjoint_geometry_t *adjoint_geo)
{
    
    fclaw_options_t *adjoint_opt = fclaw_get_options(adjoint_glob);

    /* -------------- setup overlap infomation -------------------*/
    /* compute process-local query points on the forward side */

    //overlap_forward_t forward, *c = &forward;
    c->glob = forward_glob;
    c->domain = forward_glob->domain;

    /* AMR breaks the old one-query-per-patch assumption: a forward patch may
     * overlap several adjoint patches or a different refinement level.  Store
     * one query point per cell center so ForestClaw can find the correct leaf
     * patch on the source side for each cell. */
    const fclaw_clawpatch_options_t *forward_clawpatch_opt =
        fclaw_clawpatch_get_options(forward_glob);
    c->num_cells_in_patch =
        forward_clawpatch_opt->mx*forward_clawpatch_opt->my;

    if (c->query_points != NULL)
    {
        sc_array_destroy(c->query_points);
        c->query_points = NULL;
    }  
    create_query_points(c);

    adjoint_geo->fclaw_opt = adjoint_opt;
    adjoint_geo->blocks    = adjoint_glob->domain->blocks;

    /* output the interpolation data for all query points */
    output_query_points (c);
    

}

/* We set up another overlap to map representer solution to model domain.*/
static
void setup_overlap_qhat_update(fclaw_global_t *model_glob,\
                   fclaw_global_t *representer_glob,
                   overlap_forward_metadata_t *c,
                   adjoint_geometry_t *representer_geo)
{
    fclaw_options_t *representer_opt = fclaw_get_options(representer_glob);

    /* -------------- setup overlap infomation -------------------*/
    /* compute process-local query points on the model side */

    //overlap_forward_t model, *c = &model;
    c->glob = model_glob;
    c->domain = model_glob->domain;

    /* Same AMR overlap rule as above: model and representer meshes may no
     * longer have matching patch layouts, so qhat receives one query per model
     * cell center. */
    const fclaw_clawpatch_options_t *model_clawpatch_opt =
        fclaw_clawpatch_get_options(model_glob);
    c->num_cells_in_patch =
        model_clawpatch_opt->mx*model_clawpatch_opt->my;

    if (c->query_points != NULL)
    {
        sc_array_destroy(c->query_points);
        c->query_points = NULL;
    } 
    
    create_query_points(c);

    representer_geo->fclaw_opt = representer_opt;
    representer_geo->blocks    = representer_glob->domain->blocks;

    /* output the interpolation data for all query points */
    output_query_points (c);

}

static
void get_gauges(fclaw_global_t* glob, double *q_out)
{
    fclaw_gauge_acc_t* gauge_acc = 
        (fclaw_gauge_acc_t*) fclaw_diagnostics_get_acc(glob)->gauge_accumulator;


    if (fclaw_diagnostics_get_acc(glob) == NULL)
    {
        printf("get_gauges: diagnostics accumulator is NULL\n");
        FCLAW_ASSERT(0);
    }
    if (gauge_acc == NULL)
    {
        printf("get_gauges: gauge accumulator is NULL\n");
        return; FCLAW_ASSERT(0);
    }

    fclaw_gauge_t *gauges = gauge_acc->gauges;
    if (gauges == NULL || gauge_acc->num_gauges == 0)
    {
        printf("No gauges\n");
        FCLAW_ASSERT(0);
    }

    //fclaw_gauge_t* g = &gauges[0];
    int num_gauges = gauge_acc->num_gauges;

    if (num_gauges != m)
    {
        printf("get_gauges: FATAL -- accumulator has %d gauges but mdata = %d."
               "   gauges.data is stale. Regenerate it\n",
               num_gauges, m);
        exit(1);
    }

    for(int i = 0; i < num_gauges; i++)
    {
        fclaw_gauge_t *g = &gauges[i];

        fclaw_clawpatch_gauge_data_t **gauge_buffer;
        int kmax;

        fclaw_gauges_get_buffer(glob, g, &kmax, (void***) &gauge_buffer);
        if (gauge_buffer == NULL || kmax <= 0)
        {
            printf("No gauge buffer available\n");
            FCLAW_ASSERT(0);
        }

        //for(int k = 0; k < kmax; k++)
        {
            fclaw_clawpatch_gauge_data_t *guser = gauge_buffer[kmax-1];
            if (guser == NULL)
                FCLAW_ASSERT(0);

            int meqn = guser->meqn;
            for(int mq = 0; mq < meqn; mq++)
            {

                q_out[i] = guser->qvar[mq];
            }
            printf("get_gauges: q_out[%d] = %24.16e\n", i, q_out[i]);

            //printf("%24.16e\n",*q_out);
        }

    }

    //fclaw_clawpatch_gauge_data_t *guser = gauge_buffer[0];
    
}


static
void run_program(fclaw_global_t *model_glob,\
                 fclaw_global_t **representer_glob,\
                 fclaw_global_t **adjoint_glob,\
                 overlap_forward_metadata_t *c,
                 adjoint_geometry_t *adjoint_geo)
{

    /*
    We run the model problem to get qF (initial estimate). We also store the value of qF at gauge points
    for use in computing the optimal beta.

    What follows works for a multiple observations
    */
    fclaw_options_t *fclaw_opt_model = fclaw_get_options(model_glob);
    fc2d_clawpack46_options_t *claw_opt = fc2d_clawpack46_get_options(model_glob);
    claw_opt->src_term = 0;

    chdir("model");
    fclaw_opt_model->checkpoint = 1;

    fclaw_output_checkpoint(model_glob, 0);
    
    fclaw_run(model_glob);

    /* Checkpoint model_glob*/
    int iframes = fclaw_opt_model->nout + 1;

    char model_restart_files[iframes][BUFSIZ];
    char model_partition_files[iframes][BUFSIZ];

    for (int i = 0; i < iframes; i++)
    {
        snprintf(model_restart_files[i], BUFSIZ, "fort_frame_%04d.checkpoint", i);
        snprintf(model_partition_files[i], BUFSIZ, "fort_frame_%04d.partition", i);

        
    }
    chdir("..");

    fclaw_opt_model->checkpoint = 0;

    
    /*
        Here is the setup with the checkpoints that enable us to solve for the 
        representer case. We run the adjoint problem and store it in a glob. 
        Then we map the adjoint solution to the representer domain, using the
        adjoint to initialize the representer solution and also in the source
        term. We will be checkpointing both the adjoint and representer solu-
        tions. We checkpoint the adjoint because we need it to restart the 
        representer simulation. On the other hand, we checkpoint the represen-
        ter solution because we need to map it to the model (initial estimate)
        domain for the final part of the computation, which is the evaluation 
        of the optimal estimate, qhat.
    */

    
    /* Declare the variables to checkpoint and restart the representers here for public visibility*/
    char rep_restart_files[iframes][BUFSIZ];
    char rep_partition_files[iframes][BUFSIZ];


    for (int j = 0; j < m; j++)
    {
        int obs_id = j + 1;
        SET_OBSERVATION_INDEX(&obs_id);

        /* Create separate directory for each representer's checkpoints */
        char rep_dir[BUFSIZ];
        snprintf(rep_dir, sizeof(rep_dir), "forward%d", j);
        mkdir(rep_dir, 0755);

        char adj_dir[BUFSIZ];
        snprintf(adj_dir, sizeof(adj_dir), "adjoint%d", j);
        mkdir(adj_dir, 0755);


        /* Remove stale checkpoint files before writing new ones */
        char cmd[4*BUFSIZ];
        snprintf(cmd, sizeof(cmd), "rm -f %s/fort.* %s/fort_frame_*.checkpoint %s/fort_frame_*.partition", adj_dir, adj_dir, adj_dir);
        system(cmd);
        snprintf(cmd, sizeof(cmd), "rm -f %s/fort.* %s/fort_frame_*.checkpoint %s/fort_frame_*.partition", rep_dir, rep_dir, rep_dir);
        system(cmd);

        /* Remove stale fort files in current directory*/
        system("rm -f fort.* ");

        /* Adjoint options :  */
        fclaw_options_t *fclaw_opt_adjoint = fclaw_get_options(adjoint_glob[j]);

        overlap_adjoint_metadata_t p;
        overlap_patch_user_t u1;

        p.glob = adjoint_glob[j];
        p.adjoint_domain = adjoint_glob[j]->domain;

        u1.init_flag = 1;  
        u1.meta = &p;
        u1.geo = adjoint_geo;



        /* Run adjoint once and generate all N+1 restart files */
        fclaw_global_essentialf("-----------------------------------------\n");
        fclaw_global_essentialf("========= Running adjoint simulation =========\n");
        fclaw_global_essentialf("-----------------------------------------\n");
        
        /* 
            When we call fclaw_run(adjoint_glob), ten (10) restart files (from 1 to 10) are generated. 
            However, since Output frame 0 corresponds to restart file 0, we need to get it as well. We 
            can get it by making a call to fclaw_output_checkpoint(adjoint_glob, 0). The catch is that 
            after reinitializing restart file 0, we need to save the context in order to step to the 
            time left for the forward (representer) run to proceed to Tfinal and the adjoint to get to 
            0 from the remaining timestep. 
        */
        
        fclaw_context_t *ctx = fclaw_context_get(adjoint_glob[j], "fclaw_run_outstyle1_ctx");

        fclaw_context_save(ctx);
        
        //chdir("adjoint");
        chdir(adj_dir);
        fclaw_output_checkpoint(adjoint_glob[j], 0); /*checkpoint 0 */
        fclaw_run(adjoint_glob[j]);
        //fclaw_output_checkpoint(adjoint_glob[j], fclaw_opt_adjoint->nout);
        fclaw_global_essentialf("============ Starting representer simulation from adjoint checkpoints =========== \n");
        chdir("..");

        /* Subsequent calls to the adjoint problem should not do any checkpointing or output */
        fclaw_opt_adjoint->nout = 1;
        fclaw_opt_adjoint->checkpoint = 0;

        /* Representer options*/
        fclaw_options_t *fclaw_opt = fclaw_get_options(representer_glob[j]);

        int nout = fclaw_opt->nout;
        double tfinal = fclaw_opt->tfinal;
       
        int n_curr = 0;
        

        /* We are setting the source term to 1 because we need it to solve for the representer*/
        fc2d_clawpack46_options_t *src_opt = fc2d_clawpack46_get_options(representer_glob[j]);
        src_opt->src_term = 1;


        while (n_curr < nout)
        {
            int restart_iframe;

            if (u1.init_flag == 1)
            {
                restart_iframe = nout; 
            }
            else
            {
                restart_iframe = nout - n_curr - 1;
            }

            char restart_file[BUFSIZ];
            char partition_file[BUFSIZ];
            snprintf(restart_file, BUFSIZ, "fort_frame_%04d.checkpoint", restart_iframe);
            snprintf(partition_file, BUFSIZ, "fort_frame_%04d.partition", restart_iframe);

            fclaw_global_essentialf("-----------------------------------------\n");
            fclaw_global_essentialf("Restarting from checkpoint file %s\n", restart_file);
            fclaw_global_essentialf("-----------------------------------------\n");

            /* Reinitialize adjoint*/
            //chdir("adjoint");
            chdir(adj_dir);
            
            fclaw_reinitialize_from_file(adjoint_glob[j], restart_file, partition_file);
            chdir("..");

            /*setup overlap infomation*/
            setup_overlap(representer_glob[j], adjoint_glob[j],&c[j], adjoint_geo);
            
            if (u1.init_flag == 1)
            {
                /* Map adjoint solution to the representer domain using overlap mechanism */
                fclaw_overlap_exchange(adjoint_glob[j]->domain, c[j].query_points,
                                    overlap_interpolate_patch, &u1);
                u1.init_flag = 0;
                continue;
            }

            /* Checkpoint and restart files for representers*/

            for(int i = 0; i < iframes; i++)
            {
                snprintf(rep_restart_files[i], BUFSIZ, "fort_frame_%04d.checkpoint", i);
                snprintf(rep_partition_files[i], BUFSIZ, "fort_frame_%04d.partition", i);
            }

            // Exit after one step
            int notdone = 1;
            while (notdone)
            {
                // Take one step of representer problem
                printf(">>>            representer run: t = %12.4f\n", representer_glob[j]->curr_time);

                /* Change into forward directory to store restart files there*/
                //chdir("forward");
                chdir(rep_dir);

                fclaw_opt->checkpoint = 1;

                int n_new = user_run(representer_glob[j]);
                
                chdir("..");

                fclaw_opt->checkpoint = 0;

                double t_curr_representer = representer_glob[j]->curr_time;

                //chdir("adjoint");
                chdir(adj_dir);
                fclaw_reinitialize_from_file(adjoint_glob[j], restart_file, partition_file);
                chdir("..");

                setup_overlap(representer_glob[j], adjoint_glob[j],&c[j], adjoint_geo);

                fclaw_opt_adjoint->tfinal = tfinal - t_curr_representer;
                printf(">>>            ADJOINT run: tfinal set to %12.4f\n", fclaw_opt_adjoint->tfinal);
                fclaw_run(adjoint_glob[j]);


                /* Map adjoint solution to the representer domain using overlap mechanism */
                fclaw_overlap_exchange(adjoint_glob[j]->domain, c[j].query_points,
                                       overlap_interpolate_patch, &u1);

                if (n_new > n_curr)
                {
                    n_curr = n_new;
                    notdone = 0;
                }
            }
        }
    }

#ifdef INSTANTANEOUS_VERIFICATION_ONLY
    fclaw_global_essentialf("Instantaneous verification complete; "
                            "skipping the coefficient-system solve.\n");
    for (int j = 0; j < m; j++)
    {
        adjoint_finalize(adjoint_glob[j]);
        forward_finalize(representer_glob[j]);
    }
    model_finalize(model_glob);
    return;
#endif
    
    /*Calculate beta using gauge value from model (initial estimate) and representer*/

    double *qF_gauge = new double[m];
    double *h  = new double[m];
    double *A  = new double[m*m];

    /* 1. Get Forward model (background trajectory) at all observation locations */
    printf("Calling get_gauges for model\n");
    get_gauges(model_glob, qF_gauge);
        

    /* 2. Build h = d - qF */
    for (int i = 0; i < m; i++)
    {
        const adjoint_options_t *adj_opt = adjoint_get_options(adjoint_glob[0]);
        h[i] = adj_opt->dm[i] - qF_gauge[i];
    }

    /* 3. Build representer matrix */
    for (int j = 0; j < m; j++)
    {

        double *Rj = new double[m];

        printf("Calling get_gauges for representer j = %d\n", j);
        get_gauges(representer_glob[j], Rj);

        for (int i = 0; i < m; i++)
        {
            A[i + j*m] = Rj[i];
        }

        delete[] Rj;
    }


    /* --- Representer-matrix symmetry: forward/adjoint consistency check --- 
        We have a small asymmetry*/ 
    double rmax = 0.0, asym = 0.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            double aij = A[i + j*m], aji = A[j + i*m];
            if (fabs(aij)       > rmax) rmax = fabs(aij);
            if (fabs(aij - aji) > asym) asym = fabs(aij - aji);
        }
    printf("Representer matrix: max|R_ij - R_ji| = %.3e, "
           "relative asymmetry = %.3e\n", asym, asym / rmax);

    double fro_diff = 0.0, fro_R = 0.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            double aij = A[i + j*m], aji = A[j + i*m];
            fro_diff += (aij - aji) * (aij - aji);
            fro_R    += aij * aij;
        }
    printf("frobenius norm: Representer matrix: ||R - R^T||_F / ||R||_F = %.3e\n",
           sqrt(fro_diff) / sqrt(fro_R));
    /* --------------------------------------------------------------------- */
    /* We remove the asymmetry here*/
    
    for (int i = 0; i < m; i++)
    for (int j = i + 1; j < m; j++)
    {
        double avg = 0.5 * (A[i + j*m] + A[j + i*m]);
        A[i + j*m] = avg;
        A[j + i*m] = avg;
    }


    /* 4. Add W^{-1} to diagonal */
    const adjoint_options_t *adj_opt = adjoint_get_options(adjoint_glob[0]);
    for (int i = 0; i < m; i++)
    {
        A[i + i*m] += 1.0 / adj_opt->W_eps[i];
    }

    /* Solve (R + W^-1)β = h using dgesv solver */

    int n    = m;
    int nrhs = 1;
    int lda  = m;
    int ldb  = m;
    int info;

    int *ipiv = new int[m];

    dgesv_(&n, &nrhs, A, &lda, ipiv, h, &ldb, &info);

    /* Save beta separately so h can be freed */
    double *beta = new double[m];
    if (info == 0)
    {
        for (int i = 0; i < m; i++)
        {
            beta[i] = h[i];
            printf("beta[%d] = %24.16e\n", i, beta[i]);
        }
    }
    else
    {
        printf("dgesv failed; info = %d\n", info);
    }

        
        /* Map representer solution to model domain */
        /*For n goes from 0 to nout, load checkpoint files from model glob and representer glob*/
    #if 1
        
    
    //fclaw_options_t *fclaw_opt = fclaw_get_options(representer_glob[j]);
    int nout = fclaw_opt_model->nout;

    

    for (int iframe=0; iframe<=nout; iframe++)
    {

        /* Load checkpoint n from model glob get initial estimate q_F */
        chdir("model");
        fclaw_global_essentialf("-----------------------------------------\n");
        fclaw_global_essentialf("Restarting model from checkpoint file %s\n", model_restart_files[iframe]);
        fclaw_global_essentialf("-----------------------------------------\n");

        fclaw_reinitialize_from_file(model_glob, model_restart_files[iframe], model_partition_files[iframe]);
        chdir("..");

        for (int j = 0; j < m; j++)
        {
            int obs_id = j + 1;
            SET_OBSERVATION_INDEX(&obs_id);


            /* Load checkpoint n from representer (representer glob)*/
            //chdir("forward");
            char rep_dir[BUFSIZ];
            snprintf(rep_dir, sizeof(rep_dir), "forward%d", j);
            chdir(rep_dir);
            fclaw_global_essentialf("-----------------------------------------\n");
            fclaw_global_essentialf("Restarting representer from checkpoint file %s\n", rep_restart_files[iframe]);
            fclaw_global_essentialf("-----------------------------------------\n");

            fclaw_reinitialize_from_file(representer_glob[j], rep_restart_files[iframe], rep_partition_files[iframe]);
            chdir("..");

            setup_overlap_qhat_update(model_glob, representer_glob[j],&c[j], adjoint_geo);

            overlap_patch_user_t u2;
            u2.meta = &c[j];
            u2.geo = adjoint_geo;
            u2.rep_idx = j;
            u2.beta = beta;

            /* Map representer solution to the model (initial estimate) domain using overlap mechanism */
            fclaw_overlap_exchange(representer_glob[j]->domain, c[j].query_points,
                                   overlap_interpolate_rep_to_initial_estimate, &u2);
        }
        /* Call output frame on updated model (model_q, that is, q_hat)*/
        fclaw_output_frame (model_glob, iframe);
    
    }
#endif

        /* Finalize solvers */
    for (int j = 0; j < m; j++)
    {    
        adjoint_finalize(adjoint_glob[j]);
        
        forward_finalize(representer_glob[j]); 
    }
    model_finalize(model_glob);
          
}


int
main (int argc, char **argv)
{
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);

    /* Register packages */
    //Global options like verbosity, etc
    fclaw_app_options_register_core(app, "forward_options.ini"); 

    /* Register once to read from ini file */
    fclaw_options_t *adjoint_fclaw_opt_base =
        fclaw_options_register(app, "adjoint", "adjoint_options.ini");
    fclaw_clawpatch_options_t *adjoint_clawpatch_opt_base =
        fclaw_clawpatch_2d_options_register(app, "adjoint-clawpatch", "adjoint_options.ini");
    fc2d_clawpack46_options_t *adjoint_claw46_opt_base =
        fc2d_clawpack46_options_register(app, "adjoint-clawpack46", "adjoint_options.ini");
    adjoint_options_t *adjoint_user_opt_base =
        adjoint_options_register(app, "adjoint_options.ini");

    /* representer options - register once to read from ini file */
    fclaw_options_t           *representer_fclaw_opt_base =
        fclaw_options_register(app, "forward", "forward_options.ini");
    fclaw_clawpatch_options_t *representer_clawpatch_opt_base =
        fclaw_clawpatch_2d_options_register(app, "forward-clawpatch", "forward_options.ini");
    fc2d_clawpack46_options_t *representer_claw46_opt_base =
        fc2d_clawpack46_options_register(app, "forward-clawpack46", "forward_options.ini");
    forward_options_t         *representer_user_opt_base =
        forward_options_register(app, "forward_options.ini");

    

    /* model (initial estimate) options */
    model_options_t             *model_user_opt;
    fclaw_options_t             *model_fclaw_opt;
    fclaw_clawpatch_options_t *model_clawpatch_opt;
    fc2d_clawpack46_options_t   *model_claw46_opt;

    model_fclaw_opt =                   
        fclaw_options_register(app, "model","model_options.ini");
    model_clawpatch_opt =   
        fclaw_clawpatch_2d_options_register(app, "model-clawpatch", "model_options.ini");
    model_claw46_opt =        
        fc2d_clawpack46_options_register(app, "model-clawpack46","model_options.ini");
    model_user_opt = model_options_register(app, "model_options.ini");  



    /* Read configuration file(s) */
    int first_arg;
    fclaw_exit_type_t vexit = 
        fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");
        
    if (!vexit)
    {
        /* Options have been checked and are valid */
        //int size, rank;
        m = adjoint_user_opt_base->mdata;

        /* Copy options for each j so they have independent memory */

        /*
        When registering options with the same section name, fclaw_options_register
        returns the same pointer for all j, causing shared state between globs.
        Modifications to one glob's options (e.g. nout=1, checkpoint=0) were
        inadvertently affecting all other globs, causing missing checkpoint files
        and crashes for mdata > 2.

        Fix by registering options once to correctly read from the ini file, then
        memcpy-ing independent copies for each j so that per-glob modifications
        are isolated.
        */
        fclaw_options_t           *adjoint_fclaw_opt[m];
        fclaw_clawpatch_options_t *adjoint_clawpatch_opt[m];
        fc2d_clawpack46_options_t *adjoint_claw46_opt[m];
        adjoint_options_t         *adjoint_user_opt[m];

        for (int j = 0; j < m; j++)
        {
            adjoint_fclaw_opt[j]     = new fclaw_options_t;
            adjoint_clawpatch_opt[j] = new fclaw_clawpatch_options_t;
            adjoint_claw46_opt[j]    = new fc2d_clawpack46_options_t;
            adjoint_user_opt[j]      = new adjoint_options_t;

            memcpy(adjoint_fclaw_opt[j],     adjoint_fclaw_opt_base,     sizeof(fclaw_options_t));
            memcpy(adjoint_clawpatch_opt[j], adjoint_clawpatch_opt_base, sizeof(fclaw_clawpatch_options_t));
            memcpy(adjoint_claw46_opt[j],    adjoint_claw46_opt_base,    sizeof(fc2d_clawpack46_options_t));
            memcpy(adjoint_user_opt[j],      adjoint_user_opt_base,      sizeof(adjoint_options_t));
        }

        
        fclaw_global_t **adjoint_glob = new fclaw_global_t*[m];

        for (int j = 0; j < m; j++)
        {
            adjoint_glob[j] = fclaw_global_new(app);

            fclaw_options_store            (adjoint_glob[j], adjoint_fclaw_opt[j]);
            fclaw_clawpatch_options_store  (adjoint_glob[j], adjoint_clawpatch_opt[j]);
            fc2d_clawpack46_options_store    (adjoint_glob[j], adjoint_claw46_opt[j]);
            adjoint_options_store           (adjoint_glob[j], adjoint_user_opt[j]);

            adjoint_create_domain(adjoint_glob[j]);
        }

        /* Copy representer options for each j */
        fclaw_options_t           *representer_fclaw_opt[m];
        fclaw_clawpatch_options_t *representer_clawpatch_opt[m];
        fc2d_clawpack46_options_t *representer_claw46_opt[m];
        forward_options_t         *representer_user_opt[m];

        for (int j = 0; j < m; j++)
        {
            representer_fclaw_opt[j]     = new fclaw_options_t;
            representer_clawpatch_opt[j] = new fclaw_clawpatch_options_t;
            representer_claw46_opt[j]    = new fc2d_clawpack46_options_t;
            representer_user_opt[j]      = new forward_options_t;

            memcpy(representer_fclaw_opt[j],     representer_fclaw_opt_base,     sizeof(fclaw_options_t));
            memcpy(representer_clawpatch_opt[j], representer_clawpatch_opt_base, sizeof(fclaw_clawpatch_options_t));
            memcpy(representer_claw46_opt[j],    representer_claw46_opt_base,    sizeof(fc2d_clawpack46_options_t));
            memcpy(representer_user_opt[j],      representer_user_opt_base,      sizeof(forward_options_t));
        }

        /* representer setup */
        fclaw_global_t **representer_glob = new fclaw_global_t*[m];
        for (int j = 0; j < m; j++)
        {
            representer_glob[j] = fclaw_global_new(app);

            fclaw_options_store           (representer_glob[j], representer_fclaw_opt[j]);
            fclaw_clawpatch_options_store (representer_glob[j], representer_clawpatch_opt[j]);
            fc2d_clawpack46_options_store   (representer_glob[j], representer_claw46_opt[j]);
            forward_options_store           (representer_glob[j], representer_user_opt[j]);

            forward_create_domain(representer_glob[j]);
        }
        
        /* Create model glob*/
        fclaw_global_t *model_glob = fclaw_global_new(app);

        fclaw_options_store           (model_glob, model_fclaw_opt);
        fclaw_clawpatch_options_store (model_glob, model_clawpatch_opt);
        fc2d_clawpack46_options_store   (model_glob, model_claw46_opt);
        model_options_store           (model_glob,model_user_opt);

        forward_create_domain(model_glob);


        /* initialize all the solvers before doing overlap */
        
        for (int j = 0; j < m; j++)
        {
            adjoint_initialize(adjoint_glob[j]);
            forward_initialize(representer_glob[j]); 
        }
        model_initialize(model_glob); 

        /* Set up overlap points in "c" */
        overlap_forward_metadata_t *c = new overlap_forward_metadata_t[m];

        for (int j = 0; j < m; j++)
        {
            memset(&c[j], 0, sizeof(overlap_forward_metadata_t));
        }

        /* Geometry of the adjoint problem */
        adjoint_geometry_t adjoint_geometry, *adjoint_geo = &adjoint_geometry;        

        /* initialize the adjoint geometry information that is needed for
            * mapping between the representer and the adjoint domain */
        for (int j = 0; j < m; j++)
        {
            setup_overlap(representer_glob[j], adjoint_glob[j],&c[j], adjoint_geo);

        }
        /* Run the program */
        run_program(model_glob,
                    representer_glob, 
                    adjoint_glob,
                    c,adjoint_geo);    

        for (int j = 0; j < m; j++)
        {
            /* destroy overlap points */
            if (c[j].query_points != NULL)
            {
                sc_array_destroy(c[j].query_points);
                c[j].query_points = NULL;
            }
        }    

        for (int j = 0; j < m; j++)
        {
            fclaw_global_destroy(adjoint_glob[j]);
            fclaw_global_destroy(representer_glob[j]);
        }
        fclaw_global_destroy(model_glob);

        /* Free copied options */
        for (int j = 0; j < m; j++)
        {
            delete adjoint_fclaw_opt[j];
            delete adjoint_clawpatch_opt[j];
            delete adjoint_claw46_opt[j];
            delete adjoint_user_opt[j];

            delete representer_fclaw_opt[j];
            delete representer_clawpatch_opt[j];
            delete representer_claw46_opt[j];
            delete representer_user_opt[j];
        }

        delete[] adjoint_glob;
        delete[] representer_glob;
    }

    fclaw_app_destroy (app);

    return 0;
}
