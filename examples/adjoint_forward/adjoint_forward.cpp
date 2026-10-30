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

/* Representer method (weak-constraint 4D-Var) for idealized smoke transport


   1. model         : prior run q_F (checkpointed at every output frame)
   2. for each observation j
        adjoint     : alpha_j, backward in time, forced by observation j
                      (checkpointed at every output frame)
        representer : r_j, forward in time, with r_j(0) = alpha_j(0)/W_i and
                      source alpha_j(t)/W_f.  Before each representer step
                      the adjoint is restarted from the nearest checkpoint
                      and advanced to tau = T - t.
   3. solve (R + W_eps^{-1}) beta = d - q_F, R_ij = r_j(x_i, y_i, t_i). Here 
      R is R_sym, the symmetrized R
   4. q_hat = q_F + sum_j beta_j r_j, written to pwd at every frame.

*/


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



/* ----------------------------- Directories --------------------------------- */

static
void change_dir(const char *dir)
{
    if (chdir(dir) != 0)
    {
        fclaw_abortf("adjoint_forward: could not change to directory %s\n", dir);
    }
}

/* Create dir if needed and remove output and checkpoints left by an earlier run */
static
void make_dir(const char *dir)
{
    mkdir(dir, 0755);

    char cmd[4*BUFSIZ];
    snprintf(cmd, sizeof(cmd),
             "rm -f %s/fort.* %s/fort_frame_*.checkpoint %s/fort_frame_*.partition",
             dir, dir, dir);
    if (system(cmd) != 0)
    {
        fclaw_global_essentialf("adjoint_forward: could not clear directory %s\n", dir);
    }
}

/* Restore glob from output frame iframe of the checkpoints written in dir */
static
void reinitialize(fclaw_global_t *glob, const char *dir, int iframe)
{
    char restart_file[BUFSIZ];
    char partition_file[BUFSIZ];
    snprintf(restart_file, BUFSIZ, "fort_frame_%04d.checkpoint", iframe);
    snprintf(partition_file, BUFSIZ, "fort_frame_%04d.partition", iframe);

    change_dir(dir);
    fclaw_reinitialize_from_file(glob, restart_file, partition_file);
    change_dir("..");
}

/* ------------------------------- Overlap ----------------------------------- */

/* Query points (one per cell centre) on the consumer mesh, and the geometry of
   the producer, the mesh that has the data.  Used both ways: adjoint (P) ->
   representer (c) and representer (P) -> model (C). */
static
void setup_overlap(fclaw_global_t *consumer_glob,
                   fclaw_global_t *producer_glob,
                   overlap_forward_metadata_t *c,
                   adjoint_geometry_t *producer_geo)
{
    c->glob = consumer_glob;
    c->domain = consumer_glob->domain;

    /* AMR breaks the old one-query-per-patch assumption: a consumer patch may
     * overlap several producer patches or a different refinement level.  Store
     * one query point per cell center so ForestClaw can find the correct leaf
     * patch on the producer side for each cell. */
    const fclaw_clawpatch_options_t *consumer_clawpatch_opt =
        fclaw_clawpatch_get_options(consumer_glob);
    c->num_cells_in_patch =
        consumer_clawpatch_opt->mx*consumer_clawpatch_opt->my;

    if (c->query_points != NULL)
    {
        sc_array_destroy(c->query_points);
        c->query_points = NULL;
    }
    create_query_points(c);

    producer_geo->fclaw_opt = fclaw_get_options(producer_glob);
    producer_geo->blocks    = producer_glob->domain->blocks;

    /* output the interpolation data for all query points */
    output_query_points(c);
}

/* alpha -> representer aux (and, while u->init_flag = 1, the representer itself) */
static
void map_adjoint_to_representer(fclaw_global_t *adjoint_glob,
                                fclaw_global_t *representer_glob,
                                overlap_forward_metadata_t *c,
                                adjoint_geometry_t *geo,
                                overlap_patch_user_t *u)
{
    setup_overlap(representer_glob, adjoint_glob, c, geo);
    fclaw_overlap_exchange(adjoint_glob->domain, c->query_points,
                           overlap_interpolate_patch, u);
}

/* q_model <- q_model + beta[j] r_j */
static
void map_representer_to_model(fclaw_global_t *representer_glob,
                              fclaw_global_t *model_glob,
                              overlap_forward_metadata_t *c,
                              adjoint_geometry_t *geo,
                              int j, double *beta)
{
    setup_overlap(model_glob, representer_glob, c, geo);

    overlap_patch_user_t u;
    u.meta = c;
    u.geo = geo;
    u.rep_idx = j;
    u.beta = beta;
    fclaw_overlap_exchange(representer_glob->domain, c->query_points,
                           overlap_interpolate_rep_to_initial_estimate, &u);
}

/* -------------------------------- Gauges ----------------------------------- */

/* q_out[i] = the last value recorded at gauge i.  Returns the number of gauges:
   0 if the run has none (output-gauges = F, or no gauges in gauges.data), and
   otherwise mdata, one gauge per observation. */
static
int get_gauges(fclaw_global_t* glob, double *q_out, int mdata)
{
    fclaw_diagnostics_accumulator_t *acc = fclaw_diagnostics_get_acc(glob);
    fclaw_gauge_acc_t *gauge_acc =
        (acc == NULL) ? NULL : (fclaw_gauge_acc_t*) acc->gauge_accumulator;

    if (gauge_acc == NULL || gauge_acc->gauges == NULL || gauge_acc->num_gauges == 0)
    {
        return 0;
    }

    int num_gauges = gauge_acc->num_gauges;
    if (num_gauges != mdata)
    {
        fclaw_abortf("get_gauges: the run has %d gauges but mdata = %d.  "
                     "gauges.data is stale: regenerate it.\n", num_gauges, mdata);
    }

    for (int i = 0; i < num_gauges; i++)
    {
        fclaw_gauge_t *g = &gauge_acc->gauges[i];

        fclaw_clawpatch_gauge_data_t **gauge_buffer;
        int kmax;
        fclaw_gauges_get_buffer(glob, g, &kmax, (void***) &gauge_buffer);

        /* The last value the gauge recorded */
        fclaw_clawpatch_gauge_data_t *guser = NULL;
        if (gauge_buffer != NULL && kmax > 0)
        {
            guser = gauge_buffer[kmax-1];
        }
        if (guser == NULL)
        {
            fclaw_abortf("get_gauges: gauge %d (observation %d) recorded no value; "
                         "check its time window in gauges.data against tm.\n",
                         i, i + 1);
        }

        FCLAW_ASSERT(guser->meqn == 1);
        q_out[i] = guser->qvar[0];
        printf("get_gauges: q_out[%d] = %24.16e\n", i, q_out[i]);
    }
    return num_gauges;
}

/* --------------------------------- Model ----------------------------------- */

/* Run the model (the prior: no source term) to T in dir, with a checkpoint at
   every output frame.  model_glob is not advanced again: compute_beta reads
   q_F at the gauges from it, and the q_hat stage restores it from the
   checkpoints. */
static
void run_model(fclaw_global_t *model_glob, const char *dir)
{
    fclaw_options_t *fclaw_opt = fclaw_get_options(model_glob);
    fc2d_clawpack46_options_t *claw_opt = fc2d_clawpack46_get_options(model_glob);
    claw_opt->src_term = 0;

    change_dir(dir);
    fclaw_opt->checkpoint = 1;
    fclaw_output_checkpoint(model_glob, 0);
    fclaw_run(model_glob);
    fclaw_opt->checkpoint = 0;
    change_dir("..");
}

/* ------------------------------ Coefficients ------------------------------- */

/* Solve (R_sym + W_eps^{-1}) beta = d - q_F, where q_F and R(i,j) = r_j are
   the model and the representers at the gauges, and R_sym = (R + R^T)/2.
   Returns 0, and leaves beta alone, if the runs have no gauges. */
static
int compute_beta(fclaw_global_t *model_glob,
                 fclaw_global_t **representer_glob,
                 const adjoint_options_t *adj_opt,
                 int mdata, double beta[])
{
    /* 1. The model (prior) at the observation locations */
    double *qF = new double[mdata];
    printf("Calling get_gauges for model\n");
    if (get_gauges(model_glob, qF, mdata) == 0)
    {
        delete[] qF;
        return 0;
    }

    /* 2. The representers at the observation locations:
          R_ij = r_j at gauge i, stored column by column in A */
    double *Rj = new double[mdata];
    double *A  = new double[mdata*mdata];

    for (int j = 0; j < mdata; j++)
    {
        printf("Calling get_gauges for representer %d (observation %d)\n", j, j + 1);
        if (get_gauges(representer_glob[j], Rj, mdata) == 0)
        {
            fclaw_abortf("compute_beta: the model has gauges but representer %d has "
                         "none (check output-gauges in forward_options.ini).\n", j);
        }
        for (int i = 0; i < mdata; i++)
        {
            A[i + j*mdata] = Rj[i];
        }
    }
    delete[] Rj;

    /* --- Representer-matrix symmetry: forward/adjoint consistency check ---
        We have a small asymmetry*/
    double rmax = 0.0;
    double asym = 0.0;
    double fro_diff = 0.0;
    double fro_R = 0.0;

    for (int i = 0; i < mdata; i++)
        for (int j = 0; j < mdata; j++)
        {
            double aij = A[i + j*mdata];
            double aji = A[j + i*mdata];
            
            if (fabs(aij)       > rmax) rmax = fabs(aij);
            if (fabs(aij - aji) > asym) asym = fabs(aij - aji);
            fro_diff += (aij - aji) * (aij - aji);
            fro_R    += aij * aij;
        }
    printf("Representer matrix: max|R_ij - R_ji| = %.3e, "
           "relative asymmetry = %.3e\n", asym, asym / rmax);
    printf("frobenius norm: Representer matrix: ||R - R^T||_F / ||R||_F = %.3e\n",
           sqrt(fro_diff) / sqrt(fro_R));

    /* We remove the asymmetry here */
    for (int i = 0; i < mdata; i++)
        for (int j = i + 1; j < mdata; j++)
        {
            double avg = 0.5 * (A[i + j*mdata] + A[j + i*mdata]);
            A[i + j*mdata] = avg;
            A[j + i*mdata] = avg;
        }

    /* Add W^{-1} to the diagonal; right-hand side h = d - q_F */
    double *h = new double[mdata];
    for (int i = 0; i < mdata; i++)
    {
        A[i + i*mdata] += 1.0 / adj_opt->W_eps[i];
        h[i] = adj_opt->dm[i] - qF[i];
    }

    /* Solve (R + W^-1) beta = h using dgesv; it overwrites h with the solution */
    int n    = mdata;
    int nrhs = 1;
    int lda  = mdata;
    int ldb  = mdata;
    int info;
    int *ipiv = new int[mdata];

    dgesv_(&n, &nrhs, A, &lda, ipiv, h, &ldb, &info);
    if (info != 0)
    {
        fclaw_abortf("compute_beta: dgesv failed; info = %d\n", info);
    }
    for (int i = 0; i < mdata; i++)
    {
        beta[i] = h[i];
        printf("beta[%d] = %24.16e\n", i, beta[i]);
    }

    delete[] h;
    delete[] ipiv;
    delete[] A;
    delete[] qF;
    return 1;
}

/* ------------------------------- Program ----------------------------------- */

static
void run_program(fclaw_global_t *model_glob,
                 fclaw_global_t **representer_glob,
                 fclaw_global_t **adjoint_glob,
                 overlap_forward_metadata_t *c,
                 adjoint_geometry_t *geo,
                 int mdata)
{
    /* Optimal estimate written by an earlier run */
    if (system("rm -f fort.*") != 0)
    {
        fclaw_global_essentialf("adjoint_forward: could not remove old fort.* files\n");
    }

    /* ---------------- 1. Model (prior q_F), checkpointed ------------------ */
    run_model(model_glob, "model");

    /*
        Here is the setup with the checkpoints that enable us to solve for the
        representer case. We run the adjoint problem and store it in a glob.
        Then we map the adjoint solution to the representer domain, using the
        adjoint to initialize the representer solution and also in the source
        term. We will be checkpointing both the adjoint and representer solu-
        tions. We checkpoint the adjoint because we need it to restart the
        representer simulation. On the other hand, we checkpoint the represen-
        ter solution because we need to map it to the model (prior) domain for 
        the final part of the computation, which is the evaluation of the 
        optimal estimate, qhat.
    */

    /* ------------- 2. Adjoint and representer, per observation ------------ */
    for (int j = 0; j < mdata; j++)
    {
        /* Observation obs_id = j + 1 (1-based, as in the ini file and the
           Fortran): its adjoint is adjoint_glob[j], run in adjoint<j>/, and
           its representer is representer_glob[j], run in forward<j>/. */
        int obs_id = j + 1;
        SET_OBSERVATION_INDEX(&obs_id);

        char rep_dir[BUFSIZ];
        snprintf(rep_dir, sizeof(rep_dir), "forward%d", j);
        make_dir(rep_dir);

        char adj_dir[BUFSIZ];
        snprintf(adj_dir, sizeof(adj_dir), "adjoint%d", j);
        make_dir(adj_dir);

        /* Adjoint options :  */
        fclaw_options_t *fclaw_opt_adjoint = fclaw_get_options(adjoint_glob[j]);

        overlap_adjoint_metadata_t p;
        overlap_patch_user_t u1;

        p.glob = adjoint_glob[j];
        p.adjoint_domain = adjoint_glob[j]->domain;

        u1.init_flag = 1;
        u1.meta = &p;
        u1.geo = geo;

        /* Run adjoint once and generate all N+1 restart files */
        fclaw_global_essentialf("-----------------------------------------\n");
        fclaw_global_essentialf("========= Running adjoint simulation: observation %d "
                                "(adjoint%d) =========\n", obs_id, j);
        fclaw_global_essentialf("-----------------------------------------\n");

        /*
            Adjoint checkpoints for frames 0..nout. When we call fclaw_run, we get frames 
            1..nout, and each one stores fclaw_run's context (frame number n, time t_curr, ...), 
            so restarting from frame k resumes fclaw_run at tau_k.  Since output frame 0 
            corresponds to restart file 0, we need it as well. We write it with 
            fclaw_output_checkpoint(adjoint_glob, 0). It is used in the last interval of 
            the representer run, t in (T - dT, T], where the adjoint restarts from frame 0 
            and steps through the time left, tau = T - t.  The catch is that frame 0 is written
            before fclaw_run has created a context, so we get and save an empty one first; 
            restarting from frame 0, fclaw_run then finds no stored values and starts from 
            n = 0, tau = 0.
        */

        fclaw_context_t *ctx = fclaw_context_get(adjoint_glob[j], "fclaw_run_outstyle1_ctx");

        fclaw_context_save(ctx);

        change_dir(adj_dir);
        fclaw_output_checkpoint(adjoint_glob[j], 0); /*checkpoint 0 */
        fclaw_run(adjoint_glob[j]);
        fclaw_global_essentialf("============ Starting representer simulation: observation %d "
                                "(forward%d), from adjoint%d checkpoints ===========\n",
                                obs_id, j, j);
        change_dir("..");

        /* Representer options*/
        fclaw_options_t *fclaw_opt = fclaw_get_options(representer_glob[j]);

        int nout = fclaw_opt->nout;
        double tfinal = fclaw_opt->tfinal;

        /* The adjoint is restarted from its output frames, so they must be
           the representer's output frames: same tfinal and nout. */
        if (fclaw_opt_adjoint->nout != nout ||
            fabs(fclaw_opt_adjoint->tfinal - tfinal) > 1e-12*tfinal)
        {
            fclaw_global_essentialf("adjoint_forward : the adjoint and the representer "
                                    "need the same tfinal and nout\n");
            exit(1);
        }

        /* Frame spacing of the adjoint, tau_k = k*dT. It uses the representer's tfinal
           and nout, which stay fixed. Those of the adjoint are overwritten by each
           re-advance below. 
        */
        const double dT = tfinal/nout;

        /* Re-advances of the adjoint write no checkpoints and no output (the
           adjoint frames from the full run above must not be overwritten). */
        fclaw_opt_adjoint->checkpoint = 0;
        fclaw_opt_adjoint->output = 0;

        int n_curr = 0;

        /* We are setting the source term to 1 because we need it to solve for the representer*/
        fc2d_clawpack46_options_t *src_opt = fc2d_clawpack46_get_options(representer_glob[j]);
        src_opt->src_term = 1;

        /* t = 0 : the representer starts from alpha at tau = T (frame nout) */
        fclaw_global_essentialf("-----------------------------------------\n");
        fclaw_global_essentialf("Observation %d: restarting adjoint%d from checkpoint file "
                                "fort_frame_%04d.checkpoint\n", obs_id, j, nout);
        fclaw_global_essentialf("-----------------------------------------\n");

        reinitialize(adjoint_glob[j], adj_dir, nout);
        map_adjoint_to_representer(adjoint_glob[j], representer_glob[j], &c[j], geo, &u1);
        u1.init_flag = 0;

        /* One representer step at a time.  After each step, alpha is brought
           to tau = T - t: restart the adjoint from the frame just below tau
           and advance it to tau. 
        */
        while (n_curr < nout)
        {
            /* Take one step of representer problem */
            printf(">>>            observation %d, representer %d: t = %12.4f\n",
                   obs_id, j, representer_glob[j]->curr_time);

            /* Change into forward directory to store restart files there*/
            change_dir(rep_dir);
            fclaw_opt->checkpoint = 1;

            int n_new = user_run(representer_glob[j]);

            change_dir("..");
            fclaw_opt->checkpoint = 0;

            n_curr = n_new;
            if (n_new == nout)
            {
                /* user_run wrote the final frame and took no step */
                break;
            }

            /* user_run returned n_new. The representer is now in the output 
               interval: t in (t_n, t_{n+1}] with n = n_new. Reversing time, 
               tau = T - t lies in [tau_k, tau_{k+1}) with k = nout - n_new - 1. 
               So frame k is the last adjoint checkpoint at or before tau. 
               E.g. T = 2, nout = 20: t = 0.437 gives n_new = 4, tau = 1.563
               and k = 15 (tau_15 = 1.5). 
            */
            double tau = tfinal - representer_glob[j]->curr_time;
            int    k   = nout - n_new - 1;

            reinitialize(adjoint_glob[j], adj_dir, k);

            /* fclaw_run loops while n < nout, with each outer step having length 
               dt_outer = tfinal/nout. With nout = k+1 and tfinal = (k+1)*(tau - tau_k) 
               it runs one interval, from tau_k exactly to tau. 
            */
            fclaw_opt_adjoint->nout   = k + 1;
            fclaw_opt_adjoint->tfinal = (k + 1)*(tau - k*dT);
            fclaw_run(adjoint_glob[j]);
            printf(">>>            observation %d, adjoint %d: frame %d -> tau = %12.4f (target %12.4f)\n",
                   obs_id, j, k, adjoint_glob[j]->curr_time, tau);

            /* After the advance, map adjoint to representer */
            map_adjoint_to_representer(adjoint_glob[j], representer_glob[j], &c[j], geo, &u1);
        }
    }

    /* ---------------------- 3. Coefficients beta ------------------------- */
    double *beta = new double[mdata];
    int have_beta = 0;
#ifdef INSTANTANEOUS_VERIFICATION_ONLY
    fclaw_global_essentialf("Instantaneous verification complete; "
                            "skipping the coefficient-system solve.\n");
#else
    /* Observations: dm, W_eps (the same in every adjoint glob) */
    have_beta = compute_beta(model_glob, representer_glob,
                             adjoint_get_options(adjoint_glob[0]), mdata, beta);
    if (!have_beta)
    {
        fclaw_global_essentialf("No gauges (output-gauges = F or none in gauges.data): "
                                "skipping the coefficient solve and q_hat.\n");
    }
#endif

    /* ------------- 4. Optimal estimate q_hat = q_F + sum beta_j r_j --------- */
    if (have_beta)
    {
        const int nout = fclaw_get_options(model_glob)->nout;
        for (int iframe = 0; iframe <= nout; iframe++)
        {
            /* Load checkpoint iframe of the model: prior q_F */
            fclaw_global_essentialf("-----------------------------------------\n");
            fclaw_global_essentialf("Restarting model from checkpoint file "
                                    "fort_frame_%04d.checkpoint\n", iframe);
            fclaw_global_essentialf("-----------------------------------------\n");
            reinitialize(model_glob, "model", iframe);

            for (int j = 0; j < mdata; j++)
            {
                int obs_id = j + 1;
                SET_OBSERVATION_INDEX(&obs_id);

                /* Load checkpoint iframe of the representer */
                char rep_dir[BUFSIZ];
                snprintf(rep_dir, sizeof(rep_dir), "forward%d", j);
                fclaw_global_essentialf("-----------------------------------------\n");
                fclaw_global_essentialf("Observation %d: restarting representer %d from "
                                        "checkpoint file fort_frame_%04d.checkpoint\n",
                                        obs_id, j, iframe);
                fclaw_global_essentialf("-----------------------------------------\n");
                reinitialize(representer_glob[j], rep_dir, iframe);

                map_representer_to_model(representer_glob[j], model_glob, &c[j], geo, j, beta);
            }
            /* Call output frame on updated model (model_q, that is, q_hat)*/
            fclaw_output_frame (model_glob, iframe);
        }
    }

    delete[] beta;

    /* Finalize solvers */
    for (int j = 0; j < mdata; j++)
    {
        adjoint_finalize(adjoint_glob[j]);
        forward_finalize(representer_glob[j]);
    }
    model_finalize(model_glob);
}

/* ----------------------------- Globs -----------------------------------------

   The three kinds of run are separate solvers, each with its own user
   options, domain and initialization:

       model       : model_options_t,   forward_create_domain
       adjoint     : adjoint_options_t, adjoint_create_domain
       representer : forward_options_t, forward_create_domain

   Registering a section returns one struct per section name.  If every
   adjoint (or representer) glob stored that struct, a change to one glob's
   options (e.g. nout, checkpoint) would change them all.  So each glob gets
   its own copy of the options read from the ini file (the *_parsed structs),
   and delete_glob frees the copies after destroying the glob. */

typedef enum
{
    RUN_MODEL,
    RUN_ADJOINT,
    RUN_REPRESENTER
} run_role_t;

/* A glob with its own copies of the parsed options */
static
fclaw_global_t* new_glob(fclaw_app_t *app,
                         fclaw_options_t *fclaw_parsed,
                         fclaw_clawpatch_options_t *clawpatch_parsed,
                         fc2d_clawpack46_options_t *claw46_parsed,
                         void *user_parsed,
                         run_role_t role)
{
    fclaw_options_t *fclaw_opt = new fclaw_options_t;
    memcpy(fclaw_opt, fclaw_parsed, sizeof(fclaw_options_t));

    fclaw_clawpatch_options_t *clawpatch_opt = new fclaw_clawpatch_options_t;
    memcpy(clawpatch_opt, clawpatch_parsed, sizeof(fclaw_clawpatch_options_t));

    fc2d_clawpack46_options_t *claw46_opt = new fc2d_clawpack46_options_t;
    memcpy(claw46_opt, claw46_parsed, sizeof(fc2d_clawpack46_options_t));

    fclaw_global_t *glob = fclaw_global_new(app);
    fclaw_options_store           (glob, fclaw_opt);
    fclaw_clawpatch_options_store (glob, clawpatch_opt);
    fc2d_clawpack46_options_store (glob, claw46_opt);

    switch (role)
    {
    case RUN_MODEL:
    {
        model_options_t *user_opt = new model_options_t;
        memcpy(user_opt, user_parsed, sizeof(model_options_t));
        model_options_store(glob, user_opt);
        forward_create_domain(glob);
        break;
    }
    case RUN_ADJOINT:
    {
        adjoint_options_t *user_opt = new adjoint_options_t;
        memcpy(user_opt, user_parsed, sizeof(adjoint_options_t));
        adjoint_options_store(glob, user_opt);
        adjoint_create_domain(glob);
        break;
    }
    case RUN_REPRESENTER:
    {
        forward_options_t *user_opt = new forward_options_t;
        memcpy(user_opt, user_parsed, sizeof(forward_options_t));
        forward_options_store(glob, user_opt);
        forward_create_domain(glob);
        break;
    }
    }
    return glob;
}

/* Destroy a glob made by new_glob and free its option copies */
static
void delete_glob(fclaw_global_t *glob, run_role_t role)
{
    fclaw_options_t           *fclaw_opt     = fclaw_get_options(glob);
    fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    fc2d_clawpack46_options_t *claw46_opt    = fc2d_clawpack46_get_options(glob);
    const void *user_opt = NULL;
    switch (role)
    {
    case RUN_MODEL:       user_opt = model_get_options(glob);   break;
    case RUN_ADJOINT:     user_opt = adjoint_get_options(glob); break;
    case RUN_REPRESENTER: user_opt = forward_get_options(glob); break;
    }

    fclaw_global_destroy(glob);

    delete fclaw_opt;
    delete clawpatch_opt;
    delete claw46_opt;
    switch (role)
    {
    case RUN_MODEL:       delete (model_options_t*)   user_opt; break;
    case RUN_ADJOINT:     delete (adjoint_options_t*) user_opt; break;
    case RUN_REPRESENTER: delete (forward_options_t*) user_opt; break;
    }
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
    fclaw_options_t *adjoint_fclaw_opt_parsed =
        fclaw_options_register(app, "adjoint", "adjoint_options.ini");
    fclaw_clawpatch_options_t *adjoint_clawpatch_opt_parsed =
        fclaw_clawpatch_2d_options_register(app, "adjoint-clawpatch", "adjoint_options.ini");
    fc2d_clawpack46_options_t *adjoint_claw46_opt_parsed =
        fc2d_clawpack46_options_register(app, "adjoint-clawpack46", "adjoint_options.ini");
    adjoint_options_t *adjoint_user_opt_parsed =
        adjoint_options_register(app, "adjoint_options.ini");

    /* representer options - register once to read from ini file */
    fclaw_options_t           *representer_fclaw_opt_parsed =
        fclaw_options_register(app, "forward", "forward_options.ini");
    fclaw_clawpatch_options_t *representer_clawpatch_opt_parsed =
        fclaw_clawpatch_2d_options_register(app, "forward-clawpatch", "forward_options.ini");
    fc2d_clawpack46_options_t *representer_claw46_opt_parsed =
        fc2d_clawpack46_options_register(app, "forward-clawpack46", "forward_options.ini");
    forward_options_t         *representer_user_opt_parsed =
        forward_options_register(app, "forward_options.ini");

    

    /* model (prior) options - register once to read from ini file */
    fclaw_options_t           *model_fclaw_opt_parsed =
        fclaw_options_register(app, "model", "model_options.ini");
    fclaw_clawpatch_options_t *model_clawpatch_opt_parsed =
        fclaw_clawpatch_2d_options_register(app, "model-clawpatch", "model_options.ini");
    fc2d_clawpack46_options_t *model_claw46_opt_parsed =
        fc2d_clawpack46_options_register(app, "model-clawpack46", "model_options.ini");
    model_options_t           *model_user_opt_parsed =
        model_options_register(app, "model_options.ini");



    /* Read configuration file(s) */
    int first_arg;
    fclaw_exit_type_t vexit = 
        fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");
        
    if (!vexit)
    {
        /* Options have been checked and are valid */
        const int mdata = adjoint_user_opt_parsed->mdata;   /* number of observations */

        /* One adjoint and one representer per observation, and the model */
        fclaw_global_t **adjoint_glob = new fclaw_global_t*[mdata];
        for (int j = 0; j < mdata; j++)
        {
            adjoint_glob[j] = new_glob(app, adjoint_fclaw_opt_parsed,
                                       adjoint_clawpatch_opt_parsed,
                                       adjoint_claw46_opt_parsed,
                                       adjoint_user_opt_parsed, RUN_ADJOINT);
        }

        fclaw_global_t **representer_glob = new fclaw_global_t*[mdata];
        for (int j = 0; j < mdata; j++)
        {
            representer_glob[j] = new_glob(app, representer_fclaw_opt_parsed,
                                           representer_clawpatch_opt_parsed,
                                           representer_claw46_opt_parsed,
                                           representer_user_opt_parsed,
                                           RUN_REPRESENTER);
        }

        fclaw_global_t *model_glob = new_glob(app, model_fclaw_opt_parsed,
                                              model_clawpatch_opt_parsed,
                                              model_claw46_opt_parsed,
                                              model_user_opt_parsed, RUN_MODEL);

        /* initialize all the solvers before doing overlap */
        
        for (int j = 0; j < mdata; j++)
        {
            adjoint_initialize(adjoint_glob[j]);
            forward_initialize(representer_glob[j]); 
        }
        model_initialize(model_glob); 

        /* Set up overlap points in "c" */
        overlap_forward_metadata_t *c = new overlap_forward_metadata_t[mdata];

        for (int j = 0; j < mdata; j++)
        {
            memset(&c[j], 0, sizeof(overlap_forward_metadata_t));
        }

        /* Geometry of the adjoint problem */
        adjoint_geometry_t adjoint_geometry, *adjoint_geo = &adjoint_geometry;        

        /* initialize the adjoint geometry information that is needed for
            * mapping between the representer and the adjoint domain */
        for (int j = 0; j < mdata; j++)
        {
            setup_overlap(representer_glob[j], adjoint_glob[j],&c[j], adjoint_geo);

        }
        /* Run the program */
        run_program(model_glob,
                    representer_glob,
                    adjoint_glob,
                    c, adjoint_geo, mdata);

        for (int j = 0; j < mdata; j++)
        {
            /* destroy overlap points */
            if (c[j].query_points != NULL)
            {
                sc_array_destroy(c[j].query_points);
                c[j].query_points = NULL;
            }
        }
        delete[] c;

        for (int j = 0; j < mdata; j++)
        {
            delete_glob(adjoint_glob[j], RUN_ADJOINT);
            delete_glob(representer_glob[j], RUN_REPRESENTER);
        }
        delete_glob(model_glob, RUN_MODEL);

        delete[] adjoint_glob;
        delete[] representer_glob;
    }

    fclaw_app_destroy (app);

    return 0;
}
