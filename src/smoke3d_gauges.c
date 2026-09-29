/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#include "smoke3d_gauges.h"

#include "fclaw_gauges.h"

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw_options.h>
#include <fclaw_global.h>

#include "2d/transport2d_fort.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct smoke3d_user
{
    int level;
    double tcurr;
    double qvar;  /* Store concentration */
} smoke3d_user_t;


void smoke3d_read_gauges_data(fclaw_global_t *glob, 
                              fclaw_gauge_t **gauges,
                              int *num_gauges, int *dim)
{
    /* Idea is that we may only need to change this file when updating to newer
       geoclaw code */
/* Sample gauges.data file */
/*
18              
   1   8.6001000000e+01   1.0000000000e-03   0.000000e+00   1.000000e+10
   2   8.6876000000e+01   1.0000000000e-03   0.000000e+00   1.000000e+10
   3   8.7751000000e+01   1.0000000000e-03   0.000000e+00   1.000000e+10
   ....
*/

    FILE *f_gauges_data = fopen("gauges.data","r");

    FCLAW_ASSERT(f_gauges_data != NULL);

    int max_line_len = 200;  /* Maximum line length in file  gauges.data */
    char *line = FCLAW_ALLOC(char,max_line_len);

    /* Skip header/comment/blank lines */
    do {
        fgets(line, max_line_len, f_gauges_data);
    } while (line[0] == '#' || line[0] == '\n');

    /* Read dimension line */
    //fgets(line, max_line_len, f_gauges_data);
    *dim = strtod(line,NULL);

    /* Read number of gauges */
    fgets(line, max_line_len, f_gauges_data);
    *num_gauges = strtod(line,NULL);

    printf("dimension = %d, num_gauges = %d\n", *dim, *num_gauges);

    printf("we are here now\n");

    if (*num_gauges == 0)
    {
        *gauges = NULL;
    }
    else
    {
        printf("are we here yet?\n");
        fclaw_gauges_allocate(glob,*num_gauges,gauges);
        fclaw_gauge_t *g = *gauges;

        int *num = FCLAW_ALLOC(int,   *num_gauges);
        double *xc  = FCLAW_ALLOC(double,*num_gauges);
        double *yc  = FCLAW_ALLOC(double,*num_gauges);
        double *t1  = FCLAW_ALLOC(double,*num_gauges);
        double *t2  = FCLAW_ALLOC(double,*num_gauges);
        double *min_time_increment = FCLAW_ALLOC(double,*num_gauges);


        /* Read gauge info */
        char *next;  
        for(int i = 0; i < *num_gauges; i++)
        {                        
            fgets(line,max_line_len, f_gauges_data);
            num[i] = strtod(line,&next);
            xc[i] = strtod(next,&next);
            yc[i] = strtod(next,&next);
            t1[i] = strtod(next,&next);
            t2[i] = strtod(next,NULL);
        }

        for(int i = 0; i < *num_gauges; i++)
            min_time_increment[i] = 0;

        for(int i = 0; i < *num_gauges; i++)
        {
            double zc = 0;
            fclaw_gauges_set_data(glob,&g[i],num[i],*dim,
                                 xc[i],yc[i],zc,xc[i],yc[i],zc,
                                 t1[i],t2[i],
                                 min_time_increment[i]);
        }
        FCLAW_FREE(num);
        FCLAW_FREE(xc);
        FCLAW_FREE(yc);
        FCLAW_FREE(t1);
        FCLAW_FREE(t2);
        FCLAW_FREE(min_time_increment);

    }   /* End of num_gauges > 0 loop */

    /* Finish up */
    fclose(f_gauges_data);    
    FCLAW_FREE(line);
}

/* This function can be virtualized so the user can specify their 
   gauge output */

void smoke3d_create_gauge_files(fclaw_global_t *glob, 
                                fclaw_gauge_t *gauges,
                                int num_gauges)
{

    /* -----------------------------------------------------
    Open output gauge files and add header information
    ----------------------------------------------------- */
    char filename[15];    /* gaugexxxxx.txt  + EOL */
    FILE *fp;

    //int num_eqns = 1;  /* meqn + 1 (h, hu, hv, eta) */
    for (int i = 0; i < num_gauges; i++)
    {
        int num;
        int dim;
        double xc,yc,zc,t1,t2;
        double x0, y0, z0;
        fclaw_gauges_get_data(glob,&gauges[i],&num, &dim, &x0, &y0, &z0,
                             &xc, &yc, &zc, &t1, &t2);

        sprintf(filename,"gauge%05d.txt",num);
        fp = fopen(filename, "w");
        /* This must have exactly the spacing indicated below. 
           See line 99 of $CLAW/pyclaw/src/pyclaw/gauges.py */
        if (dim == 2)
            fprintf(fp, "# gauge_id= %5d location=( %17.10e %17.10e )\n",
                    num, x0, y0);
        else if (dim == 3)
            fprintf(fp, "# gauge_id= %5d location=( %17.10e %17.10e %17.10e)\n",
                    num, x0, y0,z0);
        fprintf(fp, "# Columns: level x, y, z, time q0   q1   q2 .... aux0  aux1  aux 2...\n");
        fclose(fp);
    }
}

void smoke3d_gauge_normalize_coordinates(fclaw_global_t *glob, 
                                         fclaw_block_t *block,
                                         int blockno, 
                                         fclaw_gauge_t *g,
                                         double *xc, double *yc, double *zc)
{
    /*  
       Map gauge to normalized coordinates in a global [0,1]x[0,1]  domain.

       Gauge coordinates (g->xc,g->yc) are whatever the user supplied above 

       Return normalized (xc,yc) coordinates for gauge.
    */

    MAP_LATLONG2UNITREGION(&g->xc,&g->yc,xc,yc);

}



void smoke3d_gauge_update(fclaw_global_t* glob,
                          fclaw_block_t* block,
                          fclaw_patch_t* patch, 
                          int blockno, 
                          int patchno,
                          double tcurr, 
                          fclaw_gauge_t *g)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int num;
    int dim;
    double t1, t2;
    double x_long, y_lat, zc;
    double x0, y0, z0;
    fclaw_gauges_get_data(glob,g,&num, &dim, &x0, &y0, &z0,
                         &x_long, &y_lat, &zc, &t1, &t2);

    double x_long_low, y_lat_low;
    MAP_BRICK2LATLONG(&blockno,&xlower,&ylower,&x_long_low,&y_lat_low);

    double xupper = xlower + mx*dx;
    double yupper = ylower + my*dy;
    double x_long_hi, y_lat_hi;
    MAP_BRICK2LATLONG(&blockno,&xupper,&yupper,&x_long_hi,&y_lat_hi);

    if (!(x_long >= x_long_low && x_long <= x_long_hi))
        return;

    if (!(y_lat >= y_lat_low && y_lat <= y_lat_hi))
        return;


#if 1   
    double xc,yc;
    MAP_LATLONG2UNITREGION(&x_long,&y_lat,&xc,&yc);

#endif        


    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);


    /* Interpolate q variables and aux variables (bathy only for now)
       to gauge location */
    double qvar;
    SMOKE3D_UPDATE_GAUGE(&blockno, &mx,&my,&mbc,&meqn,&xlower,&ylower,
                         &dx,&dy,q,&xc,&yc,&qvar);
                
    /* Store qvar, avar in gauge buffers;  Anything stored will be printed
       in print_buffers */
    smoke3d_user_t *guser = FCLAW_ALLOC(smoke3d_user_t,1);

    guser->level = patch->level;
    guser->tcurr = tcurr;

    for(int m = 0; m < meqn; m++)
        guser->qvar = qvar;

    fclaw_gauges_set_buffer_entry(glob,g,guser);
}


void smoke3d_print_gauges(fclaw_global_t *glob, 
                          fclaw_gauge_t *gauge) 
{
    /* This assumes on buffers be organized as an array; entries
       start at 0 and with kmax-1 */
    smoke3d_user_t **gauge_buffer;
    int kmax;
    fclaw_gauges_get_buffer(glob,gauge,&kmax,(void***) &gauge_buffer);

    int id = fclaw_gauges_get_id(glob,gauge);

    char filename[15];  /* gaugexxxxx.txt + EOL character */
    sprintf(filename,"gauge%05d.txt",id);

    FILE *fp = fopen(filename, "a");
    for(int k = 0; k < kmax; k++)
    {
        smoke3d_user_t *guser = gauge_buffer[k];

        double q = guser->qvar;
        q = fabs(q) < 1e-99 ? 0 : q; /* For reading in Matlab */
        fprintf(fp, "%5d %15.7e %15.7e\n",
                guser->level, guser->tcurr,guser->qvar);

        FCLAW_FREE(guser);
    }
    fclose(fp);
}

#ifdef __cplusplus
}
#endif

