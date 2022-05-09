/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/
/**
 * \file mmg2d/mmg2d1.c
 * \brief Mesh adaptation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */


#include "mmg2d.h"

#ifdef USE_STARPU
#include <pthread.h>
#include "partitionning_2d.h"
#include "starpu_2d.h"

/* Mutex used for lock/unlock hash access */
pthread_mutex_t lock;

#endif

/* Mesh adaptation routine for the first stages of the algorithm: intertwine splitting
 based on patterns, collapses and swaps.
   typchk = 1 -> adaptation based on edge lengths
   typchk = 2 -> adaptation based on lengths calculated in metric met */
int MMG2D_anatri(MMG5_pMesh mesh,MMG5_pSol met,int typchk) {
  MMG5_HashP hash2;
  MMG5_Hash  hash;
  int        it,maxit,ns,nc,nsw,nns,nnc,nnsw;

  /* Main routine; intertwine split, collapse and swaps */
  nns = nnc = nnsw = 0;
  it = 0;
  maxit = 5;

#ifdef USE_STARPU
  int       color,ret,ier;
  starpu_data_handle_t handle_mesh,handle_met,handle_hash,handle_deps;
  starpu_data_handle_t handle_ns,handle_nc,handle_nsw,handle_ier;

  starpu_variable_data_register(&handle_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, sizeof(MMG5_pMesh));
  starpu_variable_data_register(&handle_met , STARPU_MAIN_RAM, (uintptr_t)met , sizeof(MMG5_pSol));

  starpu_variable_data_register(&handle_hash, STARPU_MAIN_RAM, (uintptr_t)&hash, sizeof(MMG5_Hash));

  starpu_variable_data_register(&handle_ns  , STARPU_MAIN_RAM, (uintptr_t)&ns  , sizeof(ns));
  starpu_variable_data_register(&handle_nc  , STARPU_MAIN_RAM, (uintptr_t)&nc  , sizeof(nc));
  starpu_variable_data_register(&handle_nsw , STARPU_MAIN_RAM, (uintptr_t)&nsw , sizeof(nsw));
  starpu_variable_data_register(&handle_ier , STARPU_MAIN_RAM, (uintptr_t)&ier , sizeof(ier));

  starpu_data_set_reduction_methods(handle_ns , &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_nc , &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_nsw, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_ier, &min_codelet       , &ione_codelet) ;

  starpu_data_acquire(handle_ns , STARPU_W);
  starpu_data_acquire(handle_nc , STARPU_W);
  starpu_data_acquire(handle_nsw, STARPU_W);

  /*  Creation of a useless array to have one data per task. This data is
      only used as a tool to express the task dependencie (alloc is done here as
      I think that we will use a variable number of colors in the future). */
  char *colors;
  MMG5_SAFE_MALLOC(colors,mesh->info.ncolors+1,char,return 0);

  /* Creation of an array of handle toward the colors array. */
  starpu_data_handle_t *handle_per_colors;
  MMG5_SAFE_MALLOC(handle_per_colors,mesh->info.ncolors+1,starpu_data_handle_t,return 0);

  /* Make handle point toward color data */
  for (color=0; color<=mesh->info.ncolors; color++) {
    starpu_variable_data_register(&handle_per_colors[color], STARPU_MAIN_RAM,
                                  (uintptr_t)&colors[color], sizeof(char));
  }
#endif

  do {
    if ( typchk == 2 && it == 0 ) {
// #warning Luca: check consistency with 3D
    }

#warning add mesh coloration here (each waves of remeshing unbalance the color repartition)?

    if ( !mesh->info.noinsert ) {

#ifndef USE_STARPU
      /* Memory free */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;
#endif

      /* Split long edges according to patterns */
      if ( !MMG5_hashNew(mesh,&hash,mesh->np,3*mesh->np) ) return 0;

#ifdef USE_STARPU
      ns = 0;
      starpu_data_release(handle_ns);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        /* Call anaelt-like function for partition \a color. Can be wrapped into
         * a task to be run in // (in this case, we can gain some alloc/unalloc
         * using a STARPU_SCRATCH handle and letting Starpu deal with the local
         * array \a deps ) */
        int ier = MMG2D_starpu_anaelt(mesh,&handle_mesh,&handle_met,handle_per_colors,
                                      &handle_hash,&handle_ns,typchk,color);

        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit anaelt task to starPU. Exit program.\n");
          return 0;
        }
      }

      /* Remark: adja array is not needed for anaelt_codelet but is needed to
       * compute dependencies inside the \a MMG2D_starpu_anaelt. To avoid
       * concurrency at adja deallocation, it is deallocated outside the color
       * loop. */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;

      starpu_data_acquire(handle_ns, STARPU_RW);

#else
      ns = MMG2D_anaelt(mesh,met,&hash,typchk,MMG_NOCOLOR);
#endif

      if ( ns < 0 ) {
        fprintf(stderr,"  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }
      /* Memory free */
      MMG5_DEL_MEM(mesh,hash.item);

      /* Recreate adjacencies */
      if ( !MMG2D_hashTria(mesh) ) {
        fprintf(stdout,"  ## Hashing problem. Exit program.\n");
        return 0;
      }

#ifdef USE_STARPU
      /** Even if splitting operator don't change deps, it is needed to
       * recompute them as now, 2 domains connected by an edge can't be remeshed
       * at the same time (which was allowed for splits). Note that a possible
       * improvement can be to directly compute the collapse dependencies for
       * the splitting operator (it reduces the split parallelization but splits
       * are fast) and it save 1 deps computation. */

      /** For each point, compute the list of colors that can be reach directly
       * or by a 1 edge connection. */
      if ( !MMG2D_pointColor_1edg(mesh,&hash2) ) {
        fprintf(stderr,"  ## Problem in second step of dependencies construction"
                " for moving operator."
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nc = 0;
      starpu_data_release(handle_nc);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        ret = starpu_task_insert(&colelt_codelet,
                                 STARPU_RW, handle_mesh,
                                 STARPU_RW, handle_met,
                                 STARPU_REDUX, handle_nc,
                                 STARPU_VALUE, &typchk, sizeof(typchk),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert:colelt_codelet");
      }

      MMG5_DEL_MEM(mesh,hash2.item);
      starpu_data_acquire(handle_nc,STARPU_RW);
#else
      nc = MMG2D_colelt(mesh,met,typchk,MMG_NOCOLOR);
#endif

      if ( nc < 0 ) {
        fprintf(stderr,"  ## Unable to collapse mesh. Exiting.\n");
        return 0;
      }

    }
    else {
      ns = 0;
      nc = 0;
    }

    /* Swap edges */
    if ( !mesh->info.noswap ) {

#ifdef USE_STARPU
      /** Collapse operators modifty the connection between colors so it is
       * mandatory to recompute deps here. */
      /** For each point, compute the list of colors that can be reach directly
       * or by a 1 edge connection. */
      if ( !MMG2D_pointColor_1edg(mesh,&hash2) ) {
        fprintf(stderr,"  ## Problem in second step of dependencies construction"
                " for moving operator."
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nsw = 0;
      starpu_data_release(handle_nsw);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        /* Call swpmsh-like function that compute dependencies (from points hash
         * table) and submit movtri task for partition \a color. */
        int ier = MMG2D_starpu_swpmsh(mesh,&hash2,&handle_mesh,&handle_met,
                                      handle_per_colors,
                                      &handle_nsw,typchk,color);
        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit swpmsh task to starPU. Exit program.\n");
          return 0;
        }
      }
      MMG5_DEL_MEM(mesh,hash2.item);

      starpu_data_acquire(handle_nsw,STARPU_RW);
#else
      nsw = MMG2D_swpmsh(mesh,met,typchk,MMG_NOCOLOR);
#endif
      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else nsw = 0;

    nns += ns;
    nnc += nc;
    nnsw += nsw;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nsw);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nsw >0 );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnsw,it);
  }

#ifdef USE_STARPU

  for (color=0; color<=mesh->info.ncolors; color++) {
   // printf("color %d\n",color);
    starpu_data_unregister(handle_per_colors[color]);
  }

  MMG5_SAFE_FREE(colors);
  MMG5_SAFE_FREE(handle_per_colors);

  starpu_data_release(handle_ns);
  starpu_data_release(handle_nc);
  starpu_data_release(handle_nsw);

  starpu_data_unregister(handle_mesh);
  starpu_data_unregister(handle_met);
  starpu_data_unregister(handle_hash);
  starpu_data_unregister(handle_ns);
  starpu_data_unregister(handle_nc);
  starpu_data_unregister(handle_nsw);
#endif

  return 1;
}


/* Travel triangles and split long edges according to patterns */
int MMG2D_anaelt(MMG5_pMesh mesh,MMG5_pSol met, MMG5_Hash *hash, int typchk,int color1) {
  MMG5_pTria      pt;
  MMG5_pPoint     ppt,p1,p2;
  double          len,s,o[2],no[2];
  int             ns,nc,npinit,ni,k,nt,ip1,ip2,ip,it,vx[3];
  int8_t          i,ic,i1,i2,ier;
  static int8_t   mmgWarn0=0;

  s = 0.5;
  ns = 0;
  npinit = mesh->np;

  /* Step 1: travel mesh, check edges, and tag those to be split; create the new vertices in hash */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color1) || (pt->ref < 0) ) continue;
    if ( MG_SIN(pt->tag[0]) || MG_SIN(pt->tag[1]) || MG_SIN(pt->tag[2]) )  continue;

    /* Check if pt should be cut */
    pt->flag = 0;

    /* typchk=1 : split on geometric basis and rough size considerations */
    if ( typchk == 1) {
      if ( !MMG2D_chkedg(mesh,k) ) continue;
    }
    /* typchk=2 : split on basis of edge lengths in the metric */
    else if ( typchk ==2 ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( len > MMG2D_LLONG ) MG_SET(pt->flag,i);
      }
    }

    /* mesh->info.fem : split edges which are not MG_BDY, but whose vertices are both MG_BDY */
    if ( mesh->info.fem ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];
        if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) && !(pt->tag[i] & MG_BDY) ) MG_SET(pt->flag,i);
      }
    }

    if ( !pt->flag ) continue;
    ns++;

    /* Travel edges to split and create new points */
    for (i=0; i<3; i++) {
      if ( !MG_GET(pt->flag,i) ) continue;
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip = MMG5_hashGet(hash,ip1,ip2);
      if ( ip > 0 ) continue;

      /* Geometric attributes of the new point */
      ier = MMG2D_bezierCurv(mesh,k,i,s,o,no);
      if ( !ier ) {
        MG_CLR(pt->flag,i);
        continue;
      }
      ip = MMG2D_newPt(mesh,o,pt->tag[i]);
      if ( !ip ) {
        /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            do {
                              MMG2D_delPt(mesh,mesh->np);
                            } while ( mesh->np>npinit );return -1;,
                            o,pt->tag[i]);
      }
      ppt = &mesh->point[ip];
      if ( MG_EDG(pt->tag[i]) ) {
        ppt->n[0] = no[0];
        ppt->n[1] = no[1];
      }

      /* If there is a metric in the mesh, interpolate it at the new point */
      assert ( met );
      if ( met->m )
        MMG2D_intmet(mesh,met,k,i,ip,s);

#ifdef USE_STARPU
      /* Add a mutex lock for hash table access*/
      pthread_mutex_lock(&lock);
#endif

      MMG5_hashEdge(mesh,hash,ip1,ip2,ip);

#ifdef USE_STARPU
      pthread_mutex_unlock(&lock);
#endif
    }
  }
  if ( !ns ) {
    return ns;
  }

  /* Step 2: Make flags at triangles consistent between themselves (check if adjacent triangle is split) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;
    else if ( pt->flag == 7 ) continue;
    nc = 0;

    for (i=0; i<3; i++) {
      i1 = MMG5_iprv2[i];
      i2 = MMG5_inxt2[i];
      if ( !MG_GET(pt->flag,i) && !MG_SIN(pt->tag[i]) ) {
        ip = MMG5_hashGet(hash,pt->v[i1],pt->v[i2]);
        if ( ip > 0 ) {
          MG_SET(pt->flag,i);
          nc++;
        }
      }
    }
    if ( nc > 0 ) ns++;
  }
  if ( mesh->info.ddebug && ns ) {
    fprintf(stdout,"     %d analyzed  %d proposed\n",mesh->nt,ns);
    fflush(stdout);
  }

  /* Step 3: Simulate splitting and delete points leading to an invalid configuration */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  nc = 0;
  do {
    ni = 0;
    for ( k=1; k<= mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;
      else if ( pt->flag == 0 ) continue;

      vx[0] = vx[1] =vx[2] = 0;
      pt->flag = 0;
      ic = 0;

      for (i=0; i<3; i++) {
        i1 = MMG5_iprv2[i];
        i2 = MMG5_inxt2[i];
        vx[i] = MMG5_hashGet(hash,pt->v[i1],pt->v[i2]);
        if ( vx[i] > 0 ) {
          MG_SET(pt->flag,i);
          if ( mesh->point[vx[i]].flag > 2 )  ic = 1;
        }
      }

      if ( !pt->flag )  continue;
      switch (pt->flag) {
        case 1: case 2: case 4:
          ier = MMG2D_split1_sim(mesh,met,k,vx);
          break;
        case 7:
          ier = MMG2D_split3_sim(mesh,met,k,vx);
          break;
        default:
          ier = MMG2D_split2_sim(mesh,met,k,vx);
          break;
      }
      if ( ier )  continue;

      /* An edge is invalidated in the process */
      ni++;
      if ( ic == 0 && MMG2D_dichoto(mesh,met,k,vx) ) {
        for (i=0; i<3; i++)
          if ( vx[i] > 0 )  mesh->point[vx[i]].flag++;
      }
      /* Relocate point at the center of the edge */
      else {
        for (i=0; i<3; i++) {
          if ( vx[i] > 0 ) {
            p1 = &mesh->point[pt->v[MMG5_iprv2[i]]];
            p2 = &mesh->point[pt->v[MMG5_inxt2[i]]];
            ppt = &mesh->point[vx[i]];
            ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
            ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
          }
        }
      }
    }
    nc += ni;
  }
  while ( ni >0 && ++it <20 );

  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d corrected,  %d invalid\n",nc,ni);
    fflush(stdout);
  }

  /* step 4: effective splitting */
  ns = 0;
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];

    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 )  continue;
    else if ( pt->flag == 0 )  continue;

    vx[0] = vx[1] = vx[2] = 0;
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( MG_GET(pt->flag,i) ) {
        vx[i] = MMG5_hashGet(hash,pt->v[i1],pt->v[i2]);
        if ( !vx[i] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: unable to create point on"
                    " at least 1 edge.\n Exit program.\n",__func__);
          }
          return -1;
        }
      }
    }
    if ( pt->flag == 1 || pt->flag == 2 || pt->flag == 4 ) {
      ier = MMG2D_split1(mesh,met,k,vx);
      ns++;
    }
    else if ( pt->flag == 7 ) {
      ier = MMG2D_split3(mesh,met,k,vx);
      ns++;
    }
    else {
      ier = MMG2D_split2(mesh,met,k,vx);
      ns++;
    }
    if ( !ier ) return -1;
  }


  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 ) {
#ifdef USE_STARPU
    fprintf(stdout,"     worker %d: %8d splitted\n", starpu_worker_get_id(),ns);
#else
    fprintf(stdout,"     %8d splitted\n", ns);
#endif
  }

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param vx pointer toward table of edges to split.
 * \return 1.
 *
 * Find acceptable position for splitting.
 *
 */
int MMG2D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx) {
  MMG5_pTria   pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[3][2],p[3][2];
  float        to,tp,t;
  int          ia,ib,ier,it,maxit;
  int8_t       i,i1,i2;

  pt = &mesh->tria[k];

  /* Get point on curve and along segment for edge split */
  for (i=0; i<3; i++) {
    memset(p[i],0,2*sizeof(double));
    memset(o[i],0,2*sizeof(double));
    if ( vx[i] > 0 ) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      ia = pt->v[i1];
      ib = pt->v[i2];
      pa = &mesh->point[ia];
      pb = &mesh->point[ib];
      ps = &mesh->point[vx[i]];
      o[i][0] = 0.5 * (pa->c[0] + pb->c[0]);
      o[i][1] = 0.5 * (pa->c[1] + pb->c[1]);
      p[i][0] = ps->c[0];
      p[i][1] = ps->c[1];
    }
  }
  maxit = 4;
  it = 0;
  tp = 1.0;
  to = 0.0;

  do {
    /* Compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
    switch (pt->flag) {
      case 1: case 2: case 4:
        ier = MMG2D_split1_sim(mesh,met,k,vx);
        break;
      case 7:
        ier = MMG2D_split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = MMG2D_split2_sim(mesh,met,k,vx);
        break;
    }
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );

  /* Restore coords of last valid pos. */
  if ( !ier ) {
    t = to;
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
  }
  return 1;
}

/* Travel triangles and collapse short edges */
int MMG2D_colelt(MMG5_pMesh mesh,MMG5_pSol met,int typchk,int color1) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,ll,hmin2;
  int          list[MMG2D_LONMAX+2],ilist,nc,k;
  uint8_t      i,i1,i2,open;

  nc = 0;
  hmin2 = mesh->info.hmin * mesh->info.hmin;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;

    /* Travel 3 edges of the triangle and decide whether to collapse p1->p2, based on length criterion */
    pt->flag = 0; // was here before, but probably serves for nothing

    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;

      /* Impossible to collapse a surface point onto a non surface point -- impossible to
       collapse a surface point along a non geometric edge */
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      /* Same test for REF points */
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      open = (mesh->adja[3*(k-1)+1+i] == 0) ? 1 : 0;

      /* Check length */
      if ( typchk == 1 ) {
        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        ll = ux*ux + uy*uy;
        if ( ll > hmin2 ) continue;
      }
      else {
        ll = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( ll > MMG2D_LSHRT ) continue;
      }

      /* Check whether the geometry is preserved */
      ilist = MMG2D_chkcol(mesh,met,k,i,list,typchk);

      if ( ilist > 3 || ( ilist == 3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) ) {
#ifdef USE_STARPU
    fprintf(stdout,"     worker %d: %8d vertices removed\n", starpu_worker_get_id(),nc);
#else
    fprintf(stdout,"     %8d vertices removed\n", nc);
#endif
  }

  return nc;
}

int MMG2D_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,int typchk, int color1) {
  MMG5_pTria pt;
  int        it,maxit,ns,nns,k;
  uint8_t    i;

  it = nns = 0;
  maxit = 2;
  mesh->base++;

  do {
    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;

      for (i=0; i<3; i++) {
        if ( MG_SIN(pt->tag[i]) || MG_EDG(pt->tag[i]) ) continue;
        else if ( MMG2D_chkswp(mesh,met,k,i,typchk) ) {
          ns += MMG2D_swapar(mesh,k,i);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ns > 0 && ++it<maxit );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 ) {
#ifdef USE_STARPU
    fprintf(stdout,"      worker %d: %8d edge swapped\n",starpu_worker_get_id(),nns);
#else
    fprintf(stdout,"      %8d edge swapped\n",nns);
#endif
  }

  return nns;
}


/* Mesh adaptation routine for the final stage of the algorithm: intertwine splitting
 based on patterns, collapses, swaps and vertex relocations.*/
int MMG2D_adptri(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_HashP           hash2;
  int                  maxit,it,nns,ns,nnc,nc,nnsw,nsw,nnm,nm;
  int                  typchk;
  int                  maxit_mov;
  int8_t               improve;

  nns = nnc = nnsw = nnm = it = 0;
  maxit = 5;

#ifdef USE_STARPU
  int                  color;
  int                  ret,i;
  starpu_data_handle_t handle_mesh,handle_met,handle_ns,handle_nc,handle_nsw,handle_nm;

  starpu_variable_data_register(&handle_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, sizeof(MMG5_pMesh));
  starpu_variable_data_register(&handle_met , STARPU_MAIN_RAM, (uintptr_t)met , sizeof(MMG5_pSol));

  starpu_variable_data_register(&handle_ns  , STARPU_MAIN_RAM, (uintptr_t)&ns , sizeof(ns));
  starpu_variable_data_register(&handle_nc  , STARPU_MAIN_RAM, (uintptr_t)&nc , sizeof(nc));
  starpu_variable_data_register(&handle_nsw , STARPU_MAIN_RAM, (uintptr_t)&nsw, sizeof(nsw));
  starpu_variable_data_register(&handle_nm  , STARPU_MAIN_RAM, (uintptr_t)&nm , sizeof(nsw));

  starpu_data_set_reduction_methods(handle_ns , &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_nc , &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_nsw, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods(handle_nm , &accumulate_codelet, &izero_codelet) ;

  starpu_data_acquire(handle_ns,  STARPU_W);
  starpu_data_acquire(handle_nc,  STARPU_W);
  starpu_data_acquire(handle_nsw, STARPU_W);
  starpu_data_acquire(handle_nm,  STARPU_W);

  /*  Creation of a useless array to have one data per task. This data is
      only used as a tool to express the task dependencie (alloc is done here as
      I think that we will use a variable number of colors in the future). */
  char *colors;
  MMG5_SAFE_MALLOC(colors,mesh->info.ncolors+1,char,return 0);

  /* Creation of an array of handle toward the colors array. */
  starpu_data_handle_t *handle_per_colors;
  MMG5_SAFE_MALLOC(handle_per_colors,mesh->info.ncolors+1,starpu_data_handle_t,return 0);

  /* Make handle point toward color data */
  for (color=0; color<=mesh->info.ncolors; color++) {
    starpu_variable_data_register(&handle_per_colors[color], STARPU_MAIN_RAM,
                                  (uintptr_t)&colors[color], sizeof(char));
  }
#endif

  do {
    if ( !mesh->info.noinsert ) {

#ifdef USE_STARPU
      ns = 0;
      starpu_data_release(handle_ns);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        /* Call adpspl-like function for partition \a color. Can be wrapped into
         * a task to be run in // (in this case, we can gain some
         * alloc/unalloc using a STARPU_SCRATCH handle and letting Starpu deal
         * with the local array \a deps ) */
        int ier = MMG2D_starpu_adpspl(mesh,&handle_mesh,&handle_met,handle_per_colors,
                                      &handle_ns,color);

        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit adpspl task to starPU. Exit program.\n");
          return 0;
        }
      }

      starpu_data_acquire(handle_ns, STARPU_RW);
#else
      ns = MMG2D_adpspl(mesh,met,MMG_NOCOLOR);
#endif

      if ( ns < 0 ) {
        fprintf(stderr,"  ## Problem in function adpspl."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }

#ifdef USE_STARPU
      /** Even if splitting operator don't change deps, it is needed to
       * recompute them as now, 2 domains connected by an edge can't be remeshed
       * at the same time (which was allowed for splits). Note that a possible
       * improvement can be to directly compute the collapse dependencies for
       * the splitting operator (it reduces the split parallelization but splits
       * are fast) and it save 1 deps computation. */

      /** For each point, compute the list of colors that can be reach directly
       * or by a 1 edge connection. */
      if ( !MMG2D_pointColor_1edg(mesh,&hash2) ) {
        fprintf(stderr,"  ## Problem in second step of dependencies construction"
                " for moving operator."
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nc = 0;
      starpu_data_release(handle_nc);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        /* Call adpcol-like function that compute dependencies (from points hash
         * table) and submit movtri task for partition \a color. */
        int ier = MMG2D_starpu_adpcol(mesh,&hash2,&handle_mesh,&handle_met,
                                      handle_per_colors,&handle_nc,color);

        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit adpcol task to starPU. Exit program.\n");
          return 0;
        }
      }
      MMG5_DEL_MEM(mesh,hash2.item);
      starpu_data_acquire(handle_nc,STARPU_RW);
#else
      nc = MMG2D_adpcol(mesh,met,MMG_NOCOLOR);
#endif

      if ( nc < 0 ) {
        fprintf(stderr,"  ## Problem in function adpcol."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else {
      ns = 0;
      nc = 0;
    }

    if ( !mesh->info.noswap ) {
      typchk=2;

#ifdef USE_STARPU
      /** Collapse operators modifty the connection between colors so it is
       * mandatory to recompute deps here. */
      /** For each point, compute the list of colors that can be reach directly
       * or by a 1 edge connection. */
      if ( !MMG2D_pointColor_1edg(mesh,&hash2) ) {
        fprintf(stderr,"  ## Problem in second step of dependencies construction"
                " for moving operator."
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nsw = 0;
      starpu_data_release(handle_nsw);

      for (color=1; color<=mesh->info.ncolors; color++)
      {

        /* Call swpmsh-like function that compute dependencies (from points hash
         * table) and submit movtri task for partition \a color. */
        int ier = MMG2D_starpu_swpmsh(mesh,&hash2,&handle_mesh,&handle_met,
                                      handle_per_colors,
                                      &handle_nsw,typchk,color);
        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit swpmsh task to starPU. Exit program.\n");
          return 0;
        }
      }
      MMG5_DEL_MEM(mesh,hash2.item);
      starpu_data_acquire(handle_nsw,STARPU_RW);
#else
      nsw = MMG2D_swpmsh(mesh,met,typchk,MMG_NOCOLOR);
#endif

      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Problem in function swpmsh."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nsw = 0;

    if ( !mesh->info.nomove ) {
      maxit_mov = 1;
      improve = 0;

#ifdef USE_STARPU

      /** For each point, compute the list of colors that can be reach directly
       * or by a 1 edge connection. */
      if ( !MMG2D_pointColor_1edg(mesh,&hash2) ) {
        fprintf(stderr,"  ## Problem in second step of dependencies construction"
                " for moving operator."
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nm = 0;
      starpu_data_release(handle_nm);

      for (color=1; color<=mesh->info.ncolors; color++)
      {
        /* Call movtri-like function that compute dependencies (from points hash
         * table) and submit movtri task for partition \a color. */
        int ier = MMG2D_starpu_movtri(mesh,&hash2,&handle_mesh,&handle_met,handle_per_colors,
                                      &handle_nm,maxit_mov,improve,color);

        if ( ier < 1 ) {
          fprintf(stderr,"  ## Unable to submit movtri task to starPU. Exit program.\n");
          return 0;
        }
      }

      starpu_data_acquire(handle_nm,STARPU_RW);
#else
      nm = MMG2D_movtri(mesh,met,maxit_mov,improve,MMG_NOCOLOR);
#endif

      if ( nm < 0 ) {
        fprintf(stderr,"  ## Problem in function movtri. "
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nm = 0;

    nns  += ns;
    nnc  += nc;
    nnsw += nsw;
    nnm  += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nsw+nm > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nsw,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && (nc+ns+nsw+nm > 0) );

  /* Last iterations of vertex relocation only */
  if ( !mesh->info.nomove ) {
    maxit_mov = 5;
    improve = 1;

#ifdef USE_STARPU
    nm = 0;
    starpu_data_release(handle_nm);

    for (color=1; color<=mesh->info.ncolors; color++)
    {
      /* Call movtri-like function that compute dependencies and submit movtri
       * task for partition \a color. Possible improvement: in fact we don't
       * need to recompute deps as those ones previously computed are still
       * valid so we can try to get it from previous call of starpu_movtri. */
      int ier = MMG2D_starpu_movtri(mesh,&hash2,&handle_mesh,&handle_met,handle_per_colors,
                                    &handle_nm,maxit_mov,improve,color);

      if ( ier < 1 ) {
        fprintf(stderr,"  ## Unable to submit movtri task to starPU. Exit program.\n");
        return 0;
      }
    }
    starpu_data_acquire(handle_nm,STARPU_RW);
    MMG5_DEL_MEM(mesh,hash2.item);

#else
    nm = MMG2D_movtri(mesh,met,maxit_mov,improve,MMG_NOCOLOR);
#endif

    if ( nm < 0 ) {
      fprintf(stderr,"  ## Problem in function movtri. Unable to complete mesh."
              " Exit program.\n");
      return 0;
    }
    nnm += nm;
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",nns,nnc,nnsw,nnm,it);
  }

#ifdef USE_STARPU

  for (color=0; color<=mesh->info.ncolors; color++) {
    starpu_data_unregister(handle_per_colors[color]);
  }

  MMG5_SAFE_FREE(colors);
  MMG5_SAFE_FREE(handle_per_colors);

  starpu_data_release(handle_ns);
  starpu_data_release(handle_nc);
  starpu_data_release(handle_nsw);
  starpu_data_release(handle_nm);

  starpu_data_unregister(handle_mesh);
  starpu_data_unregister(handle_met);
  starpu_data_unregister(handle_ns);
  starpu_data_unregister(handle_nc);
  starpu_data_unregister(handle_nsw);
  starpu_data_unregister(handle_nm);
#endif

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param color1 color to treat (useless without starPU).
 *
 * \return -1 if failed or number of new points.
 *
 * Analysis and splitting routine for edges in the final step of the algorithm;
 * edges are only splitted on a one-by-one basis
 *
 */
int MMG2D_adpspl(MMG5_pMesh mesh,MMG5_pSol met, int color1) {
  MMG5_pTria         pt;
  double             lmax,len;
  int                k,ns,nt,ip,ier;
  int8_t             i,i1,i2,imax;

  ns = 0;

  /*loop until nt to avoid the split of new triangle*/
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;

    imax = -1;
    lmax = -1.0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }

    if ( lmax < MMG2D_LOPTL ) continue;
    else if ( MG_SIN(pt->tag[imax]) ) continue;

    /* Check the feasibility of splitting */
    ip = MMG2D_chkspl(mesh,met,k,imax);

    /* Lack of memory; abort the routine */
    if ( ip < 0 ){
      return ns;
    }
    else if ( ip > 0 ) {

      ier = MMG2D_split1b(mesh,k,imax,ip);

      /* Lack of memory; abort the routine */
      if ( !ier ) {
        MMG2D_delPt(mesh,ip);
        return ns;
      }
      ns += ier;
    }
  }

  return ns;
}

/* Analysis and collapse routine for edges in the final step of the algorithm */
int MMG2D_adpcol(MMG5_pMesh mesh,MMG5_pSol met, int color1) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            len;
  int               k,nc,ilist,list[MMG2D_LONMAX+2];
  int8_t            i,i1,i2,open;

  nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;

    /* Check edge length, and attempt collapse */
    pt->flag = 0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;

      open = ( mesh->adja[3*(k-1)+1+i] == 0 ) ? 1 : 0;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > MMG2D_LOPTS ) continue;

      ilist = MMG2D_chkcol(mesh,met,k,i,list,2);

      if ( ilist > 3 || ( ilist==3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  return nc;
}


int MMG2D_movtri(MMG5_pMesh mesh,MMG5_pSol met,int maxit,int8_t improve, int color1) {
  MMG5_pTria           pt;
  MMG5_pPoint          p0;
  int                  base,k,nnm,nm,ns,it,ilist,list[MMG2D_LONMAX+2];
  int8_t               i,ier;

  it = nnm = 0;
  base = 0;

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];

      if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 )  continue;

      for (i=0; i<3; i++) {
        p0 = &mesh->point[pt->v[i]];
        if ( p0->flag == base || MG_SIN(p0->tag) || p0->tag & MG_NOM ) continue;

        ilist = MMG2D_boulet(mesh,k,i,list);

        if ( MG_EDG(p0->tag) ) {
          ier = MMG2D_movedgpt(mesh,met,ilist,list,improve);
          if ( ier ) ns++;
        }
        else {
          if ( met->size == 3 && met->m )
            ier = MMG2D_movintpt_ani(mesh,met,ilist,list,improve);
          else
            ier = MMG2D_movintpt(mesh,met,ilist,list,improve);
        }

        if ( ier ) {
          nm++;
          p0->flag = base;
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug ) {
#ifdef USE_STARPU
      fprintf(stdout,"     worker %d: %8d moved, %d geometry\n",starpu_worker_get_id(),nm,ns);
#else
      fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
#endif
    }
  }
  while ( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm > 0 ) {
#ifdef USE_STARPU
    fprintf(stdout,"     worker %d: %8d vertices moved, %d iter.\n",starpu_worker_get_id(),nnm,it);
#else
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);
#endif
  }

  return nnm;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 1 if success, 0 if strongly fail.
 *
 * Mesh adaptation -- new version of mmg2d1.c
 *
 **/
int MMG2D_mmg2d1n(MMG5_pMesh mesh,MMG5_pSol met) {

#ifdef USE_STARPU
  pthread_mutex_init(&mesh->lock,NULL);

  /* Stage 0: mesh coloration with metis*/
  int status;
  int ret;

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** MESH PARTITIONNING\n");

  status = MMG2D_part_meshElts(mesh);

  /* StarPU configuration: set the sceduling policy */
  struct starpu_conf conf;
  starpu_conf_init(&conf);
  conf.sched_policy_name = "eager";

  /* StarPU initialization method */
  ret = starpu_init(&conf);
  if (ret == -ENODEV) return 0;
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

  /* STARPU task profiling info */
  starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
#endif

  /* Stage 1: creation of a geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !MMG2D_anatri(mesh,met,1) ) {
    fprintf(stderr,"  ## Unable to split mesh-> Exiting.\n");
    return 0;
  }

  /* Stage 2: creation of a computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  if ( !MMG2D_defsiz(mesh,met) ) {
    fprintf(stderr,"  ## Metric undefined. Exit program.\n");
    return 0;
  }

  MMG5_gradation_info(mesh);
  if ( mesh->info.hgrad > 0. ) {
    if (!MMG2D_gradsiz(mesh,met) ) {
      fprintf(stderr,"  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }

  if ( mesh->info.hgradreq > 0. ) {
    MMG2D_gradsizreq(mesh,met);
  }

  if ( !MMG2D_anatri(mesh,met,2) ) {
    fprintf(stderr,"  ## Unable to proceed adaptation. Exit program.\n");
    return 0;
  }

  /* Stage 3: fine mesh improvements */
  if ( !MMG2D_adptri(mesh,met) ) {
    fprintf(stderr,"  ## Unable to make fine improvements. Exit program.\n");
    return 0;
  }

#ifdef USE_STARPU
  pthread_mutex_destroy(&mesh->lock);

  /* Terminate StarPU */
  starpu_shutdown();
#endif

  return 1;
}
