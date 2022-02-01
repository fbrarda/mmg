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
#include <starpu.h>
#include "metis_mmg.h"
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include "libmmg2d.h"
#include <stdbool.h>
#include "mmgcommon.h"

struct starpu_codelet colelt_codelet =
{
  .cpu_funcs = {MMG2D_starpu_colelt},
  .cpu_funcs_name = {"MMG2D_starpu_colelt"},
  .nbuffers = 3,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "colelt"
};


struct starpu_codelet swpmsh_codelet =
{
  .cpu_funcs = {MMG2D_starpu_swpmsh},
  .cpu_funcs_name = {"MMG2D_starpu_swpmsh"},
  .nbuffers = 3,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "swpmsh"
};


struct starpu_codelet anaelt_codelet =
{
  .cpu_funcs = {MMG2D_starpu_anaelt},
  .cpu_funcs_name = {"MMG2D_starpu_anaelt"},
  .nbuffers = 4,
  .modes = {STARPU_RW, STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "anaelt"
};

struct starpu_codelet movtri_codelet =
{
  .cpu_funcs = {MMG2D_starpu_movtri},
  .cpu_funcs_name = {"MMG2D_starpu_movtri"},
  .nbuffers = 3,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "movtri"
};

struct starpu_codelet adpspl_codelet =
{
  .cpu_funcs = {MMG2D_starpu_adpspl},
  .cpu_funcs_name = {"MMG2D_starpu_adpspl"},
  .nbuffers = 3,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "adpspl"
};

struct starpu_codelet adpcol_codelet =
{
  .cpu_funcs = {MMG2D_starpu_adpcol},
  .cpu_funcs_name = {"MMG2D_starpu_adpcol"},
  .nbuffers = 3,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "adpcol"
};

struct starpu_codelet hashTria_codelet =
{
  .cpu_funcs = {MMG2D_starpu_hashTria},
  .cpu_funcs_name = {"MMG2D_starpu_hashTria"},
  .nbuffers = 1,
  .modes = {STARPU_RW},
  .specific_nodes = 1,
  .nodes = {STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "hashTria"
};

void izero_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);

  *a = 0;
}

/*
 *    Codelet to neutral element initializer
 */
static struct starpu_codelet izero_codelet =
{
  .cpu_funcs = {izero_cpu},
  .cpu_funcs_name = {"izero_cpu"},
  .modes = {STARPU_W},
  .nbuffers = 1,
  .name = "izero"
};


void accumulate_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);
  int *b = (int *)STARPU_VARIABLE_GET_PTR(descr[1]);

  *a = *a + *b;
}

/*
 *    Codelet to perform the reduction of two elements
 */
static struct starpu_codelet accumulate_codelet =
{
  .cpu_funcs = {accumulate_cpu},
  .cpu_funcs_name = {"redux_cpu_func"},
  .modes = {STARPU_RW|STARPU_COMMUTE, STARPU_R},
  .nbuffers = 2,
  .name = "redux"
};


void print_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);

  printf("print *a %d\n",*a);
}
/*
 *    Codelet to perform the reduction of two elements
 */
static struct starpu_codelet print_codelet =
{
  .cpu_funcs = {print_cpu},
  .cpu_funcs_name = {"print_cpu_func"},
  .modes = {STARPU_R},
  .nbuffers = 1,
  .name = "redux"
};

/***********************************************************/

static starpu_pthread_mutex_t mutex;
pthread_mutex_t count_mutex;
pthread_mutex_t count_mutex1;

/**************************************************************************************************/

/*identifier les triangles qui sont adjacents à un triangle de la tâche, et les inclure dans les boucles des step 2,3,4 */
//static inline
int MMG2D_ADJEOK( MMG5_pMesh mesh, int k, int color1 ) {
  int *adja;
  int isAdja;
  int i;

  /** Step 1: mesh adjacency creation: adja allocation */
  /*if ( (!mesh->adja) && (1 != MMG2D_hashTria(mesh) ) ) {
	  fprintf(stderr,"\n  ## Error: %s: unable to create "
	  "adjacency table.\n",__func__);
	  return 0;
    }*/


  adja = &mesh->adja[3*(k-1) + 1];
  //pt   = &mesh->tria[k];
  isAdja = 0;
  for( i = 0; i < 3; i++ ) {
    if( MMG2D_EOK( &mesh->tria[(adja[i]/3)], color1 ) ) {
      isAdja = 1;
      break;
    }
  }

  return isAdja;
}


int MMG2D_overlap(MMG5_pMesh mesh, idx_t* part)
{

  fprintf(stdout,"  -- Overlapping  \n");

  int *adja;

  int status, i,k;
  int NbAdj;
  int *nadjcy;

  MMG5_pTria pt;
  MMG5_pTria pt1;

  for( k = 1; k <= mesh->nt; k++ ) {
    adja = &mesh->adja[3*(k-1) + 1];
    pt   = &mesh->tria[k];

    //Test if MG_EOK is valid
    if ( MG_EOK(pt) )
    {
      for( i = 0; i < 3; i++ ) {
        pt1=&mesh->tria[(adja[i]/3)];

        if( (pt1->color2) > (pt->color1)) {

          //fprintf(stdout,"----Elt = %d, color_elt = %d , adja = %d, color_adj= %d \n", k, pt ->color1, (adja[i]/3), pt1->color2 );

          pt1->color2 = pt->color1;
          pt->ref = pt->color1;
          pt1->ref = pt1->color2;
        }

      }

    }

  }

  return 1;

}

/* Mesh adaptation routine for the first stages of the algorithm: intertwine splitting
   based on patterns, collapses and swaps.
   typchk = 1 -> adaptation based on edge lengths
   typchk = 2 -> adaptation based on lengths calculated in metric met */
int MMG2D_anatri(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  MMG5_Hash hash;
  int       it,maxit,ns,nc,nsw,nns,nnc,nnsw;
  int       ret;
  int       color, i;
  starpu_data_handle_t vector_mesh, vector_met, vector_hash;
  starpu_data_handle_t handle_ns, handle_nc, handle_nsw;

  //fprintf(stdout," Mesh Computation: Hello anatri function ----------\n");

  starpu_vector_data_register(&vector_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, 1, sizeof(MMG5_pMesh));
  starpu_vector_data_register(&vector_met, STARPU_MAIN_RAM, (uintptr_t)met, 1, sizeof(MMG5_pSol));

  starpu_variable_data_register(&vector_hash, STARPU_MAIN_RAM, (uintptr_t)&hash, sizeof(MMG5_Hash));

  starpu_variable_data_register(&handle_ns, STARPU_MAIN_RAM, (uintptr_t)&ns, sizeof(ns));
  starpu_variable_data_register(&handle_nc, STARPU_MAIN_RAM, (uintptr_t)&nc, sizeof(nc));
  starpu_variable_data_register(&handle_nsw, STARPU_MAIN_RAM, (uintptr_t)&nsw, sizeof(nsw));

  starpu_data_set_reduction_methods ( handle_ns, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods ( handle_nc, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods ( handle_nsw, &accumulate_codelet, &izero_codelet) ;

  /* Main routine; intertwine split, collapse and swaps */
  nns = nnc = nnsw = 0;
  it = 0;
  maxit = 5;

  starpu_data_acquire(handle_ns, STARPU_W);
  starpu_data_acquire(handle_nc, STARPU_W);
  starpu_data_acquire(handle_nsw, STARPU_W);

  do {
    if ( typchk == 2 && it == 0 )  mesh->info.fem = 1;

    if ( !mesh->info.noinsert ) {
      /* Memory free */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;

      /* Split long edges according to patterns */
      //ns = MMG2D_anaelt(mesh,met,typchk, color);

      if ( !MMG5_hashNew(mesh,&hash,mesh->np,3*mesh->np) ) return 0;
      //starpu_pthread_mutex_init(&mutex, NULL);

      ns = 0;
      starpu_data_release(handle_ns);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&anaelt_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_RW, vector_hash,
                                 STARPU_REDUX, handle_ns,
                                 STARPU_VALUE, &typchk, sizeof(typchk),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
        //ns = MMG2D_anaelt(mesh,met,typchk, color);
      }
      starpu_data_acquire(handle_ns, STARPU_RW);

      //starpu_pthread_mutex_destroy(&mutex);

      if ( ns < 0 ) {
        fprintf(stderr,"  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }

      MMG5_DEL_MEM(mesh,hash.item);

      /* Recreate adjacencies */
      if ( !MMG2D_hashTria(mesh) ) {
        fprintf(stdout,"  ## Hashing problem. Exit program.\n");
        return 0;
      }

      /* fprintf(stdout," Begin insert hashtria codelet \n"); */

      /* ret = starpu_task_insert(&hashTria_codelet, */
      /*                STARPU_RW, vector_mesh,  */
      /*                            0); */

      /* STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");  */

      /* fprintf(stdout," End insert hashTria codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

      /* starpu_task_wait_for_all(); */

      /* Collapse short edges */
      //nc = MMG2D_colelt(mesh,met,typchk,color);

      //fprintf(stdout," Begin insert colelt codelet \n");

      nc = 0;
      starpu_data_release(handle_nc);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&colelt_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_nc,
                                 STARPU_VALUE, &typchk, sizeof(typchk),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
      }

      /* fprintf(stdout," End insert colelt codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

      starpu_data_acquire(handle_nc,STARPU_RW);
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
      // nsw = MMG2D_swpmsh(mesh,met,typchk,color);

      // fprintf(stdout," Begin insert swpmsh codelet \n");

      nsw = 0;
      starpu_data_release(handle_nsw);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&swpmsh_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_nsw,
                                 STARPU_VALUE, &typchk, sizeof(typchk),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

      }
      /* fprintf(stdout," End insert swpmsh codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

      starpu_data_acquire(handle_nsw,STARPU_RW);
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

  starpu_data_release(handle_ns);
  starpu_data_release(handle_nc);
  starpu_data_release(handle_nsw);

  starpu_data_unregister(vector_mesh);
  starpu_data_unregister(vector_met);
  starpu_data_unregister(vector_hash);
  starpu_data_unregister(handle_ns);
  starpu_data_unregister(handle_nc);
  starpu_data_unregister(handle_nsw);

  return 1;
}


void MMG2D_starpu_anaelt(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_ns, *handle_hash ;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int8_t typchk;
  int color;
  int *ns;
  MMG5_Hash *hash;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);

  handle_hash = (struct starpu_variable_interface *) buffers[2];
  hash = (MMG5_Hash *) STARPU_VARIABLE_GET_PTR(handle_hash);

  handle_ns = (struct starpu_variable_interface *) buffers[3];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *ns = MMG2D_anaelt(mesh,met,hash,typchk,color);
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

  /*workerid is an integer between 0 and starpu_worker_get_count() - 1 */
  /*Return the identifier of the current worker, i.e the one associated to the calling thread*/
  int worker_id = starpu_worker_get_id();
  enum starpu_worker_archtype type = starpu_worker_get_type(worker_id);
  char workername[128];

  //starpu_worker_get_name(worker_id, workername, 128);

  //fprintf(stdout," --Call anaelt codelet =, number of CPU workers= %d, worker_id= %d, Worker_name %s: , hash_address= %p \n", starpu_worker_get_count(), starpu_worker_get_id(), workername, hash );

  /*******************************************/


  /*******************************************/

  //if ( !MMG5_hashNew(mesh,hash,mesh->np,3*mesh->np) ) return 0;

  //int   k;
  //int hsiz=mesh->np; int hmax=3*mesh->np;

  /* adjust hash table params */
  /*hash->siz  = hsiz+1;
    hash->max  = hmax + 2;
    hash->nxt  = hash->siz;

    pthread_mutex_lock(&count_mutex);

    MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hedge),"hash table",
    return 0);

    fprintf(stdout,"hash.C,, worker_id= %d \n",starpu_worker_get_id() );

    MMG5_SAFE_CALLOC(hash->item,(hash->max+1),MMG5_hedge,return 0);

    pthread_mutex_unlock(&count_mutex);

    for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;*/

  /* Step 1: travel mesh, check edges, and tag those to be split; create the new vertices in hash */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    //fprintf(stdout,"----Elt = %d, color_elt = %d  \n", k, pt ->color1 );

    if ( !MMG2D_EOK(pt,color1) || (pt->ref < 0) ) continue;
    if ( MG_SIN(pt->tag[0]) || MG_SIN(pt->tag[1]) || MG_SIN(pt->tag[2]) )  continue;

    /* Check if pt should be cut */
    pt->flag = 0;

    /* typchk=1 : split on geometric basis and rough size considerations */
    if ( typchk == 1) {
      if ( !MMG2D_chkedg(mesh,k) ) continue;
    }

    /* typchk =2 : split on basis of edge lengths in the metric */
    else if ( typchk ==2 ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( len > MMG2D_LLONG )
        {MG_SET(pt->flag,i);


        }
      }
    }

    /* mesh->info.fem : split edges which are not MG_BDY, but whose vertices are both MG_BDY */
    if ( mesh->info.fem ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]]; //point
        p2 = &mesh->point[pt->v[i2]];

        if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) && !(pt->tag[i] & MG_BDY) ) {
          MG_SET(pt->flag,i);
          fprintf(stdout," worker_id= %d, Worker_name %s: , hash_address= %p \n",starpu_worker_get_id(), workername, hash );
        }
      }
    }

    if ( !pt->flag ) continue;
    ns++;

    /* Travel edges to split and create new points */
    for (i=0; i<3; i++) {
      if ( !MG_GET(pt->flag,i) ) continue;

      i1  = MMG5_inxt2[i]; //i1 sommet du triangle k
      i2  = MMG5_iprv2[i];//i2 sommet du triangle k
      ip1 = pt->v[i1]; //v[3] est le tableau de sommet du triangle k
      ip2 = pt->v[i2];
      ip = MMG5_hashGet(hash,ip1,ip2); //parcourir les arretes à decouper

      if ( ip > 0 ) continue;

      /* Geometric attributes of the new point */
      ier = MMG2D_bezierCurv(mesh,k,i,s,o,no);

      if ( !ier ) {
        MG_CLR(pt->flag,i);
        continue;
      }

      //pthread_mutex_lock(&count_mutex);
      ip = MMG2D_newPt(mesh,o,pt->tag[i]);  //tag contains binary flag

      //pthread_mutex_unlock(&count_mutex);

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


      // TODO: ici on met lock;
      /* Add point to the hashing structure */
      //starpu_worker_lock(worker_id);
      MMG5_hashEdge(mesh,hash,ip1,ip2,ip);
      //starpu_worker_unlock(worker_id);
    }
  }
  if ( !ns ) {
    MMG5_DEL_MEM(mesh,hash->item);
    return ns;
  }

  /* Step 2: Make flags at triangles consistent between themselves (check if adjacent triangle is split) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    //if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;
    if ( !MMG2D_EOK(pt,color1) || pt->ref < 0)
      continue;
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
      //if ( !MMG2D_EOK(pt,color1) || pt->ref < 0 ) continue;MG_EOK(pt)
      if ( !MMG2D_EOK(pt,color1) || pt->ref < 0)
        continue;
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
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0  )
    fprintf(stdout,"     worker %d: %8d splitted\n", starpu_worker_get_id(),ns);

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
    /* compute new position */
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

  /* restore coords of last valid pos. */
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


void MMG2D_starpu_colelt(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int8_t typchk;
  int color;

  int *nc;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nc = MMG2D_colelt(mesh,met,typchk,color);
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

  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     worker %d: %8d vertices removed\n", starpu_worker_get_id(),nc);

  return nc;
}

/* Travel triangles and swap edges to improve quality */
void MMG2D_starpu_swpmsh(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nsw;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int8_t typchk;
  int color;

  int *nsw;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);

  handle_nsw = (struct starpu_variable_interface *) buffers[2];
  nsw = (int *)STARPU_VARIABLE_GET_PTR(handle_nsw);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nsw = MMG2D_swpmsh(mesh,met,typchk,color);

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
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"      worker %d: %8d edge swapped\n",starpu_worker_get_id(),nns);

  return nns;
}


/* Mesh adaptation routine for the final stage of the algorithm: intertwine splitting
   based on patterns, collapses, swaps and vertex relocations.*/
int MMG2D_adptri(MMG5_pMesh mesh,MMG5_pSol met) {

  int it,nns,ns,nnc,nc,nnsw,nsw,nnm,nm;
  int typchk;
  int color;
  int ret, i;
  int maxit;
  int8_t improve;

  //fprintf(stdout," Mesh Adaptation: Hello adptri function-----------\n");

  starpu_data_handle_t vector_mesh,vector_met, handle_ns, handle_nc,handle_nsw,handle_nm;

  nns = nnc = nnsw = nnm = it = 0;

  starpu_vector_data_register(&vector_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, 1, sizeof(MMG5_pMesh));
  starpu_vector_data_register(&vector_met, STARPU_MAIN_RAM, (uintptr_t)met, 1, sizeof(MMG5_pSol));

  starpu_variable_data_register(&handle_ns, STARPU_MAIN_RAM, (uintptr_t)&ns, sizeof(ns));
  starpu_variable_data_register(&handle_nc, STARPU_MAIN_RAM, (uintptr_t)&nc, sizeof(nc));
  starpu_variable_data_register(&handle_nsw, STARPU_MAIN_RAM, (uintptr_t)&nsw, sizeof(nsw));
  starpu_variable_data_register(&handle_nm, STARPU_MAIN_RAM, (uintptr_t)&nm, sizeof(nsw));

  starpu_data_set_reduction_methods ( handle_ns, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods ( handle_nc, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods ( handle_nsw, &accumulate_codelet, &izero_codelet) ;
  starpu_data_set_reduction_methods ( handle_nm, &accumulate_codelet, &izero_codelet) ;

  starpu_data_acquire(handle_ns,  STARPU_W);
  starpu_data_acquire(handle_nc,  STARPU_W);
  starpu_data_acquire(handle_nsw, STARPU_W);
  starpu_data_acquire(handle_nm,  STARPU_W);

  do {
    if ( !mesh->info.noinsert ) {
      //ns = MMG2D_adpspl(mesh,met,color);

      //fprintf(stdout," Begin insert adpspl codelet \n");

      ns = 0;
      starpu_data_release(handle_ns);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&adpspl_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_ns,
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

      }

      starpu_data_acquire(handle_ns, STARPU_RW);

      /* fprintf(stdout," End insert adpspl codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

      if ( ns < 0 ) {
        fprintf(stderr,"  ## Problem in function adpspl."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }

      //fprintf(stdout," Begin insert adpcol codelet \n");

      //nc = MMG2D_adpcol(mesh,met,color);

      nc = 0;
      starpu_data_release(handle_nc);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&adpcol_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_nc,
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

      }
      starpu_data_acquire(handle_nc,STARPU_RW);

      /* fprintf(stdout," End insert adpcol codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

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
      //nsw = MMG2D_swpmsh(mesh,met,2,color);

      //fprintf(stdout," Begin insert swpmsh codelet \n");

      typchk=2;

      nsw = 0;
      starpu_data_release(handle_nsw);

      for (color=0; color< mesh->ncolors; color++)
      {

        ret = starpu_task_insert(&swpmsh_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_nsw,
                                 STARPU_VALUE, &typchk, sizeof(typchk),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

      }
      starpu_data_acquire(handle_nsw,STARPU_RW);

      /* fprintf(stdout," End insert swpmsh1 codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Problem in function swpmsh."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nsw = 0;

    if ( !mesh->info.nomove ) {
      //nm = MMG2D_movtri(mesh,met,1,0,color);

      //fprintf(stdout," Begin insert movtri1 codelet \n");

      maxit=1;
      improve=0;

      nm = 0;
      starpu_data_release(handle_nm);

      for (color=0; color< mesh->ncolors; color++)
      {
        ret = starpu_task_insert(&movtri_codelet,
                                 STARPU_RW, vector_mesh,
                                 STARPU_RW, vector_met,
                                 STARPU_REDUX, handle_nm,
                                 STARPU_VALUE, &maxit, sizeof(maxit),
                                 STARPU_VALUE, &improve, sizeof(improve),
                                 STARPU_VALUE, &color, sizeof(color),
                                 0);

        STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

      }

      starpu_data_acquire(handle_nm,STARPU_RW);

      /* fprintf(stdout," End insert movtri1 codelet \n"); */
      /* fprintf(stdout," -------------------------\n"); */
      /* fprintf(stdout," -------------------------\n"); */

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
    maxit=5;
    improve=1;
    //nm = MMG2D_movtri(mesh,met,5,1,color);

    //fprintf(stdout," Begin insert movtri2 codelet \n");

    nm = 0;
    starpu_data_release(handle_nm);

    for (color=0; color< mesh->ncolors; color++)
    {
      ret = starpu_task_insert(&movtri_codelet,
                               STARPU_RW, vector_mesh,
                               STARPU_RW, vector_met,
                               STARPU_REDUX, handle_nm,
                               STARPU_VALUE, &maxit, sizeof(maxit),
                               STARPU_VALUE, &improve, sizeof(improve),
                               STARPU_VALUE, &color, sizeof(color),
                               0);

      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

    }
    starpu_data_acquire(handle_nm,STARPU_RW);

    /* fprintf(stdout," End insert movtri2 codelet \n"); */
    /* fprintf(stdout," -------------------------\n"); */
    /* fprintf(stdout," -------------------------\n"); */

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

  starpu_data_release(handle_ns);
  starpu_data_release(handle_nc);
  starpu_data_release(handle_nsw);
  starpu_data_release(handle_nm);

  starpu_data_unregister(vector_mesh);
  starpu_data_unregister(vector_met);
  starpu_data_unregister(handle_ns);
  starpu_data_unregister(handle_nc);
  starpu_data_unregister(handle_nsw);
  starpu_data_unregister(handle_nm);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return -1 if failed or number of new points.
 *
 * Analysis and splitting routine for edges in the final step of the algorithm;
 * edges are only splitted on a one-by-one basis
 *
 */

void MMG2D_starpu_adpspl(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_ns;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int color;
  int *ns;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);


  handle_ns = (struct starpu_variable_interface *) buffers[2];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg,&color);

  *ns = MMG2D_adpspl(mesh,met,color);

}

int MMG2D_adpspl(MMG5_pMesh mesh,MMG5_pSol met, int color1) {
  MMG5_pTria         pt;
  double             lmax,len;
  int                k,ns,nt,ip,ier;
  int8_t             i,i1,i2,imax;

  ns = 0;

  //fprintf(stdout," ----------Call adpspl codelet -----------\n");

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

void MMG2D_starpu_adpcol(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int color;
  int *nc;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg,&color);

  *nc = MMG2D_adpspl(mesh,met,color);

}

/* Analysis and collapse routine for edges in the final step of the algorithm */
int MMG2D_adpcol(MMG5_pMesh mesh,MMG5_pSol met, int color1) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            len;
  int               k,nc,ilist,list[MMG2D_LONMAX+2];
  int8_t            i,i1,i2,open;

  //fprintf(stdout," ----------Call adpcol codelet -----------\n");
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

/* Analyze points to relocate them according to a quality criterion */
void MMG2D_starpu_movtri(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nm;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int maxit,improve,color;

  int *nm;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);

  handle_nm = (struct starpu_variable_interface *) buffers[2];
  nm = (int *)STARPU_VARIABLE_GET_PTR(handle_nm);

  starpu_codelet_unpack_args(cl_arg, &maxit, &improve,&color);

  *nm = MMG2D_movtri(mesh,met,maxit,improve,color);

}

int MMG2D_movtri(MMG5_pMesh mesh,MMG5_pSol met,int maxit,int8_t improve, int color1) {
  MMG5_pTria           pt;
  MMG5_pPoint          p0;
  int                  base,k,nnm,nm,ns,it,ilist,list[MMG2D_LONMAX+2];
  int8_t               i,ier;

  it = nnm = 0;
  base = 0;

  //fprintf(stdout," ----------Call movtri codelet -----------\n");

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
    if ( mesh->info.ddebug )
      fprintf(stdout,"     worker %d: %8d moved, %d geometry\n",starpu_worker_get_id(),nm,ns);
  }
  while ( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm > 0 )
    fprintf(stdout,"     worker %d: %8d vertices moved, %d iter.\n",starpu_worker_get_id(),nnm,it);

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

  /* Stage 0: mesh coloration with metis*/

  //fprintf(stdout,"  -- Call Metis for partioning  \n");

  idx_t *part;

  //int color1;

  /** Allocate the table part */
  MMG5_SAFE_CALLOC(part,mesh->nt,idx_t,return 0);
  int status, i;
  int status1;

  MMG5_pTria pt;

  //fprintf(stdout,"  --Begin Call Metis---- \n");

  status=MMG_part_meshElts2metis( mesh, part, (idx_t)mesh->ncolors );

  for (i=0; i< mesh->nt; i++)
  {

	  pt= &mesh->tria[i+1];
    pt->color1 = part[i];
    pt->color2 = pt->color1;
    pt->ref = pt->color1;

    //fprintf(stdout,"----nbelts= %d, --color1= %d \n",i, pt->color1);
  }

  //fprintf(stdout,"  --END Call Metis------ \n");

  /*save result*/
  if ( MMG2D_saveMesh(mesh, "metis.mesh") != 1 )
    exit(EXIT_FAILURE);


  //status1=MMG2D_overlap(mesh, part);

  /*save result*/
  /* if ( MMG2D_saveMesh(mesh, "metise.mesh") != 1 )
     exit(EXIT_FAILURE);*/

  /*************************************/

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

  if ( !MMG2D_anatri(mesh,met,2)  ) {
    fprintf(stderr,"  ## Unable to proceed adaptation. Exit program.\n");
    return 0;
  }

  /* Stage 3: fine mesh improvements */
  if ( !MMG2D_adptri(mesh,met) ) {
    fprintf(stderr,"  ## Unable to make fine improvements. Exit program.\n");
    return 0;
  }




  return 1;
}
