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
 * \file mmg2d/starpu_2d.c
 * \brief Mesh adaptation functions.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Mariem Makni (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 * Starpu codelets and functions.
 *
 */

#ifdef USE_STARPU

#include "starpu_2d.h"

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

/*
 *    Codelet to neutral element initializer
 */
struct starpu_codelet izero_codelet =
{
  .cpu_funcs = {izero_cpu},
  .cpu_funcs_name = {"izero_cpu"},
  .modes = {STARPU_W},
  .nbuffers = 1,
  .name = "izero"
};

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Integer initializtion to 0 for StarPU codelet.
 *
 */
void izero_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);

  *a = 0;
}

/*
 *    Codelet to perform the reduction of two elements using sum operator
 */
struct starpu_codelet accumulate_codelet =
{
  .cpu_funcs = {accumulate_cpu},
  .cpu_funcs_name = {"redux_sum_func"},
  .modes = {STARPU_RW|STARPU_COMMUTE, STARPU_R},
  .nbuffers = 2,
  .name = "redux_sum"
};

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Reduction by sum for StarPU codelet.
 *
 */
void accumulate_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);
  int *b = (int *)STARPU_VARIABLE_GET_PTR(descr[1]);

  *a = *a + *b;
}

/*
 *    Codelet element with value 1 initializer
 */
struct starpu_codelet ione_codelet =
{
  .cpu_funcs = {ione_cpu},
  .cpu_funcs_name = {"ione_cpu"},
  .modes = {STARPU_W},
  .nbuffers = 1,
  .name = "ione"
};

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Integer initializtion to 1 for StarPU codelet.
 *
 */
void ione_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);

  *a = 1;
}

/*
 *    Codelet to perform the reduction of two elements using min operator
 */
struct starpu_codelet min_codelet =
{
  .cpu_funcs = {min_cpu},
  .cpu_funcs_name = {"redux_min_func"},
  .modes = {STARPU_RW|STARPU_COMMUTE, STARPU_R},
  .nbuffers = 2,
  .name = "redux_min"
};

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Reduction by min for StarPU codelet.
 *
 */
void min_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);
  int *b = (int *)STARPU_VARIABLE_GET_PTR(descr[1]);

  *a = MG_MIN ( *a , *b );
}

/*
 *    Codelet to call print_cpu func.
 */
struct starpu_codelet print_codelet =
{
  .cpu_funcs = {print_cpu},
  .cpu_funcs_name = {"print_cpu_func"},
  .modes = {STARPU_R},
  .nbuffers = 1,
  .name = "print"
};

void print_cpu(void *descr[], void *cl_arg)
{
  (void)cl_arg;
  int *a = (int *)STARPU_VARIABLE_GET_PTR(descr[0]);
}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to anaelt function to be called inside StarPU codelet.
 *
 */
void MMG2D_starpu_anaelt(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_ns, *handle_hash ;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;
  int *ns;
  MMG5_Hash *hash;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_hash = (struct starpu_variable_interface *) buffers[2];
  hash = (MMG5_Hash *) STARPU_VARIABLE_GET_PTR(handle_hash);

  handle_ns = (struct starpu_variable_interface *) buffers[3];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *ns += MMG2D_anaelt(mesh,met,hash,typchk,color);
}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to colelt function to be called inside StarPU codelet.
 *
 */
void MMG2D_starpu_colelt(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;

  int *nc;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nc += MMG2D_colelt(mesh,met,typchk,color);
}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to swpmsh function to be called inside StarPU codelet.
 *
 */
void MMG2D_starpu_swpmsh(void *buffers[], void *cl_arg) {

  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nsw;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;

  int *nsw;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_nsw = (struct starpu_variable_interface *) buffers[2];
  nsw = (int *)STARPU_VARIABLE_GET_PTR(handle_nsw);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nsw += MMG2D_swpmsh(mesh,met,typchk,color);

}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to adpspl function to be called inside StarPU codelet.
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
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_ns = (struct starpu_variable_interface *) buffers[2];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg,&color);

  *ns += MMG2D_adpspl(mesh,met,color);

}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to adpcol function to be called inside StarPU codelet.
 *
 */
void MMG2D_starpu_adpcol(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int color;
  int *nc;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg,&color);

  *nc += MMG2D_adpcol(mesh,met,color);
}

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

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to movtri function to be called inside StarPU codelet.
 *
 */
void MMG2D_starpu_movtri(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met;
  struct starpu_vector_interface *vect_mesh, *vect_met;
  struct starpu_variable_interface *handle_nm;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int maxit,improve,color;

  int *nm;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  vect_met = (struct starpu_vector_interface *) buffers[1];
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  assert ( (nx_met = STARPU_VECTOR_GET_NX(vect_met)) == 1 );

  handle_nm = (struct starpu_variable_interface *) buffers[2];
  nm = (int *)STARPU_VARIABLE_GET_PTR(handle_nm);

  starpu_codelet_unpack_args(cl_arg, &maxit, &improve,&color);

  *nm += MMG2D_movtri(mesh,met,maxit,improve,color);

}

/*
 * Codelet to compute dependencies between splitting tasks.
 */
struct starpu_codelet spldep_codelet =
{
  .cpu_funcs = {MMG2D_starpu_spldep},
  .cpu_funcs_name = {"MMG2D_starpu_spldep"},
  .nbuffers = 2,
  .modes = {STARPU_R,STARPU_W},
  .specific_nodes = 2,
  .nodes = {STARPU_SPECIFIC_NODE_CPU,STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "spldep"
};


/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Compute split dependencies: we have dependencies between two colors if they
 * are adjacent throug an edge.
 *
 */
void MMG2D_starpu_spldep(void *buffers[], void *cl_arg) {

  int nx_mesh;
  struct starpu_vector_interface *vect_mesh;
  struct starpu_variable_interface *handle_ier;

  MMG5_pMesh   mesh;
#warning don't build on MSVC
  starpu_tag_t deps[mesh->ncolors];
  int color;
  int k;
  int *ier;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  handle_ier = (struct starpu_variable_interface *) buffers[1];
  ier = (int *)STARPU_VARIABLE_GET_PTR(handle_ier);

  starpu_codelet_unpack_args(cl_arg, &color);

  assert ( mesh->adja );

  memset(deps,0x0,mesh->ncolors*sizeof(starpu_tag_t));

  *ier = 1;
  for ( k=1; k<=mesh->nt; ++k ) {
    MMG5_pTria pt1, pt2;

    pt1 = &mesh->tria[k];

    if ( !MMG2D_EOK(pt1,color) ) continue;

    /* Triangle has to be treated (it has the color treated by the task and is
     * used) */
    int i;
    for ( i=0; i<3; ++i ) {
      /* Search the color of adjacent triangles */
      int k2 = mesh->adja [ 3*(k-1)+1+i ]/3;
      if ( !k2 ) continue;

      pt2 = &mesh->tria[k2];

      if ( pt2->color1 != color ) {
        /* Mark the dependency */
        deps[pt2->color1] = pt2->color1;
      }
    }
  }

 #warning can we pass "null tags" inside the array that will be ignored??? otherwise, deps array has to be packed
  starpu_tag_declare_deps_array((starpu_tag_t)color, mesh->ncolors, &deps[1]);
}


/*
 * Codelet to compute dependencies between move tasks.
 */
struct starpu_codelet movdep_codelet =
{
  .cpu_funcs = {MMG2D_starpu_movdep},
  .cpu_funcs_name = {"MMG2D_starpu_movdep"},
  .nbuffers = 2,
  .modes = {STARPU_R,STARPU_W},
  .specific_nodes = 2,
  .nodes = {STARPU_SPECIFIC_NODE_CPU,STARPU_SPECIFIC_NODE_CPU},
  .where = STARPU_CPU,
  .name = "movdep"
};


/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Compute move dependencies: we have dependencies between two colors if they
 * are connected by an edge (i.e. it exist an edge such as one of the extremity
 * of the edge belongs to one color and the other extremity to the other color).
 *
 */
void MMG2D_starpu_movdep(void *buffers[], void *cl_arg) {

  int nx_mesh;
  struct starpu_vector_interface *vect_mesh;
  struct starpu_variable_interface *handle_ier;

  MMG5_pMesh   mesh;
  starpu_tag_t deps[MMG_NDEPSMAX];
  int color;
  int k;
  int *ier;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  assert ( (nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh)) == 1 );

  handle_ier = (struct starpu_variable_interface *) buffers[1];
  ier = (int *)STARPU_VARIABLE_GET_PTR(handle_ier);

  starpu_codelet_unpack_args(cl_arg, &color);

  assert ( mesh->adja );

  int ndep = 0;
  memset(deps,0x0,MMG_NDEPSMAX*sizeof(starpu_tag_t));

  *ier = 1;
  for ( k=1; k<=mesh->nt; ++k ) {
    MMG5_pTria pt1, pt2;

    pt1 = &mesh->tria[k];

    if ( !MMG2D_EOK(pt1,color) ) continue;

    /* Triangle has to be treated (it has the color treated by the task and is
     * used) */
    int i;
    for ( i=0; i<3; ++i ) {
      /* Search the color of adjacent triangles */
      int k2 = mesh->adja [ 3*(k-1)+1+i ]/3;
      if ( !k2 ) continue;

      pt2 = &mesh->tria[k2];

      if ( pt2->color1 != color ) {
        if ( ndep >= MMG_NDEPSMAX ) {
          fprintf(stderr,"  # Error: %s: %d: number of dependencies exceed the"
                  " maximal authorized value.\n",__func__,__LINE__);

          *ier = 0;
          return;
        }

        /* Mark the dependency */
        deps[pt2->color1] = pt2->color1;
      }
    }
  }

  starpu_tag_declare_deps_array((starpu_tag_t)color, ndep, &deps[1]);
}

/**
 * Compute deps for move operator: for now each color travels the entire mesh to
 * build its deps. A future improvement can be to build in // the primary hash
 * table (\hash), to synchronize threads, then to build the second hash table,
 * to synchronize threads and last to build the color deps.
 */
int MMG2D_movdep(MMG5_pMesh mesh, int color) {
  MMG5_pTria   pt;
  int          k, i;
#warning don't build on MSVC
  starpu_tag_t deps[mesh->ncolors];

  /* Store list of colors per points using a hashtable */
  MMG5_HashP *hash;

  if ( !MMG5_hashNew( mesh,&hash,mesh->np, 3.01*mesh->np) ) return 0;

  for ( k=1; k<=mesh->nt; ++k ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for ( i=0; i<3; ++i ) {
      MMG5_hashPoint(mesh,&hash,pt->v[i],pt->color1);
    }
  }

  /* Append to each point the colors of points that are
   * connected to the current point by an edge. */

  if ( !MMG5_hashNew( mesh,&hash2,mesh->np, 6.01*mesh->np) ) return 0;

  for ( k=1; k<=mesh->nt; ++k ) {

    if ( !MG_EOK(pt) ) continue;

    pt = &mesh->tria[k];
    for ( i=0; i<3; ++i ) {
      int ip1  = pt->v[i];
      int ip2  = pt->v[MMG5_inxt2[i]];

      assert ( hash->item[ip1] && hash->item[ip2] && "no color");

      /* Append colors of point 2 to point 1 */
      MMG5_hpoint  *ph = &hash->item[ip2];
      while ( ph->nxt ) {
        MMG5_hashPoint(mesh,&hash2,ip1,ph->data);
        ph = &hash->item[ph->nxt];
      }

      /* Append colors of point 1 to point 2 */
      MMG5_hpoint  *ph = &hash->item[ip1];
      while ( ph->nxt ) {
        MMG5_hashPoint(mesh,&hash2,ip2,ph->data);
        ph = &hash->item[ph->nxt];
      }
    }
  }

  /* Declare dependencies for color \a color */
  for ( k=1; k<=mesh->nt; ++k ) {

    if ( !MMG2D_EOK(pt) ) continue;

    pt = &mesh->tria[k];
    for ( i=0; i<3; ++i ) {
      ip1 = pt->v[i];
      depx[]
    }
  }

  return 1;
}

#endif
