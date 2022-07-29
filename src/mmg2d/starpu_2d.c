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
  .cpu_funcs = {MMG2D_hashTria_task},
  .cpu_funcs_name = {"MMG2D_hashTria_task"},
  .nbuffers = 1,
  .modes = {STARPU_RW},
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
  .modes = {STARPU_RW, STARPU_R},
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
  .modes = {STARPU_RW, STARPU_R},
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
  .cpu_funcs = {MMG2D_anaelt_task},
  .cpu_funcs_name = {"MMG2D_anaelt_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_RW, STARPU_RW, STARPU_RW,STARPU_REDUX},//,STARPU_DATA_MODE_ARRAY},
  .where = STARPU_CPU,
  .name = "anaelt"
};

/**
 * \param mesh pointer toward the mesh
 * \return 1 if success, 0 if fail
 *
 * Create adjacency relations between the triangles dein the mesh
 *
 * \warning unused: insertion inside the hash table can't be done in parallel using tasks
 */

void MMG2D_hashTria_task(void *buffers[], void *cl_arg) {

  int nx_mesh;
  struct starpu_vector_interface *vect_mesh;

  MMG5_pMesh mesh;

  int ret;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);

  ret=MMG2D_hashTria(mesh);

  if (ret<0)
  {
    fprintf(stdout,"  ## Hashing problem. Exit program.\n");

  }
}

/**
 * \param buffers Codelet buffers (to unpack)
 * \param cl_arg Codelet arguments (to unpack)
 *
 * Wrapper to anaelt function to be called inside StarPU codelet.
 *
 */
void MMG2D_anaelt_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_ns, *handle_hash;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;
  int *ns;
  MMG5_Hash *hash;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_hash = (struct starpu_variable_interface *) buffers[2];
  hash = (MMG5_Hash *)STARPU_VARIABLE_GET_PTR(handle_hash);

  handle_ns = (struct starpu_variable_interface *) buffers[3];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *ns += MMG2D_anaelt(mesh,met,hash,typchk,color);
}

struct starpu_codelet colelt_codelet =
{
  .cpu_funcs = {MMG2D_colelt_task},
  .cpu_funcs_name = {"MMG2D_colelt_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX,STARPU_DATA_MODE_ARRAY},
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
void MMG2D_colelt_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;

  int *nc;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nc += MMG2D_colelt(mesh,met,typchk,color);
}

struct starpu_codelet swpmsh_codelet =
{
  .cpu_funcs = {MMG2D_swpmsh_task},
  .cpu_funcs_name = {"MMG2D_swpmsh_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX,STARPU_DATA_MODE_ARRAY},
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
void MMG2D_swpmsh_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_nsw;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int typchk;
  int color;

  int *nsw;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_nsw = (struct starpu_variable_interface *) buffers[2];
  nsw = (int *)STARPU_VARIABLE_GET_PTR(handle_nsw);

  starpu_codelet_unpack_args(cl_arg, &typchk, &color);

  *nsw += MMG2D_swpmsh(mesh,met,typchk,color);
}

struct starpu_codelet adpspl_codelet =
{
  .cpu_funcs = {MMG2D_adpspl_task},
  .cpu_funcs_name = {"MMG2D_adpspl_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX,STARPU_DATA_MODE_ARRAY},
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
void MMG2D_adpspl_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_ns;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int color;
  int *ns;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_ns = (struct starpu_variable_interface *) buffers[2];
  ns = (int *)STARPU_VARIABLE_GET_PTR(handle_ns);

  starpu_codelet_unpack_args(cl_arg,&color);

  *ns += MMG2D_adpspl(mesh,met,color);
}

struct starpu_codelet adpcol_codelet =
{
  .cpu_funcs = {MMG2D_adpcol_task},
  .cpu_funcs_name = {"MMG2D_adpcol_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_RW, STARPU_RW, STARPU_REDUX,STARPU_DATA_MODE_ARRAY},
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
void MMG2D_adpcol_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_nc;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int color;
  int *nc;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_nc = (struct starpu_variable_interface *) buffers[2];
  nc = (int *)STARPU_VARIABLE_GET_PTR(handle_nc);

  starpu_codelet_unpack_args(cl_arg,&color);

  *nc += MMG2D_adpcol(mesh,met,color);
}

struct starpu_codelet movtri_codelet =
{
  .cpu_funcs = {MMG2D_movtri_task},
  .cpu_funcs_name = {"MMG2D_movtri_task"},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .modes = {STARPU_R, STARPU_R, STARPU_REDUX,STARPU_DATA_MODE_ARRAY},
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
void MMG2D_movtri_task(void *buffers[], void *cl_arg) {

  struct starpu_variable_interface  *handle_mesh, *handle_met;
  struct starpu_variable_interface *handle_nm;

  MMG5_pMesh mesh;
  MMG5_pSol met;

  int maxit,improve,color;

  int *nm;

  handle_mesh = (struct starpu_variable_interface  *) buffers[0];
  mesh = (MMG5_pMesh)STARPU_VARIABLE_GET_PTR(handle_mesh);

  handle_met = (struct starpu_variable_interface  *) buffers[1];
  met = ( MMG5_pSol)STARPU_VARIABLE_GET_PTR(handle_met);

  handle_nm = (struct starpu_variable_interface *) buffers[2];
  nm = (int *)STARPU_VARIABLE_GET_PTR(handle_nm);

  starpu_codelet_unpack_args(cl_arg, &maxit, &improve,&color);

  *nm += MMG2D_movtri(mesh,met,maxit,improve,color);

}

/**
 * \param mesh pointer toward the mesh structure
 * \param deps array to store colors dependencies (has to be of size ncolors+1)
 * \param color color on which we compute deps
 *
 * \return \a ndeps, the number of dependencies for current color.
 *
 * Compute the list of colors with which the color \a color has dependency for
 * split operator: there is a dependency between two colors if they are adjacent.
 *
 */
int MMG2D_spldeps ( MMG5_pMesh mesh,int *deps,int color) {
  int k;

  /** Step 1: mark colors adjacent to current color */
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

  /** Step 2: pack the array of dependencies and count number of dependencies */
  k     = 0;
  int ndeps = mesh->info.ncolors;

  /* Search for the last mark color */
  while ( !deps[ndeps] && ndeps > 0 ) --ndeps;

  while ( ++k < ndeps ) {
    if ( !deps[k] ) {
      /* copy dependency stored in last position in the empty place */
      deps[k]     = deps[ndeps];
      deps[ndeps] = 0;

      /* Search the last marked color */
      while ( !deps[ndeps] ) --ndeps;
    }
  }

  // debug
  /* printf(" Color %d -- %d deps:\n",color,ndeps); */
  /* for ( i=1; i<=ndeps; ++i ) { */
  /*   printf(" %d",deps[i]); */
  /* } */
  /* printf("\n"); */

  return ndeps;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_per_color array to have 1 starpu handle per color
 * \param handle_hash starpu handle toward the hash table of edges
 * \param handle_ns starpu_handle toward the number of splits
 * \param typchk 1 if adaptation based on euclidean edge lengths
 * 2 if edge lengths computed in the metric
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the splitting operator and submit
 * anaelt task for this color to starPU.
 *
 * \remark For the splitting operator: there is a dependency between two
 * colors if they are adjacent.
 *
 * \todo Improvement: move temporary arrays outside the function to not
 * reallocate them at each call.
 */
int MMG2D_starpu_anaelt ( MMG5_pMesh mesh,starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_hash,
                          starpu_data_handle_t *handle_ns,
                          int typchk,int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_spldeps (mesh,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* For now, currunt color is not listed as a dependency; add it */
  task_handles[0].handle = handle_per_colors[color];
  task_handles[0].mode   = STARPU_W|STARPU_COMMUTE;
  /* Now, list the other deps */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&anaelt_codelet,
                               STARPU_SEQUENTIAL_CONSISTENCY,0,
                               STARPU_RW, *handle_mesh,
                               STARPU_RW, *handle_met,
                               STARPU_RW, *handle_hash,
                               STARPU_REDUX, *handle_ns,
                               //STARPU_DATA_MODE_ARRAY, task_handles, ndeps+1,
                               STARPU_VALUE, &typchk, sizeof(typchk),
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert:anaelt_codelet");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash point hash table listing the colors to which each point is
 * directly connected or connected through 1 edge.
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_per_color array to have 1 starpu handle per color
 * \param handle_nc starpu_handle toward the number of moves
 * \param typchk 1 if edge lengths are computed using euclidean norm, 2 if
 * computed using input metric.
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the collapse operator and submit
 * colelt task for this color to starPU.
 *
 * \remark For the collapse operator: there is a dependency between two
 * colors if an edge that connects points belonging to each colors exist.
 *
 */
int MMG2D_starpu_colelt ( MMG5_pMesh mesh,MMG5_HashP *hash,
                          starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_nc,
                          int typchk,int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_1edgdeps (mesh,hash,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* This time, current color is already listed in computed deps so we don't
   * need the 0 slot of the deps array */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&colelt_codelet,
                               STARPU_RW, *handle_mesh,
                               STARPU_RW, *handle_met,
                               STARPU_REDUX, *handle_nc,
                               STARPU_DATA_MODE_ARRAY, task_handles+1, ndeps,
                               STARPU_VALUE, &typchk, sizeof(typchk),
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_ns starpu_handle toward the number of splits
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the splitting operator and submit
 * adpspl task for this color to starPU.
 *
 * \remark For the splitting operator: there is a dependency between two
 * colors if they are adjacent.
 *
 * \todo Improvement: move temporary arrays outside the function to not
 * reallocate them at each call.
 */
int MMG2D_starpu_adpspl ( MMG5_pMesh mesh,starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_ns,
                          int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_spldeps (mesh,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* For now, currunt color is not listed as a dependency: add it */
  task_handles[0].handle = handle_per_colors[color];
  task_handles[0].mode   = STARPU_W|STARPU_COMMUTE;
  /* Now, list the other deps */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&adpspl_codelet,
                               STARPU_RW, *handle_mesh,
                               STARPU_RW, *handle_met,
                               STARPU_REDUX, *handle_ns,
                               STARPU_DATA_MODE_ARRAY, task_handles, ndeps+1,
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert:adpspl_codelet");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash point hash table (to fill by the list of colors of each point)
 *
 * \return \a 1 if success, 0 if fail.
 *
 * List the colors to which each vertex belongs (using a hash table)
 *
 * \todo A future improvement can be to build in // this hashtable and to
 * synchronize threads.
 *
 */
int MMG2D_pointColor(MMG5_pMesh mesh,MMG5_HashP *hash) {
  MMG5_pTria   pt;
  int          k, i;
  MMG5_hpoint  *ph;

  /** Step 1: Store list of colors to which each point belong (primary hash
   * table) */
  if ( !MMG5_hashPNew( mesh,hash,mesh->np, 3.01*mesh->np) ) return 0;

  for ( k=1; k<=mesh->nt; ++k ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for ( i=0; i<3; ++i ) {
      MMG5_hashPoint(mesh,hash,pt->v[i],pt->color1);
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash primary point hash table (contains the list of colors of each point)
 * \param hash secondary point hash table (to fill by the list of colors to
 * which each point is connected either directly, or through an edge)
 *
 * \return \a 1 if success, 0 if fail.
 *
 * List the colors to which each vertex belongs or is connected through 1 edge
 * (using a hash table)
 *
 * \todo A future improvement can be to build in // this hashtable and to
 * synchronize threads.
 *
 */
int MMG2D_pointColor_to_1edgColor(MMG5_pMesh mesh,MMG5_HashP *hash,MMG5_HashP *hash2) {
  MMG5_pTria   pt;
  int          k, i;
  MMG5_hpoint  *ph;

  if ( !MMG5_hashPNew( mesh,hash2,mesh->np, 6.01*mesh->np) ) return 0;

  for ( k=1; k<=mesh->nt; ++k ) {

    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for ( i=0; i<3; ++i ) {
      int ip1  = pt->v[i];
      int ip2  = pt->v[MMG5_inxt2[i]];

      assert ( hash->item[ip1].data && hash->item[ip2].data && "no color");

      /* Append colors of point 2 to point 1 */
      ph = &hash->item[ip2];
      while ( ph->nxt ) {
        MMG5_hashPoint(mesh,hash2,ip1,ph->data);
        ph = &hash->item[ph->nxt];
      }

      /* Append colors of point 1 to point 2 */
      ph = &hash->item[ip1];
      while ( ph->nxt ) {
        MMG5_hashPoint(mesh,hash2,ip2,ph->data);
        ph = &hash->item[ph->nxt];
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash->item);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash secondary point hash table (to fill by the list of colors to
 * which each point is connected either directly, or through an edge)
 *
 * \return \a 1 if success, 0 if fail.
 *
 * List the colors to which each vertex belongs or is connected through 1 edge
 * (using a hash table)
 *
 * \todo A future improvement can be to build in // this hashtable and to
 * synchronize threads.
 *
 */
int MMG2D_pointColor_1edg(MMG5_pMesh mesh,MMG5_HashP *hash2) {
  MMG5_HashP  hash;

  /** Step 1: Store list of colors to which each point belongs (primary hash
   * table) */
  if ( !MMG2D_pointColor(mesh,&hash) ) {
    fprintf(stderr,"  ## Problem in first step of dependencies construction"
            " for moving operator."
            "Unable to complete mesh. Exit program.\n");
    return 0;
  }

  /** Step 2: Append to each point the colors of points that are
   * connected to the current point by an edge (secondary hash table). */
  if ( !MMG2D_pointColor_to_1edgColor(mesh,&hash,hash2) ) {
    fprintf(stderr,"  ## Problem in second step of dependencies construction"
            " for moving operator."
            "Unable to complete mesh. Exit program.\n");
    return 0;
  }
  return 1;
}



/**
 * \param mesh pointer toward the mesh structure
 * \param hash2 point hash table listing the colors to which each point is
 * directly connected or connected through 1 edge.
 * \param deps array to store colors dependencies (has to be of size ncolors+1)
 * \param color color on which we compute deps
 *
 * \return \a ndeps, the number of dependencies for current color.
 *
 * Compute deps for move operator: for now each color travels the entire mesh to
 * build its deps.
 *
 */
int MMG2D_1edgdeps(MMG5_pMesh mesh,MMG5_HashP *hash2,int *deps,int color) {
  MMG5_pTria   pt;
  int          k, i;
  MMG5_hpoint  *ph;

  /** Step 1: Fill \a deps[i] by \a i if color \a color has a dependency with color \a
   * i */
  for ( k=1; k<=mesh->nt; ++k ) {

    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) ) continue;

    pt = &mesh->tria[k];
    for ( i=0; i<3; ++i ) {
      int ip1 = pt->v[i];
      ph = &hash2->item[ip1];
      deps[ph->data] = ph->data;
      while ( ph->nxt ) {
        ph = &hash2->item[ph->nxt];
        deps[ph->data] = ph->data;
      }
    }
  }

  /** Step 2: pack the array of dependencies and count number of dependencies */
  k     = 0;
  int ndeps = mesh->info.ncolors;

  /* Search for the last mark color */
  while ( !deps[ndeps] && ndeps > 0 ) --ndeps;

  while ( ++k < ndeps ) {
    if ( !deps[k] ) {
      /* copy dependency stored in last position in the empty place */
      deps[k]     = deps[ndeps];
      deps[ndeps] = 0;

      /* Search the last marked color */
      while ( !deps[ndeps] ) --ndeps;
    }
  }

  // debug
  /* printf(" Color %d -- %d deps:\n",color,ndeps); */
  /* for ( i=1; i<=ndeps; ++i ) { */
  /*   printf(" %d",deps[i]); */
  /* } */
  /* printf("\n"); */

  return ndeps;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash point hash table listing the colors to which each point is
 * directly connected or connected through 1 edge.
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_per_color array to have 1 starpu handle per color
 * \param handle_nc starpu_handle toward the number of moves
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the collapse operator and submit
 * adpcol task for this color to starPU.
 *
 * \remark For the collapse operator: there is a dependency between two
 * colors if an edge that connects points belonging to each colors exist.
 *
 */
int MMG2D_starpu_adpcol ( MMG5_pMesh mesh,MMG5_HashP *hash,
                          starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_nc,
                          int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_1edgdeps (mesh,hash,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* This time, current color is already listed in computed deps so we don't
   * need the 0 slot of the deps array */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&adpcol_codelet,
                               STARPU_RW, *handle_mesh,
                               STARPU_RW, *handle_met,
                               STARPU_REDUX, *handle_nc,
                               STARPU_DATA_MODE_ARRAY, task_handles+1, ndeps,
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash point hash table listing the colors to which each point is
 * directly connected or connected through 1 edge.
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_per_color array to have 1 starpu handle per color
 * \param handle_nsw starpu_handle toward the number of moves
 * \param typchk 1 if edge lengths are computed using euclidean norm, 2 if
 * computed using input metric.
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the swapping operator and submit
 * movtri task for this color to starPU.
 *
 * \remark For the swapping operator: there is a dependency between two
 * colors if an edge that connects points belonging to each colors exist.
 *
 */
int MMG2D_starpu_swpmsh ( MMG5_pMesh mesh,MMG5_HashP *hash,
                          starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_nsw,
                          int typchk,int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_1edgdeps (mesh,hash,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* This time, current color is already listed in computed deps so we don't
   * need the 0 slot of the deps array */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&swpmsh_codelet,
                               STARPU_RW, *handle_mesh,
                               STARPU_RW, *handle_met,
                               STARPU_REDUX, *handle_nsw,
                               STARPU_DATA_MODE_ARRAY, task_handles+1, ndeps,
                               STARPU_VALUE, &typchk, sizeof(typchk),
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param hash point hash table listing the colors to which each point is
 * directly connected or connected through 1 edge.
 * \param handle_mesh starpu handle toward the mesh structure
 * \param handle_met starpu handle toward the metric structure
 * \param handle_per_color array to have 1 starpu handle per color
 * \param handle_nm starpu_handle toward the number of moves
 * \param maxit_mov maximal number of iteration of moves.
 * \param improve 1 for last waves of move
 * \param color color to process
 *
 * \return 1 if success, 0 otherwise
 *
 * Compute dependencies of color \a color for the moving operator and submit
 * movtri task for this color to starPU.
 *
 * \remark For the moving operator: there is a dependency between two
 * colors if an edge that connects points belonging to each colors exist.
 *
 */
int MMG2D_starpu_movtri ( MMG5_pMesh mesh,MMG5_HashP *hash,
                          starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_nm,
                          int maxit_mov,int8_t improve,int color ) {

  starpu_data_handle_t handle_deps;

  int i,ndeps,*deps;

  /** Step 1: Find colors that have dependencies with current color */
  MMG5_SAFE_CALLOC(deps,mesh->info.ncolors+1,int,return 0);
  ndeps = MMG2D_1edgdeps (mesh,hash,deps,color);

  /** Step 2: Create task_handles list from dependencies */
  struct starpu_data_descr *task_handles;
  MMG5_SAFE_CALLOC(task_handles,ndeps+1,struct starpu_data_descr,return 0);

  /* This time, current color is already listed in computed deps so we don't
   * need the 0 slot of the deps array */
  for ( i=1; i<=ndeps; i++ ) {
    task_handles[i].handle = handle_per_colors[deps[i]];
    task_handles[i].mode   = STARPU_W|STARPU_COMMUTE;
  }
  MMG5_SAFE_FREE(deps);

  /** Step 3: Insert starpu task */
  int ret = starpu_task_insert(&movtri_codelet,
                               STARPU_R, *handle_mesh,
                               STARPU_R, *handle_met,
                               STARPU_REDUX, *handle_nm,
                               STARPU_DATA_MODE_ARRAY, task_handles+1, ndeps,
                               STARPU_VALUE, &maxit_mov, sizeof(maxit_mov),
                               STARPU_VALUE, &improve, sizeof(improve),
                               STARPU_VALUE, &color, sizeof(color),
                               0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  MMG5_SAFE_FREE(task_handles);

  return 1;
}

#endif
