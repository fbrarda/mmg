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
 * \file mmg2d/starpu_2d.h
 * \brief Mesh adaptation functions.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Mariem Makni (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 * Starpu codelets and functions.
 *
 */

#ifndef STARPU_2D_H

#define STARPU_2D_H

#include "starpu.h"
#include "mmg2d.h"

#define MMG2D_NCOLOR_MAX 1048576

/** Main functions */
int MMG2D_starpu_anaelt ( MMG5_pMesh mesh,starpu_data_handle_t *handle_mesh,
                          starpu_data_handle_t *handle_met,
                          starpu_data_handle_t *handle_per_colors,
                          starpu_data_handle_t *handle_hash,
                          starpu_data_handle_t *handle_ns,
                          int typchk,int color );


/** Task wrappers */
void MMG2D_hashTria_task(void *buffers[], void *cl_arg);
void MMG2D_anaelt_task(void *buffers[], void *cl_arg);
void MMG2D_colelt_task(void *buffers[], void *cl_arg);
void MMG2D_swpmsh_task(void *buffers[], void *cl_arg);
void MMG2D_swpmsh_task(void *buffers[], void *cl_arg);
void MMG2D_adpspl_task(void *buffers[], void *cl_arg);
void MMG2D_adpcol_task(void *buffers[], void *cl_arg);
void MMG2D_movtri_task(void *buffers[], void *cl_arg);

/** Task dependencies computation */
void MMG2D_spldep_task(void *buffers[], void *cl_arg);

/** Tools */
int MMG2D_spldeps ( MMG5_pMesh mesh,int *deps,int color);

/** Codelets for needed functions */
extern struct starpu_codelet colelt_codelet;
extern struct starpu_codelet swpmsh_codelet;
extern struct starpu_codelet anaelt_codelet;
extern struct starpu_codelet movtri_codelet;
extern struct starpu_codelet adpspl_codelet;
extern struct starpu_codelet adpcol_codelet;
extern struct starpu_codelet hashTria_codelet;

/** Codelets for dependencies */
extern struct starpu_codelet spldep_codelet;


/** For reduction */
void izero_cpu(void *descr[], void *cl_arg);
void ione_cpu(void *descr[], void *cl_arg);
void accumulate_cpu(void *descr[], void *cl_arg);
void min_cpu(void *descr[], void *cl_arg);

extern struct starpu_codelet accumulate_codelet;
extern struct starpu_codelet izero_codelet;
extern struct starpu_codelet min_codelet;
extern struct starpu_codelet ione_codelet;


/** Dbg tool */
void print_cpu(void *descr[], void *cl_arg);
extern struct starpu_codelet print_codelet;

#endif
