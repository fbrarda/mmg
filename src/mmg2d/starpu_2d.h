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

/** Main functions */
void MMG2D_starpu_anaelt(void *buffers[], void *cl_arg);
void MMG2D_starpu_colelt(void *buffers[], void *cl_arg);
void MMG2D_starpu_swpmsh(void *buffers[], void *cl_arg);
void MMG2D_starpu_swpmsh(void *buffers[], void *cl_arg);
void MMG2D_starpu_adpspl(void *buffers[], void *cl_arg);
void MMG2D_starpu_adpcol(void *buffers[], void *cl_arg);
void MMG2D_starpu_movtri(void *buffers[], void *cl_arg);

/** Task dependencies computation */
void MMG2D_starpu_spldep(void *buffers[], void *cl_arg);

/** Codelets for main functions */
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
void accumulate_cpu(void *descr[], void *cl_arg);

extern struct starpu_codelet accumulate_codelet;
extern struct starpu_codelet izero_codelet;

/*
 *    Codelet to neutral element initializer
 */
extern struct starpu_codelet izero_codelet;


/** Dbg tool */
void print_cpu(void *descr[], void *cl_arg);
extern struct starpu_codelet print_codelet;

#endif
