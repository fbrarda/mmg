/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

#ifdef USE_STARPU

/**
 * \brief definitions for mesh partitionning.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef PARTITIONNING_H

#define PARTITIONNING_H

#include "mmgcommon.h"

#ifdef USE_SCOTCH_PARTITIONNER
#include <scotch.h>

#define CHECK_SCOTCH2(t,m,e) if(0!=t){perror(m);e;}

#else
#include <metis.h>

#define CHECK_METIS(t,m,e) do                                           \
  {                                                                     \
  int metis_ier_private = t;                                            \
  if ( metis_ier_private != METIS_OK ) {                                \
    switch ( metis_ier_private ) {                                      \
    case METIS_ERROR_INPUT:                                             \
      fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );        \
      break;                                                            \
    case METIS_ERROR_MEMORY:                                            \
      fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" ); \
      break;                                                            \
    case METIS_ERROR:                                                   \
      fprintf(stderr, "METIS_ERROR: generic error\n" );                 \
      break;                                                            \
    default:                                                            \
      fprintf(stderr, "METIS_ERROR: update your METIS error handling\n" ); \
      break;                                                            \
    }                                                                   \
    perror(m);                                                          \
    e;                                                                  \
  }                                                                     \
  } while(0)

#endif

int MMG5_part_meshElts( MMG_PART_INT,MMG_PART_INT*,MMG_PART_INT*,MMG_PART_INT,MMG_PART_INT*);

#endif

#endif
