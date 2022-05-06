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
 * \brief Mesh partitionning
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "partitionning.h"

/**
 * \param nelt number of nodes in graph (number of mesh elements)
 * \param xadj pointer toward the position of adjacents elements in adjncy
 * \param adjncy pointer toward the list of the adjtcent of each elt
 * \param npart number of partitions to create.
 * \param part output partition array
 *
 * \return  1 if success, 0 if fail
 *
 *  Mesh partitionning into \a npart partitions from already built graph.
 *
 * \remark the MMG_PART_INT type is defined by CMake to idx_t type with metis or
 * SCOTCH_num type with scotch.
 */
int MMG5_part_meshElts( MMG_PART_INT nelt,MMG_PART_INT*xadj,MMG_PART_INT*adjncy,
                        MMG_PART_INT npart,MMG_PART_INT*part) {
  int status = 1;

#ifdef USE_SCOTCH_PARTITIONNER
  SCOTCH_Graph      grafdat;

  CHECK_SCOTCH2 ( SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit", status=0);
  CHECK_SCOTCH2 ( SCOTCH_graphBuild (&grafdat, 0, nelt, xadj, NULL, NULL, NULL, xadj[nelt], adjncy, NULL),
                 "SCOTCH_graphBuild", status=0);

  CHECK_SCOTCH2 ( SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck",
                  fprintf(stderr,"\n  ## Error: %s: error detected by SCOTCH_grafCheck.",__func__);
                  return 0; )

  SCOTCH_Strat      grafstrat;
  CHECK_SCOTCH2 ( SCOTCH_stratInit (&grafstrat), "SCOTCH_stratInit", status=0 );

  /** Call scotch and get the partition array */
  CHECK_SCOTCH2 ( SCOTCH_graphPart (&grafdat, npart, &grafstrat, part),
                  "SCOTCH_graphPart", status=0 );

  SCOTCH_graphExit (&grafdat);
  SCOTCH_stratExit (&grafstrat);

#else
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;

  METIS_SetDefaultOptions(options);

  /** Call metis and get the partition array */
  if( npart >= 8 ) {

    /* Ask for partition contiguity */
    options[METIS_OPTION_CONTIG] = 1;
    CHECK_METIS ( METIS_PartGraphKway( &nelt,&ncon,xadj,adjncy,NULL,NULL, NULL,&npart,
                                       NULL,NULL,NULL,&objval, part ),
                  "METIS_PartGraphKway",status=0);
  }
  else
  {
    CHECK_METIS ( METIS_PartGraphRecursive( &nelt,&ncon,xadj,adjncy,NULL,NULL, NULL,&npart,
                                            NULL,NULL,NULL,&objval, part ),
                  "METIS_PartGraphRecursive",status=0);
  }
#endif

  return status;
}

#endif
