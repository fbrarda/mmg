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

/**
 * \file metis_MMG.c
 * \brief Partition mesh using metis
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "metis_mmg.h"

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjtcents in adjncy
 * \param adjncy pointer toward the list of the adjtcent of each elt
 *
 * \return  1 if success, 0 if fail
 *
 * Build the metis graph with the mesh elements as metis nodes.
 *
 * \warning the mesh must be packed
 *
 */
int MMG_graph_meshElts2metis( MMG5_pMesh mesh,idx_t **xadj,idx_t **adjncy ) {
  MMG5_pTria pt;
  idx_t      nadjncy;
  int        *adja;
  int        j,k,iadr,jel,count,nbAdj,ier;

  /** Step 1: mesh adjacency creation */
  /*if ( (!mesh->adja) && (1 != MMG2D_hashTria(mesh) ) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to create "
    "adjacency table.\n",__func__);
    return 0;
    }*/

  /* create tria adjacency */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
            __func__);
    return 0;
  }

  /** Step 2: build the metis graph */
  /* allocate xadj */
  MMG5_SAFE_CALLOC((*xadj),mesh->nt+1,idx_t,return 0);

  /** 1) Count the number of adjacent of each elements and fill xadj */
  (*xadj)[0] = 0;
  nadjncy = 0;
  for( k = 1; k <= mesh->nt; k++ ) {
    nbAdj = 0;
    adja = &mesh->adja[3*(k-1) + 1];

    for( j = 0; j < 3; j++ ) {
      if ( adja[j] ) {
        nbAdj++;
      }
    }

    nadjncy+= nbAdj;
    (*xadj)[k] = nadjncy;
  }

  /** 2) List the adjacent of each elts in adjncy */
  ier = 1;
  ++nadjncy;
  MMG5_SAFE_CALLOC((*adjncy), nadjncy, idx_t, ier=0;);
  if( !ier ) {
    MMG5_DEL_MEM(mesh, (*xadj) );
    MMG5_DEL_MEM(mesh, (*adjncy));
    return ier;
  }

  count = 0;
  for( k = 1; k <= mesh->nt; k++ ) {
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[3*(k-1) + 1];
    pt   = &mesh->tria[k];
    for ( j = 0; j < 3; j++ ) {
      jel = adja[j] / 3;
      if ( !jel ) continue;

      if (adja[j]/3)
      {
        (*adjncy)[count]   = (adja[j]/3)-1;
        count++;
      }
    }
  }

  return ier;
}

/**
 * \param mesh pointer toward the pmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param nproc number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the first mesh into nprocs
 *
 */
int MMG_part_meshElts2metis( MMG5_pMesh mesh, idx_t* part, idx_t nproc )
{
  idx_t      *xadj,*adjncy,*vwgt;
  idx_t      nelt = mesh->nt;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier = 0;
  int        status = 1;

  xadj = adjncy = vwgt = NULL;

  METIS_SetDefaultOptions(options);

  /** Build the graph */
  if ( !MMG_graph_meshElts2metis(mesh,&xadj,&adjncy ) )
    return 0;

  /** Call metis and get the partition array */
  if( nproc >= 8 ) {
    options[METIS_OPTION_CONTIG] = 1;
    ier = METIS_PartGraphKway( &nelt,&ncon,xadj,adjncy,NULL,NULL, NULL,&nproc,
                               NULL,NULL,NULL,&objval, part );
  }
  else
  {
    ier = METIS_PartGraphRecursive( &nelt,&ncon,xadj,adjncy,NULL,NULL, NULL,&nproc,
                                    NULL,NULL,NULL,&objval, part );
  }

  if ( ier != METIS_OK ) {
    switch ( ier ) {
    case METIS_ERROR_INPUT:
      fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );
      break;
    case METIS_ERROR_MEMORY:
      fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" );
      break;
    case METIS_ERROR:
      fprintf(stderr, "METIS_ERROR: generic error\n" );
      break;
    default:
      fprintf(stderr, "METIS_ERROR: update your METIS error handling\n" );
      break;
    }
    status = 0;
  }


  /*deallocate xadj et adjncy */
  MMG5_DEL_MEM(mesh, adjncy);
  MMG5_DEL_MEM(mesh, xadj);

  return status;
}
