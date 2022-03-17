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
#include "metis_mmg3d.h"

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
  MMG5_pTetra pt;
  idx_t      nadjncy;
  int        *adja;
  int        j,k,iadr,jel,count,nbAdj,ier;

  /** Step 1: mesh adjacency creation */
  /*if ( (!mesh->adja) && (1 != MMG3D_hashTetra(mesh) ) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to create "
    "adjacency table.\n",__func__);
    return 0;
    }*/

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  /** Step 2: build the metis graph */
  /* allocate xadj */
  MMG5_SAFE_CALLOC((*xadj),mesh->ne+1,idx_t,return 0);

  /** 1) Count the number of adjacent of each elements and fill xadj */
  (*xadj)[0] = 0;
  nadjncy = 0;
  for (k=1; k<=mesh->ne; k++) {
    nbAdj = 0;
    adja = &mesh->adja[4*(k-1) + 1];

    for( j = 0; j < 4; j++ ) {
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
  for (k=1; k<=mesh->ne; k++) {
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[4*(k-1) + 1];
    pt   = &mesh->tetra[k];
    for ( j=0; j<4; j++ ) {
      jel = adja[j] / 4;
      if ( !jel ) continue;

      if (adja[j]/4)
      {
        (*adjncy)[count]   = (adja[j]/4)-1;
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
int MMG3D_part_meshElts2metis( MMG5_pMesh mesh)
{
  idx_t      *xadj,*adjncy,*vwgt;
  idx_t      nelt = mesh->ne;
  idx_t*     part;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier = 0;
  int        status = 1;
  int        i;
  MMG5_pTetra pt;

  xadj = adjncy = vwgt = NULL;

  METIS_SetDefaultOptions(options);

  fprintf(stdout,"----nbelts= %d", mesh->ne);
  /** Build the graph */
  if ( !MMG_graph_meshElts2metis(mesh,&xadj,&adjncy ) )
    return 0;

  idx_t nproc=mesh->info.ncolors;
  
  /* Allocate the table part */
  MMG5_SAFE_CALLOC(part,mesh->ne,idx_t,return 0);
  
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
  
  for (i=0; i< mesh->ne; i++) {
    pt= &mesh->tetra[i+1];
    pt->color1 = part[i];
    pt->ref = pt->color1;
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
 
  /*deallocate xadj, adjncy et part */
  MMG5_DEL_MEM(mesh, adjncy);
  MMG5_DEL_MEM(mesh, xadj);
  MMG5_DEL_MEM(mesh, part);

  return status;
}
