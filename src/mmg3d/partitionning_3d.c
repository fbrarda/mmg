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
 * \brief Mesh partitionning
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "partitionning_3d.h"

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjtcents in adjncy
 * \param adjncy pointer toward the list of the adjtcent of each elt
 *
 * \return  1 if success, 0 if fail
 *
 * Build mesh graph with the mesh elements as graph nodes.
 *
 * \warning the mesh must be packed
 *
 */
int MMG3D_build_meshEltsGraph( MMG5_pMesh mesh,MMG_PART_INT **xadj,MMG_PART_INT **adjncy ) {
  MMG5_pTetra  pt;
  MMG_PART_INT nadjncy;
  int          *adja;
  int          j,k,iadr,jel,count,nbAdj,ier;

  /** Step 1: mesh adjacency creation */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  /** Step 2: build the metis graph */
  /* allocate xadj */
  MMG5_SAFE_CALLOC((*xadj),mesh->ne+1,MMG_PART_INT,return 0);

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
  MMG5_SAFE_CALLOC((*adjncy), nadjncy, MMG_PART_INT, ier=0;);
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

      if (adja[j]/4) {
        (*adjncy)[count]   = (adja[j]/4)-1;
        count++;
      }
    }
  }

  return ier;
}

/**
 * \param mesh pointer toward the pmesh structure
 *
 * \return  1 if success, 0 if fail
 *
 *  Mesh partitionning into mesh->ncolors partitions.
 *
 */
int MMG3D_part_meshElts( MMG5_pMesh mesh)
{
  MMG_PART_INT       *xadj,*adjncy,*vwgt;
  MMG_PART_INT       nelt = mesh->ne;
  MMG_PART_INT       npart= mesh->info.ncolors;
  MMG_PART_INT       *part;
  int                i;

  xadj = adjncy = vwgt = NULL;

  /** Build the graph */
  if ( !MMG3D_build_meshEltsGraph(mesh,&xadj,&adjncy ) )
    return 0;

  /** Mesh partitioning */
  int status = MMG5_part_meshElts( nelt,xadj,adjncy,npart,part);

 if ( status==1 ) {
    int i;
    for (i=0; i< mesh->ne; i++) {
      MMG5_pTetra pt= &mesh->tetra[i+1];
      pt->color1 = part[i];
    }
 }

  /*deallocate xadj, adjncy et part */
  MMG5_DEL_MEM(mesh, adjncy);
  MMG5_DEL_MEM(mesh, xadj);
  MMG5_DEL_MEM(mesh, part);

  return status;
}
