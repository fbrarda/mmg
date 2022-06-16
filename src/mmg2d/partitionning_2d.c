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
#include "partitionning_2d.h"
#include "libmmg2d.h"

#ifdef USE_STARPU

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of adjacents elements in adjncy
 * \param adjncy pointer toward the list of the adjtcent of each elt
 *
 * \return  1 if success, 0 if fail
 *
 * Build the mesh graph with the mesh elements as graph nodes.
 *
 * \warning the mesh must be packed
 *
 */
int MMG2D_build_meshEltsGraph( MMG5_pMesh mesh,MMG_PART_INT **xadj,MMG_PART_INT **adjncy ) {
  MMG5_pTria   pt;
  MMG_PART_INT nadjncy;
  int          *adja;
  int          j,k,iadr,jel,count,nbAdj,ier;

  /** Step 1: mesh adjacency creation */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
            __func__);
    return 0;
  }

  /** Step 2: build the graph */
  /* allocate xadj */
  MMG5_SAFE_CALLOC((*xadj),mesh->nt+1,MMG_PART_INT,return 0);

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
  MMG5_SAFE_CALLOC((*adjncy), nadjncy, MMG_PART_INT, ier=0;);
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
 * \param mesh pointer toward the MMG5_pMesh structure
 *
 * \return  1 if success, 0 if fail
 *
 *  Mesh partitionning into mesh->ncolors partitions.
 *
 * \remark the MMG_PART_INT type is defined by CMake to idx_t type with metis or
 * SCOTCH_num type with scotch.
 */
int MMG2D_part_meshElts( MMG5_pMesh mesh )
{
  MMG_PART_INT *xadj,*adjncy;
  MMG_PART_INT nelt = mesh->nt;
  MMG_PART_INT npart= mesh->info.ncolors;
  MMG_PART_INT *part;

  xadj = adjncy = NULL;

  /** Build the graph */
  if ( !MMG2D_build_meshEltsGraph(mesh,&xadj,&adjncy ) )
    return 0;

  /** Mesh partitioning */
  /* Allocate partition table */
  MMG5_SAFE_CALLOC(part,nelt,MMG_PART_INT,return 0);

  int status = MMG5_part_meshElts( nelt,xadj,adjncy,npart,part);

  if ( status==1 ) {
    int i;
    /* Mesh partitionning succeed */
    // triangles
     for (i=mesh->nt-1; i>=0; i--) {
       MMG5_pTria pt = &mesh->tria[i+1];
       if (!mesh->initlltria[part[i]]){ 
	       pt->nxt = 0;
       	       mesh->lastlltria[part[i]] = i+1;
       } else
	       pt->nxt = mesh->initlltria[part[i]];
       mesh->initlltria[part[i]] = i+1;
       pt->color1 = part[i]+1;
       pt->idx = i+1;
       // vertices
       for (int j=0; j<3; j++) {
         MMG5_pPoint ppt =&mesh->point[pt->v[j]];
	 if (ppt->color1) continue;
         if (!mesh->initllpoint[part[i]]){
           ppt->nxt = 0;
       	   mesh->lastllpoint[part[i]] = pt->v[j];
	 } else
           ppt->nxt = mesh->initllpoint[part[i]];
         mesh->initllpoint[part[i]] = pt->v[j];
         ppt->color1 = part[i]+1;
         ppt->idx = pt->v[j];
       } 
     }
     for (i=1; i<=mesh->info.ncolors; i++) {
       MMG5_pPoint ppt = &mesh->point[mesh->initllpoint[i-1]];
       MMG5_pTria  pt = &mesh->tria[mesh->initlltria[i-1]];
       ppt->prv = 0;
       pt->prv = 0;
       while (ppt->nxt){
         int prv = ppt->idx;
	 ppt = &mesh->point[ppt->nxt];
	 ppt->prv = prv;
//	 printf(" previous = %d\n current = %d, color =%d\n next = %d\n",ppt->prv,ppt->idx,ppt->color1,ppt->nxt);
       }
       while (pt->nxt){
         int prv = pt->idx;
	 pt = &mesh->tria[pt->nxt];
	 pt->prv = prv;
	// printf("previous = %d\n current = %d\n next = %d\n",pt->prv,pt->idx,pt->nxt);
       }
      // printf("New color");
     }
  }

  /* Use environment variable to choose to save partition as it is a
   * developper+debug tool */
  if ( getenv("MMG_SAVE_PART") ) {
    /* Save partition */
    int i;
    for (i=0; i< mesh->nt; i++) {
      MMG5_pTria pt= &mesh->tria[i+1];
      part[i] = pt->ref;
      pt->ref = pt->color1;
    }

    char *basename = MMG5_Get_basename(mesh->namein);
    char *filename;

    MMG5_SAFE_CALLOC(filename,strlen(basename)+11,char,return 0);
    strncpy(filename,basename,strlen(basename));
    strcat(filename,"_part.mesh");

    MMG2D_saveMesh(mesh,filename);

    free(basename);basename=0;
    MMG5_SAFE_FREE(filename);

    for (i=0; i< mesh->nt; i++) {
      MMG5_pTria pt= &mesh->tria[i+1];
      pt->ref = part[i];
    }
  }

  /*deallocate xadj, adjncy et part */
  MMG5_DEL_MEM(mesh, adjncy);
  MMG5_DEL_MEM(mesh, xadj);
  MMG5_DEL_MEM(mesh, part);

  return status;
}

#endif
