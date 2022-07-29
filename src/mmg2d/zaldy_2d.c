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
 * \file mmg2d/zaldy_2d.c
 * \brief Memory management.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */
#include "mmg2d.h"


/* Create a new vertex in the mesh, and return its number */
int MMG2D_newPt(MMG5_pMesh mesh,double c[2],int16_t tag, int color) {
  MMG5_pPoint  ppt;
  int          curpt;

// warning concurrency access lead to reset npnil (use of point curpt) while trying to access it, as a 0 npnil means that we don't have anymore memory, it raises a realloc issue
  // For now, we solve this using locks. These locks have to be removed once the
  // list will be parallelized
 // MMG5_LOCK(&mesh->lock); // Comment locks to use correctly the parallel newPt

#ifdef USE_STARPU
  int worker_id = starpu_worker_get_id();
  int max = mesh->np+mesh->npmax4t*(worker_id+1);
#else
  int worker_id = 0;
  int max = mesh->npmax;
#endif
  
  int prev_id = mesh->lastllpoint[color];
  if ( !mesh->npnil[worker_id] ) {
    return 0;
  }
  curpt = mesh->npnil[worker_id];
  if ( curpt > max) return 0;
  ppt = &mesh->point[curpt];
  memcpy(ppt->c,c,2*sizeof(double));
  ppt->tag &= ~MG_NUL;
  // is this needed?
  //mesh->npnil[starpu_worker_get_id()] = ppt->tmp;
  mesh->npnil[worker_id] = curpt + 1;
  ppt->tmp = 0;
  ppt->tag = tag;
  ppt->nxt = 0;
  ppt->prv = prev_id;
  mesh->point[prev_id].nxt = curpt;
  mesh->np4t[worker_id]++;

//#else
//  if ( !mesh->npnil[0] ) {
//    return 0;
//  }
//  curpt = mesh->npnil[0];
//  if ( mesh->npnil[0] > mesh->np )  mesh->np = mesh->npnil[0];
//  ppt = &mesh->point[curpt];
//  memcpy(ppt->c,c,2*sizeof(double));
//  ppt->tag &= ~MG_NUL;
//  mesh->npnil[0] = ppt->tmp;
//  ppt->tmp = 0;
//  ppt->tag = tag;
//#endif

  return curpt;
}

/* Delete a point in the mesh and update the garbage collector accordingly */
void MMG2D_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;

  ppt = &mesh->point[ip];
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag = MG_NUL;
#ifdef USE_STARPU
  int worker_id = starpu_worker_get_id();
#else
  int worker_id = 0;
#endif

  ppt->tmp = mesh->npnil[worker_id]; //useful?
  mesh->npnil[worker_id] = ip;
  int id_prv = mesh->point[ip].prv;
  int id_nxt = mesh->point[ip].nxt;
  mesh->point[id_prv].nxt = id_nxt;
  mesh->point[id_nxt].prv = id_prv;
  mesh->np4t[worker_id]--;
// #else
//  ppt->tmp = mesh->npnil[0];
//  mesh->npnil[0] = ip;
//  if ( ip == mesh->np )  mesh->np--;
// #endif
}

void MMG5_delEdge(MMG5_pMesh mesh,int iel) {
  MMG5_pEdge    pt;

  pt = &mesh->edge[iel];
  if ( !pt->a ) {
    fprintf(stdout,"  ## INVALID EDGE.\n");
    return;
  }
  memset(pt,0,sizeof(MMG5_Edge));
  pt->b = mesh->nanil[0];
  mesh->nanil[0] = iel;
  if ( iel == mesh->na )  mesh->na--;
}

/* Create a new triangle in the mesh and return its address */
int MMG2D_newElt(MMG5_pMesh mesh, int color) {
  int     curiel;

// warning concurrency access lead to reset ntnil (use of tria curiel) while trying to access it, as a 0 npnil means that we don't have anymore memory, it raises a realloc issue
  // For now, we solve this using locks. These locks have to be removed once the
  // list will be parallelized

#ifdef USE_STARPU
  int worker_id = starpu_worker_get_id();
  int max = mesh->nt+mesh->nemax4t*(worker_id+1);
#else
  int worker_id = 0;
  int max = mesh->ntmax;
#endif

  int prev_id = mesh->lastlltria[color];
  if ( !mesh->nenil[worker_id] ) {
    return 0;
  }
  curiel = mesh->nenil[worker_id];
  if ( curiel > max) return 0;
  mesh->nenil[worker_id] = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0; // can be deleted?
  mesh->tria[curiel].ref = 0;
  mesh->tria[curiel].color1 = color;
  mesh->tria[curiel].prv = prev_id;
  mesh->tria[curiel].nxt = 0;
  mesh->tria[prev_id].nxt = curiel;
  mesh->nenil[worker_id] = curiel + 1;
// #else
//  if ( !mesh->nenil[0] ) {
//    return 0;
//  }
//  curiel = mesh->nenil[0];
//  if ( mesh->nenil[0] > mesh->nt )  mesh->nt = mesh->nenil[0];
//  mesh->nenil[0] = mesh->tria[curiel].v[2];
//  mesh->tria[curiel].v[2] = 0;
//  mesh->tria[curiel].ref = 0;
// #endif
  mesh->tria[curiel].base = 0;
  mesh->tria[curiel].edg[0] = 0;
  mesh->tria[curiel].edg[1] = 0;
  mesh->tria[curiel].edg[2] = 0;
  mesh->ne4t[worker_id]++;
  return curiel;
}

/* Delete a triangle in the mesh and update the garbage collector accordingly */
int MMG2D_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;
  int           iadr;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT.\n");
    return 0;
  }
  memset(pt,0,sizeof(MMG5_Tria));

#ifdef USE_STARPU
  int worker_id = starpu_worker_get_id();
#else
  int worker_id = 0;
#endif

  pt->v[2] = mesh->nenil[worker_id];
  pt->qual = 0.0;
  iadr = (iel-1)*3 + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,3*sizeof(int));
  
  mesh->nenil[worker_id] = iel;
  pt->color1 = MG_NUL; //needed?
  int id_prv = mesh->tria[iel].prv;
  int id_nxt = mesh->tria[iel].nxt;
  mesh->tria[id_prv].nxt = id_nxt;
  mesh->tria[id_nxt].prv = id_prv;
  mesh->ne4t[worker_id]--;

//#else
//  pt->v[2] = mesh->nenil[0];
//  pt->qual = 0.0;
//  iadr = (iel-1)*3 + 1;
//  if ( mesh->adja )
//    memset(&mesh->adja[iadr],0,3*sizeof(int));
//
//  mesh->nenil[0] = iel;
//  if ( iel == mesh->nt )  mesh->nt--;
//#endif
  return 1;
}


/* check if n elets available */
int MMG5_getnElt(MMG5_pMesh mesh,int n) {
  int     curiel;

  if ( !mesh->nenil[0] )  return 0;
  curiel = mesh->nenil[0];
  do {
    curiel = mesh->tria[curiel].v[2];
  }
  while (--n);

  return n == 0;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * Initialize the linked lists of points and triangles based on their color.
 * This routine works wheather StarPU is used or not.
 *
 */
int MMG2D_linkInit(MMG5_pMesh mesh, int *part) {
  int i,j,idx;
  // triangles
  for (i=mesh->nt-1; i>=0; i--) {
    MMG5_pTria pt = &mesh->tria[i+1];
    if (!mesh->initlltria[part[i]]){ // if it's true, it means we are in the last element of the color part[i]
      pt->nxt = 0; // so we put its nxt field to 0
      mesh->lastlltria[part[i]] = i+1; // and the element is saved as last element of the chain of the color part[i]
    } else
      pt->nxt = mesh->initlltria[part[i]]; // if pt is not the last element we put as nxt the element previously saved of the same color. 
    mesh->initlltria[part[i]] = i+1; // we update the first element of the chain of color part[i].
    pt->color1 = part[i]+1;
    // vertices
    for (j=0; j<3; j++) {
      MMG5_pPoint ppt =&mesh->point[pt->v[j]];
      if (ppt->flag) continue;  // Francesco: this may be cause of error (before there was the check on the color
      if (!mesh->initllpoint[part[i]]){
        ppt->nxt = 0;
        mesh->lastllpoint[part[i]] = pt->v[j];
      } else
        ppt->nxt = mesh->initllpoint[part[i]];
      mesh->initllpoint[part[i]] = pt->v[j];
      ppt->flag++;
      // ppt->color1 = part[i]+1;
    }
  }
  for (i=0; i<mesh->np; i++)
    mesh->point[i].flag = 0;
  for (i=1; i<=mesh->info.ncolors; i++) {
    MMG5_pPoint ppt = &mesh->point[mesh->initllpoint[i-1]];
    MMG5_pTria  pt = &mesh->tria[mesh->initlltria[i-1]];
    ppt->prv = 0;
    pt->prv = 0;
    while (ppt->nxt){
      if (!ppt->prv) {
        idx = mesh->initllpoint[i-1];
      } else {
        idx = mesh->point[ppt->prv].nxt;
      }
      mesh->np4t[i-1]++;
      ppt = &mesh->point[ppt->nxt];
      ppt->prv = idx;
      // printf(" previous = %d\n current = %d, color =%d\n next = %d\n",ppt->prv,ppt->idx,ppt->color1,ppt->nxt);
    }
    mesh->np4t[i-1]++;
    while (pt->nxt){
      if (!pt->prv) {
        idx = mesh->initlltria[i-1];
      } else {
        idx = mesh->tria[ppt->prv].nxt;
      }
      mesh->ne4t[i-1]++;
      pt = &mesh->tria[pt->nxt];
      pt->prv = idx;
      // printf("previous = %d\n current = %d\n next = %d\n",pt->prv,pt->idx,pt->nxt);
    }
    mesh->ne4t[i-1]++;
    // printf("New color");
  }
}


/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * Set the memMax value to its "true" value (50% of the RAM or memory asked by
 * user) and perform memory repartition for the -m option.  If -m is not given,
 * memMax is the detected RAM. If -m is provided, check the user option and set
 * memMax to the available RAM if the user ask for too much memory. Last,
 * perform the memory repartition between the mmg arrays with respect to the
 * memMax value.
 *
 * \remark Here, mesh->npmax/ntmax must be setted.
 *
 */
static inline
int MMG2D_memOption_memSet(MMG5_pMesh mesh) {
  size_t   usedMem,avMem,reservedMem,npadd;
  int      ctri,bytes;

  MMG5_memOption_memSet(mesh);

  /** init allocation need MMG5_MEMMIN B */
  reservedMem = MMG5_MEMMIN + mesh->nquad*sizeof(MMG5_Quad);

  /** Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->nt+1)*sizeof(MMG5_Tria) + (3*mesh->nt+1)*sizeof(int)
    + (mesh->na+1)*sizeof(MMG5_Edge) + (mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %zu MB of memory ",__func__,
            mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu MB minimum\n",
            usedMem/MMG5_MILLION+1);
    return 0;
  }

  ctri = 2;

  /** Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
   * point+tria+edges+adjt+ aniso sol */
  bytes = sizeof(MMG5_Point) +
    2*sizeof(MMG5_Tria) + 3*2*sizeof(int)
    + 0.2*sizeof(MMG5_Edge) + 3*sizeof(double);

  avMem = mesh->memMax-usedMem;

  npadd = avMem/(2*bytes);
#ifdef USE_STARPU
  /* We want to use the max possible memory because realloc is not available */
  if ( mesh->info.mem < 0 ) {
    mesh->npmax = MG_MAX(mesh->npmax,mesh->np+npadd);
    mesh->ntmax = MG_MAX(mesh->ntmax,ctri*npadd+mesh->nt);
    mesh->namax = MG_MAX(mesh->namax,ctri*npadd+mesh->na);
  }
  else {
    mesh->npmax = MG_MAX(1.5*mesh->np,mesh->np+npadd);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,ctri*npadd+mesh->nt);
    mesh->namax = MG_MAX(mesh->na,ctri*npadd+mesh->na);
  }

#else
  /* Realloc is available, we don't want to allocate too large arrays in a first
   * guess */
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,ctri*npadd+mesh->nt);
  mesh->namax = MG_MIN(mesh->namax,ctri*npadd+mesh->na);
#endif
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (MB)    %zu\n",
            mesh->memMax/MMG5_MILLION);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MMG2D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  MMG2D_NTMAX    %d\n",mesh->ntmax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the -m option
 *
 */
int MMG2D_memOption(MMG5_pMesh mesh) {

  mesh->memMax = MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,MMG2D_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,MMG2D_NEMAX);
  mesh->namax = mesh->na;

  return  MMG2D_memOption_memSet(mesh);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh.
 *
 */
int MMG2D_setMeshSize_alloc( MMG5_pMesh mesh ) {
  int k;

#ifdef USE_STARPU
  pthread_mutex_init(&mesh->lock,NULL);

  /* StarPU configuration: set the sceduling policy */
  struct starpu_conf conf;
  starpu_conf_init(&conf);
  conf.sched_policy_name = "eager";

  /* StarPU initialization method */
  int ret = starpu_init(&conf);
  if (ret == -ENODEV) return 0;
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

  /* STARPU task profiling info */
  starpu_profiling_status_set(STARPU_PROFILING_ENABLE);

  /* Keep space so each thread can simulate operators in different memory slots. */
  int nbadd_pos = starpu_worker_get_count();
  mesh->nemax4t = (mesh->ntmax-mesh->nt)/nbadd_pos;
  mesh->npmax4t = (mesh->npmax-mesh->np)/nbadd_pos;
  mesh->namax4t = (mesh->namax-mesh->na)/nbadd_pos;

#else
  int nbadd_pos = 0;
#endif

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,(mesh->npmax+nbadd_pos)+1,MMG5_Point,return 0);
  /* Now shift pointer so 0 position allow to access to \a nbadd_pos. Previous
   * memory slots will be used by the threads to simulate operators */
  mesh->point += nbadd_pos;

  MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return 0);
  MMG5_SAFE_CALLOC(mesh->tria,(mesh->ntmax+nbadd_pos)+1,MMG5_Tria,return 0);

  /* Now shift pointer so 0 position allow to access to \a nbadd_pos. Previous
   * memory slots will be used by the threads to simulate operators */
  mesh->tria += nbadd_pos;

  if ( mesh->nquad ) {
    MMG5_ADD_MEM(mesh,(mesh->nquad+1)*sizeof(MMG5_Quad),"initial quadrilaterals",return 0);
    MMG5_SAFE_CALLOC(mesh->quadra,(mesh->nquad+1),MMG5_Quad,return 0);
  }

  mesh->namax = mesh->na;
  if ( mesh->na ) {
    MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"initial edges",return 0);
    MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge,return 0);
  }

#ifndef USE_STARPU
  MMG5_SAFE_CALLOC(mesh->npnil,1,int,return MMG5_STRONGFAILURE);
  MMG5_SAFE_CALLOC(mesh->nenil,1,int,return MMG5_STRONGFAILURE);
  MMG5_SAFE_CALLOC(mesh->nanil,1,int,return MMG5_STRONGFAILURE);
  /* keep track of empty links */
  mesh->npnil[0] = mesh->np + 1;
  mesh->nenil[0] = mesh->nt + 1;
  mesh->nanil[0] = 0;

  for (k=mesh->npnil[0]; k<mesh->npmax-1; k++) {
    /* Set tangent field of point to 0 */
    mesh->point[k].n[0] = 0;
    mesh->point[k].n[1] = 0;
    mesh->point[k].n[2] = 0;
    /* link */
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil[0]; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;
#else
    /* keep track of empty links */
  MMG5_SAFE_CALLOC(mesh->npnil,starpu_worker_get_count(),int,return MMG5_STRONGFAILURE);
  MMG5_SAFE_CALLOC(mesh->nenil,starpu_worker_get_count(),int,return MMG5_STRONGFAILURE);
  MMG5_SAFE_CALLOC(mesh->nanil,starpu_worker_get_count(),int,return MMG5_STRONGFAILURE);

  for(int i=0; i<starpu_worker_get_count(); i++) {
    mesh->npnil[i] = mesh->np + i*mesh->npmax4t + 1;
    mesh->nenil[i] = mesh->nt + i*mesh->nemax4t + 1;
    mesh->nanil[i] = 0;

    for (k=mesh->npnil[i]; k<mesh->npmax4t-1; k++) {
      /* Set tangent field of point to 0 */
      mesh->point[k].n[0] = 0;
      mesh->point[k].n[1] = 0;
      mesh->point[k].n[2] = 0;
      /* link */
      mesh->point[k].tmp  = k+1;
    }

    for (k=mesh->nenil[i]; k<mesh->nemax4t-1; k++)
      mesh->tria[k].v[2] = k+1;
  }
#endif


  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * allocate main structure
 *
 */
int MMG2D_zaldy(MMG5_pMesh mesh) {

  if ( !MMG2D_memOption(mesh) )  return 0;

  return  MMG2D_setMeshSize_alloc(mesh);
}
