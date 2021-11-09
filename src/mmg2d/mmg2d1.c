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
 * \file mmg2d/mmg2d1.c
 * \brief Mesh adaptation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"
#include <starpu.h>
#include "metis_mmg.h"
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include "libmmg2d.h"


struct starpu_codelet colelt_codelet =
{
        .cpu_funcs = {MMG2D_starpu_colelt},
        .cpu_funcs_name = {"MMG2D_starpu_colelt"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "colelt"
};


struct starpu_codelet swpmsh_codelet =
{
        .cpu_funcs = {MMG2D_starpu_swpmsh},
        .cpu_funcs_name = {"MMG2D_starpu_swpmsh"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "swpmsh"
};


struct starpu_codelet anaelt_codelet =
{
        .cpu_funcs = {MMG2D_starpu_anaelt},
        .cpu_funcs_name = {"MMG2D_starpu_anaelt"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "anaelt"
};

struct starpu_codelet movtri_codelet =
{
        .cpu_funcs = {MMG2D_starpu_movtri},
        .cpu_funcs_name = {"MMG2D_starpu_movtri"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "movtri"
};

struct starpu_codelet adpspl_codelet =
{
        .cpu_funcs = {MMG2D_starpu_adpspl},
        .cpu_funcs_name = {"MMG2D_starpu_adpspl"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "adpspl"
};

struct starpu_codelet adpcol_codelet =
{
        .cpu_funcs = {MMG2D_starpu_adpcol},
        .cpu_funcs_name = {"MMG2D_starpu_adpcol"},
        .nbuffers = 3,
        .modes = {STARPU_RW, STARPU_RW, STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU, STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "adpcol"
};

struct starpu_codelet hashTria_codelet =
{
        .cpu_funcs = {MMG2D_starpu_hashTria},
        .cpu_funcs_name = {"MMG2D_starpu_hashTria"},
        .nbuffers = 1,
        .modes = {STARPU_RW}, 
        .specific_nodes = 1,
        .nodes = {STARPU_SPECIFIC_NODE_CPU},
        .where = STARPU_CPU,
        .name = "hashTria"
};


/* Mesh adaptation routine for the first stages of the algorithm: intertwine splitting
 based on patterns, collapses and swaps.
   typchk = 1 -> adaptation based on edge lengths
   typchk = 2 -> adaptation based on lengths calculated in metric met */
int MMG2D_anatri(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {

  int      it,maxit,ns,nc,nsw,nns,nnc,nnsw;

  nns = nnc = nnsw = 0;
  ns=0; nc=0; nsw=0;
  it = 0;
  maxit = 5;
  int ret;
  
  int color, i;
  int SIZE=1;
  int vect_ns[SIZE];
  int vect_nc[SIZE];
  int vect_nsw[SIZE];
  
  
  starpu_data_handle_t vector_mesh, vector_met, vector_ns, vector_nc, vector_nsw;
  
   for (i=0; i< SIZE; i++)
      {
        vect_ns[0]=ns;
        vect_nc[0]=nc;
        vect_nsw[0]=nsw;

       }
  
 
  starpu_vector_data_register(&vector_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, 1, sizeof(MMG5_pMesh)); 
  starpu_vector_data_register(&vector_met, STARPU_MAIN_RAM, (uintptr_t)met, 1, sizeof(MMG5_pSol));
  starpu_vector_data_register(&vector_ns, STARPU_MAIN_RAM, (uintptr_t)vect_ns, SIZE, sizeof(vect_ns[0]));
  starpu_vector_data_register(&vector_nc, STARPU_MAIN_RAM, (uintptr_t)vect_nc, SIZE, sizeof(vect_nc[0]));
  starpu_vector_data_register(&vector_nsw, STARPU_MAIN_RAM, (uintptr_t)vect_nsw, SIZE, sizeof(vect_nsw[0]));
  
  /* Main routine; intertwine split, collapse and swaps */
  do {
    if ( typchk == 2 && it == 0 )  mesh->info.fem = 1;

    if ( !mesh->info.noinsert ) {
      /* Memory free */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;
     
      /* Split long edges according to patterns */
     //ns = MMG2D_anaelt(mesh,met,typchk, color);
     
     
     fprintf(stdout," vector_ns_1----------= %d .\n", vect_ns[0]);
     
     for (color=0; color< 1; color++)
     {
      ret = starpu_task_insert(&anaelt_codelet,
                     STARPU_RW, vector_mesh, 
                     STARPU_RW, vector_met, 
                     STARPU_RW, vector_ns, 
                     STARPU_VALUE, &typchk, sizeof(typchk), 
		      STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      }
      
      starpu_task_wait_for_all();
            
      fprintf(stdout," vector_ns_11111= %d .\n", vect_ns[0]);
      
      ns=vect_ns[0];
      if ( ns < 0 ) {
        fprintf(stderr,"  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }

      /* Recreate adjacencies */
      /*if ( !MMG2D_hashTria(mesh) ) {
        fprintf(stdout,"  ## Hashing problem. Exit program.\n");
        return 0;
      }*/
      //for (color=0; color< mesh->ncolors; color++)
      //{
      ret = starpu_task_insert(&hashTria_codelet,
                     STARPU_RW, vector_mesh, 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      //} 
       
      /* Collapse short edges */    
      //nc = MMG2D_colelt(mesh,met,typchk,color);

     for (color=0; color< 1; color++)
      {
      ret = starpu_task_insert(&colelt_codelet,
                     STARPU_RW, vector_mesh, 
                     STARPU_RW, vector_met, 
                     STARPU_RW, vector_nc, 
                     STARPU_VALUE, &typchk, sizeof(typchk), 
		      STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");  
      }
      starpu_task_wait_for_all();
      
      fprintf(stdout," vector_nc_1= %d .\n", vect_nc[0]);
      
      nc=vect_nc[0];      
      
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Unable to collapse mesh. Exiting.\n");
        return 0;
      }

    }
    else {
      ns = 0;
      nc = 0;
    }

    /* Swap edges */
    if ( !mesh->info.noswap ) {
     // nsw = MMG2D_swpmsh(mesh,met,typchk,color);
      
     for (color=0; color< 1; color++)
     {
      ret = starpu_task_insert(&swpmsh_codelet,
                     STARPU_RW, vector_mesh, 
                     STARPU_RW, vector_met, 
                     STARPU_RW, vector_nsw, 
                     STARPU_VALUE, &typchk, sizeof(typchk), 
		      STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      }
      starpu_task_wait_for_all();
      fprintf(stdout," vector_nsw_1= %d .\n", vect_nsw[0]);
      nsw=vect_nsw[0];      
      
      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else nsw = 0;

    nns += ns;
    nnc += nc;
    nnsw += nsw;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nsw);
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nsw >0 );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnsw,it);
  }
  
  starpu_data_unregister(vector_mesh);
  starpu_data_unregister(vector_met);
  starpu_data_unregister(vector_ns);
  starpu_data_unregister(vector_nc);
  starpu_data_unregister(vector_nsw);
  return 1;
}


void MMG2D_starpu_anaelt(void *buffers[], void *cl_arg) {
 
  int nx_mesh, nx_met, nx_ns;
  struct starpu_vector_interface *vect_mesh, *vect_met,*vect_ns ;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  int8_t typchk;
  int color,i;
  int *val_ns;
  //int ns;

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  vect_ns = (struct starpu_vector_interface *) buffers[2];
  nx_ns = STARPU_VECTOR_GET_NX(vect_ns);
  val_ns = (int *)STARPU_VECTOR_GET_PTR(vect_ns);
               
  //ns=val_ns[0];      
    
  starpu_codelet_unpack_args(cl_arg, &typchk, &color);
  //fprintf(stdout, "Hello world, the id_task is %d\n \n \n", color);
    
  //fprintf(stdout, "typchk= %d\n \n \n", typchk);
 
  val_ns[0]=MMG2D_anaelt(mesh,met,typchk,color);
  
  /*if ( ns < 0 ) {
     fprintf(stderr,"  ## Unable to anaelt mesh. Exiting.\n");
  }*/
  

}
/* Travel triangles and split long edges according to patterns */
int MMG2D_anaelt(MMG5_pMesh mesh,MMG5_pSol met,int typchk,int color) {

  
  MMG5_pTria      pt;
  MMG5_pPoint     ppt,p1,p2;
  MMG5_Hash       hash;
  double          len,s,o[2],no[2];
  int             ns,nc,npinit,ni,k,nt,ip1,ip2,ip,it,vx[3];
  int8_t          i,ic,i1,i2,ier;
  static int8_t   mmgWarn0=0;
  
  s = 0.5;
  ns = 0;
  npinit = mesh->np;

  if ( !MMG5_hashNew(mesh,&hash,mesh->np,3*mesh->np) ) return 0;

  /* Step 1: travel mesh, check edges, and tag those to be split; create the new vertices in hash */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || (pt->ref < 0) ) continue;
    if ( MG_SIN(pt->tag[0]) || MG_SIN(pt->tag[1]) || MG_SIN(pt->tag[2]) )  continue;

    /* Check if pt should be cut */
    pt->flag = 0;

    /* typchk=1 : split on geometric basis and rough size considerations */
    if ( typchk == 1) {
      if ( !MMG2D_chkedg(mesh,k) ) continue;
    }
    /* typchk =2 : split on basis of edge lengths in the metric */
    else if ( typchk ==2 ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( len > MMG2D_LLONG ) MG_SET(pt->flag,i);
      }
    }

    /* mesh->info.fem : split edges which are not MG_BDY, but whose vertices are both MG_BDY */
    if ( mesh->info.fem ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];
        if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) && !(pt->tag[i] & MG_BDY) ) MG_SET(pt->flag,i);
      }
    }

    if ( !pt->flag ) continue;
    ns++;

    /* Travel edges to split and create new points */
    for (i=0; i<3; i++) {
      if ( !MG_GET(pt->flag,i) ) continue;
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip = MMG5_hashGet(&hash,ip1,ip2);
      if ( ip > 0 ) continue;

      /* Geometric attributes of the new point */
      ier = MMG2D_bezierCurv(mesh,k,i,s,o,no);
      if ( !ier ) {
        MG_CLR(pt->flag,i);
        continue;
      }
      ip = MMG2D_newPt(mesh,o,pt->tag[i]);
      if ( !ip ) {
        /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            do {
                              MMG2D_delPt(mesh,mesh->np);
                            } while ( mesh->np>npinit );return -1;,
                            o,pt->tag[i]);
      }
      ppt = &mesh->point[ip];
      if ( MG_EDG(pt->tag[i]) ) {
        ppt->n[0] = no[0];
        ppt->n[1] = no[1];
      }

      /* If there is a metric in the mesh, interpolate it at the new point */
      assert ( met );
      if ( met->m )
        MMG2D_intmet(mesh,met,k,i,ip,s);

      /* Add point to the hashing structure */
      MMG5_hashEdge(mesh,&hash,ip1,ip2,ip);
    }
  }
  if ( !ns ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return ns;
  }

  /* Step 2: Make flags at triangles consistent between themselves (check if adjacent triangle is split) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;
    else if ( pt->flag == 7 ) continue;
    nc = 0;

    for (i=0; i<3; i++) {
      i1 = MMG5_iprv2[i];
      i2 = MMG5_inxt2[i];
      if ( !MG_GET(pt->flag,i) && !MG_SIN(pt->tag[i]) ) {
        ip = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( ip > 0 ) {
          MG_SET(pt->flag,i);
          nc++;
        }
      }
    }
    if ( nc > 0 ) ns++;
  }
  if ( mesh->info.ddebug && ns ) {
    fprintf(stdout,"     %d analyzed  %d proposed\n",mesh->nt,ns);
    fflush(stdout);
  }

  /* Step 3: Simulate splitting and delete points leading to an invalid configuration */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  nc = 0;
  do {
    ni = 0;
    for ( k=1; k<= mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;
      else if ( pt->flag == 0 ) continue;

      vx[0] = vx[1] =vx[2] = 0;
      pt->flag = 0;
      ic = 0;

      for (i=0; i<3; i++) {
        i1 = MMG5_iprv2[i];
        i2 = MMG5_inxt2[i];
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( vx[i] > 0 ) {
          MG_SET(pt->flag,i);
          if ( mesh->point[vx[i]].flag > 2 )  ic = 1;
        }
      }

      if ( !pt->flag )  continue;
      switch (pt->flag) {
        case 1: case 2: case 4:
          ier = MMG2D_split1_sim(mesh,met,k,vx);
          break;
        case 7:
          ier = MMG2D_split3_sim(mesh,met,k,vx);
          break;
        default:
          ier = MMG2D_split2_sim(mesh,met,k,vx);
          break;
      }
      if ( ier )  continue;

      /* An edge is invalidated in the process */
      ni++;
      if ( ic == 0 && MMG2D_dichoto(mesh,met,k,vx) ) {
        for (i=0; i<3; i++)
          if ( vx[i] > 0 )  mesh->point[vx[i]].flag++;
      }
      /* Relocate point at the center of the edge */
      else {
        for (i=0; i<3; i++) {
          if ( vx[i] > 0 ) {
            p1 = &mesh->point[pt->v[MMG5_iprv2[i]]];
            p2 = &mesh->point[pt->v[MMG5_inxt2[i]]];
            ppt = &mesh->point[vx[i]];
            ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
            ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
          }
        }
      }
    }
    nc += ni;
  }
  while ( ni >0 && ++it <20 );

  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d corrected,  %d invalid\n",nc,ni);
    fflush(stdout);
  }

  /* step 4: effective splitting */
  ns = 0;
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || pt->ref < 0 )  continue;
    else if ( pt->flag == 0 )  continue;

    vx[0] = vx[1] = vx[2] = 0;
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( MG_GET(pt->flag,i) ) {
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( !vx[i] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: unable to create point on"
                    " at least 1 edge.\n Exit program.\n",__func__);
          }
          return -1;
        }
      }
    }
    if ( pt->flag == 1 || pt->flag == 2 || pt->flag == 4 ) {
      ier = MMG2D_split1(mesh,met,k,vx);
      ns++;
    }
    else if ( pt->flag == 7 ) {
      ier = MMG2D_split3(mesh,met,k,vx);
      ns++;
    }
    else {
      ier = MMG2D_split2(mesh,met,k,vx);
      ns++;
    }
    if ( !ier ) return -1;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);
  MMG5_DEL_MEM(mesh,hash.item);

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param vx pointer toward table of edges to split.
 * \return 1.
 *
 * Find acceptable position for splitting.
 *
 */
int MMG2D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx) {
  MMG5_pTria   pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[3][2],p[3][2];
  float        to,tp,t;
  int          ia,ib,ier,it,maxit;
  int8_t       i,i1,i2;

  pt = &mesh->tria[k];

  /* Get point on curve and along segment for edge split */
  for (i=0; i<3; i++) {
    memset(p[i],0,2*sizeof(double));
    memset(o[i],0,2*sizeof(double));
    if ( vx[i] > 0 ) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      ia = pt->v[i1];
      ib = pt->v[i2];
      pa = &mesh->point[ia];
      pb = &mesh->point[ib];
      ps = &mesh->point[vx[i]];
      o[i][0] = 0.5 * (pa->c[0] + pb->c[0]);
      o[i][1] = 0.5 * (pa->c[1] + pb->c[1]);
      p[i][0] = ps->c[0];
      p[i][1] = ps->c[1];
    }
  }
  maxit = 4;
  it = 0;
  tp = 1.0;
  to = 0.0;

  do {
    /* compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
    switch (pt->flag) {
      case 1: case 2: case 4:
        ier = MMG2D_split1_sim(mesh,met,k,vx);
        break;
      case 7:
        ier = MMG2D_split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = MMG2D_split2_sim(mesh,met,k,vx);
        break;
    }
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );

  /* restore coords of last valid pos. */
  if ( !ier ) {
    t = to;
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
  }
  return 1;
}


void MMG2D_starpu_colelt(void *buffers[], void *cl_arg) {
 
  int nx_mesh, nx_met,nx_nc;
  struct starpu_vector_interface *vect_mesh, *vect_met, *vect_nc;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  int8_t typchk;
  int color;
  
  int *val_nc;
  //int nc;
 
  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  vect_nc = (struct starpu_vector_interface *) buffers[2];
  nx_nc = STARPU_VECTOR_GET_NX(vect_nc);
  val_nc = (int *)STARPU_VECTOR_GET_PTR(vect_nc);
     
  //nc=val_nc[0];  
  starpu_codelet_unpack_args(cl_arg, &typchk, &color);
 
  val_nc[0]=MMG2D_colelt(mesh,met,typchk,color);

 }

/* Travel triangles and collapse short edges */
int MMG2D_colelt(MMG5_pMesh mesh,MMG5_pSol met,int typchk,int color) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,ll,hmin2;
  int          list[MMG2D_LONMAX+2],ilist,nc,k;
  uint8_t      i,i1,i2,open;

  nc = 0;
  hmin2 = mesh->info.hmin * mesh->info.hmin;
 
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;

    /* Travel 3 edges of the triangle and decide whether to collapse p1->p2, based on length criterion */
    pt->flag = 0; // was here before, but probably serves for nothing

    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;

      /* Impossible to collapse a surface point onto a non surface point -- impossible to
       collapse a surface point along a non geometric edge */
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      /* Same test for REF points */
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      open = (mesh->adja[3*(k-1)+1+i] == 0) ? 1 : 0;

      /* Check length */
      if ( typchk == 1 ) {
        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        ll = ux*ux + uy*uy;
        if ( ll > hmin2 ) continue;
      }
      else {
        ll = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( ll > MMG2D_LSHRT ) continue;
      }

      /* Check whether the geometry is preserved */
      ilist = MMG2D_chkcol(mesh,met,k,i,list,typchk);

      if ( ilist > 3 || ( ilist == 3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     %8d vertices removed\n",nc);

  return nc;
}

/* Travel triangles and swap edges to improve quality */
void MMG2D_starpu_swpmsh(void *buffers[], void *cl_arg) {
 
  int nx_mesh, nx_met, nx_nsw;
  struct starpu_vector_interface *vect_mesh, *vect_met, *vect_nsw;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  int8_t typchk;
  int rett;
  int color;
  
  int *val_nsw;
  //int nsw;
  
  fprintf(stdout," hellooooo swpm \n");

  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  vect_nsw = (struct starpu_vector_interface *) buffers[2];
  nx_nsw = STARPU_VECTOR_GET_NX(vect_nsw);
  val_nsw = (int *)STARPU_VECTOR_GET_PTR(vect_nsw);
          
  //nsw=val_nsw[0];
  starpu_codelet_unpack_args(cl_arg, &typchk, &color);
  
  fprintf(stdout," typchk= %d  \n", typchk);
  
  val_nsw[0]=MMG2D_swpmsh(mesh,met,typchk,color);
  
  
  }

int MMG2D_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,int typchk, int color) {

  MMG5_pTria pt;
  int        it,maxit,ns,nns,k;
  uint8_t    i;

  it = nns = 0;
  maxit = 2;
  mesh->base++;

  do {
    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;

      for (i=0; i<3; i++) {
        if ( MG_SIN(pt->tag[i]) || MG_EDG(pt->tag[i]) ) continue;
        else if ( MMG2D_chkswp(mesh,met,k,i,typchk) ) {
          ns += MMG2D_swapar(mesh,k,i);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ns > 0 && ++it<maxit );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return nns;
}


/* Mesh adaptation routine for the final stage of the algorithm: intertwine splitting
 based on patterns, collapses, swaps and vertex relocations.*/
int MMG2D_adptri(MMG5_pMesh mesh,MMG5_pSol met) {

  int it,nns,ns,nnc,nc,nnsw,nsw,nnm,nm;

  nns = nnc = nnsw = nnm = it = 0;
  ns=0; nc=0; nsw=0; nm=0;
  int typchk;
  int color;
  int ret, i;
  int maxit;
  int8_t improve;
    
  int SIZE=1;
  int vect_ns[SIZE];
  int vect_nc[SIZE];
  int vect_nsw[SIZE];
  int vect_nm[SIZE];
  
  starpu_data_handle_t vector_mesh,vector_met, vector_ns, vector_nc,vector_nsw,vector_nm;
  
  
  for (i=0; i< SIZE; i++)
      {
        vect_ns[0]=ns;
        vect_nc[0]=nc;
        vect_nsw[0]=nsw;
        vect_nm[0]=nm;

       }
  
  
  starpu_vector_data_register(&vector_mesh, STARPU_MAIN_RAM, (uintptr_t)mesh, 1, sizeof(MMG5_pMesh)); 
  starpu_vector_data_register(&vector_met, STARPU_MAIN_RAM, (uintptr_t)met, 1, sizeof(MMG5_pSol));
  starpu_vector_data_register(&vector_ns, STARPU_MAIN_RAM, (uintptr_t)vect_ns, SIZE, sizeof(vect_ns[0]));
  starpu_vector_data_register(&vector_nc, STARPU_MAIN_RAM, (uintptr_t)vect_nc, SIZE, sizeof(vect_nc[0]));
  starpu_vector_data_register(&vector_nsw, STARPU_MAIN_RAM, (uintptr_t)vect_nsw, SIZE, sizeof(vect_nsw[0]));
  starpu_vector_data_register(&vector_nm, STARPU_MAIN_RAM,(uintptr_t)vect_nm, SIZE, sizeof(vect_nm[0]));
  
   do {

    if ( !mesh->info.noinsert ) {
      //ns = MMG2D_adpspl(mesh,met,color);
      
     fprintf(stdout," vector_ns_2-------------= %d .\n", vect_ns[0]);
     for (color=0; color< 1; color++)
     {
      ret = starpu_task_insert(&adpspl_codelet,
                STARPU_RW, vector_mesh, 
                STARPU_RW, vector_met, 
                STARPU_RW, vector_ns, 
		 STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      }
      
      starpu_task_wait_for_all();
      fprintf(stdout," vector_ns_222222222222222= %d .\n", vect_ns[0]);
      
      ns=vect_ns[0];       
      if ( ns < 0 ) {
        fprintf(stderr,"  ## Problem in function adpspl."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }

      //nc = MMG2D_adpcol(mesh,met,color);
      
     for (color=0; color< 1; color++)
     {
      ret = starpu_task_insert(&adpcol_codelet,
                STARPU_RW, vector_mesh, 
                STARPU_RW, vector_met, 
                 STARPU_RW, vector_nc, 
		 STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      }
      starpu_task_wait_for_all();
      
      fprintf(stdout," vector_nc_2= %d .\n", vect_nc[0]);
      nc=vect_nc[0];  
          
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Problem in function adpcol."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else {
      ns = 0;
      nc = 0;
    }
    fprintf(stdout," helloooooooooooooo  \n");
    if ( !mesh->info.noswap ) {
      //nsw = MMG2D_swpmsh(mesh,met,2,color);
      
    fprintf(stdout," helloooooooooooooo111111&  \n");
    
    typchk=2;
      
    for (color=0; color< 1; color++)
     {

     ret = starpu_task_insert(&swpmsh_codelet,
                STARPU_RW, vector_mesh, 
                STARPU_RW, vector_met, 
                STARPU_RW, vector_nsw, 
                STARPU_VALUE, &typchk, sizeof(typchk), 
		 STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
     STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      }
      fprintf(stdout," helloooooooooooooo12222222221&  \n");
      starpu_task_wait_for_all();
      
      fprintf(stdout," vector_nsw_2= %d .\n", vect_nsw[0]);
      
      nsw=vect_nsw[0];
      
      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Problem in function swpmsh."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
     fprintf(stdout," byeeeeeeee \n");
      nsw = 0;

    if ( !mesh->info.nomove ) {
      //nm = MMG2D_movtri(mesh,met,1,0,color);
      
    maxit=1;
    improve=0;
      
    for (color=0; color< 1; color++)
     {
      ret = starpu_task_insert(&movtri_codelet,
                STARPU_RW, vector_mesh, 
                STARPU_RW, vector_met, 
                STARPU_RW, vector_nm, 
                STARPU_VALUE, &maxit, sizeof(maxit), 
                STARPU_VALUE, &improve, sizeof(improve), 
		 STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
      STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
      
      }
      starpu_task_wait_for_all();
      
      fprintf(stdout," vector_nm_2= %d .\n", vect_nm[0]);
      nm=vect_nm[0];      
      
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Problem in function movtri. "
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nm = 0;

    nns  += ns;
    nnc  += nc;
    nnsw += nsw;
    nnm  += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nsw+nm > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nsw,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && (nc+ns+nsw+nm > 0) );

  /* Last iterations of vertex relocation only */
  if ( !mesh->info.nomove ) {
    maxit=5;
   improve=1;
    //nm = MMG2D_movtri(mesh,met,5,1,color);  
     
   for (color=0; color< 1; color++)
     {
    ret = starpu_task_insert(&movtri_codelet,
                STARPU_RW, vector_mesh, 
                STARPU_RW, vector_met, 
                STARPU_RW, vector_nm, 
                STARPU_VALUE, &maxit, sizeof(maxit), 
                STARPU_VALUE, &improve, sizeof(improve), 
		 STARPU_VALUE, &color, sizeof(color), 
                                 0);
                                 
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit"); 
    
    }
    starpu_task_wait_for_all();
    fprintf(stdout," vector_nm_2= %d .\n", vect_nm[0]);
    nm=vect_nm[0];      
    
    if ( nm < 0 ) {
      fprintf(stderr,"  ## Problem in function movtri. Unable to complete mesh."
              " Exit program.\n");
      return 0;
    }
    nnm += nm;
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved, %d iter. \n",nns,nnc,nnsw,nnm,it);
  }
  
  starpu_data_unregister(vector_mesh);
  starpu_data_unregister(vector_met);
  starpu_data_unregister(vector_ns);
  starpu_data_unregister(vector_nc);
  starpu_data_unregister(vector_nsw);
  starpu_data_unregister(vector_nm);
  
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return -1 if failed or number of new points.
 *
 * Analysis and splitting routine for edges in the final step of the algorithm;
 * edges are only splitted on a one-by-one basis
 *
 */
 
 void MMG2D_starpu_adpspl(void *buffers[], void *cl_arg) {


  int nx_mesh, nx_met, nx_ns;
  struct starpu_vector_interface *vect_mesh, *vect_met, *vect_ns;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  int color;
  int *val_ns;
  //int ns;
 
  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  
  vect_ns = (struct starpu_vector_interface *) buffers[2];
  nx_ns = STARPU_VECTOR_GET_NX(vect_ns);
  val_ns = (int *)STARPU_VECTOR_GET_PTR(vect_ns);
     
  //ns=val_ns;  
  starpu_codelet_unpack_args(cl_arg,&color);

 
  val_ns[0]=MMG2D_adpspl(mesh,met,color);
  
  /*if ( ns < 0 ) {
     fprintf(stderr,"  ## Unable to adpspl mesh. Exiting.\n");
   } */
 
 }
int MMG2D_adpspl(MMG5_pMesh mesh,MMG5_pSol met, int color) {
  MMG5_pTria         pt;
  double             lmax,len;
  int                k,ns,nt,ip,ier;
  int8_t             i,i1,i2,imax;

  ns = 0;

  /*loop until nt to avoid the split of new triangle*/
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;

    imax = -1;
    lmax = -1.0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }

    if ( lmax < MMG2D_LOPTL ) continue;
    else if ( MG_SIN(pt->tag[imax]) ) continue;

    /* Check the feasibility of splitting */
    ip = MMG2D_chkspl(mesh,met,k,imax);

    /* Lack of memory; abort the routine */
    if ( ip < 0 ){
      return ns;
    }
    else if ( ip > 0 ) {

      ier = MMG2D_split1b(mesh,k,imax,ip);

      /* Lack of memory; abort the routine */
      if ( !ier ) {
        MMG2D_delPt(mesh,ip);
        return ns;
      }
      ns += ier;
    }
  }

  return ns;
}

void MMG2D_starpu_adpcol(void *buffers[], void *cl_arg) {

 
  int nx_mesh, nx_met, nx_nc;
  struct starpu_vector_interface *vect_mesh, *vect_met, *vect_nc;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  //int nc;
  int color;
  int *val_nc;
  //int nc;
 
  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  vect_nc = (struct starpu_vector_interface *) buffers[2];
  nx_nc = STARPU_VECTOR_GET_NX(vect_nc);
  val_nc = (int *)STARPU_VECTOR_GET_PTR(vect_nc);
     
  //nc=val_nc[0];
  starpu_codelet_unpack_args(cl_arg,&color);
 
  val_nc[0]=MMG2D_adpspl(mesh,met,color);
  
 }

/* Analysis and collapse routine for edges in the final step of the algorithm */
int MMG2D_adpcol(MMG5_pMesh mesh,MMG5_pSol met, int color) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            len;
  int               k,nc,ilist,list[MMG2D_LONMAX+2];
  int8_t            i,i1,i2,open;
  
  nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;

    /* Check edge length, and attempt collapse */
    pt->flag = 0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;

      open = ( mesh->adja[3*(k-1)+1+i] == 0 ) ? 1 : 0;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > MMG2D_LOPTS ) continue;

      ilist = MMG2D_chkcol(mesh,met,k,i,list,2);

      if ( ilist > 3 || ( ilist==3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  return nc;
}

/* Analyze points to relocate them according to a quality criterion */
void MMG2D_starpu_movtri(void *buffers[], void *cl_arg) {

 
  int nx_mesh, nx_met, nx_nm;
  struct starpu_vector_interface *vect_mesh, *vect_met, *vect_nm;
  
  MMG5_pMesh mesh;
  MMG5_pSol met;
  
  int maxit,improve,color;
  
  int *val_nm;
  //int nm;
 
  vect_mesh = (struct starpu_vector_interface *) buffers[0];
  nx_mesh = STARPU_VECTOR_GET_NX(vect_mesh);
  mesh = (MMG5_pMesh)STARPU_VECTOR_GET_PTR(vect_mesh);
  
  vect_met = (struct starpu_vector_interface *) buffers[1];
  nx_met = STARPU_VECTOR_GET_NX(vect_met);
  met = ( MMG5_pSol)STARPU_VECTOR_GET_PTR(vect_met);
  
  vect_nm = (struct starpu_vector_interface *) buffers[2];
  nx_nm = STARPU_VECTOR_GET_NX(vect_nm);
  val_nm = (int *)STARPU_VECTOR_GET_PTR(vect_nm);
  //nm=val_nm[0];
     
     
  starpu_codelet_unpack_args(cl_arg, &maxit, &improve,&color);
 
  val_nm[0]=MMG2D_movtri(mesh,met,maxit,improve,color);
  
}

int MMG2D_movtri(MMG5_pMesh mesh,MMG5_pSol met,int maxit,int8_t improve, int color) {
  MMG5_pTria           pt;
  MMG5_pPoint          p0;
  int                  base,k,nnm,nm,ns,it,ilist,list[MMG2D_LONMAX+2];
  int8_t               i,ier;

  it = nnm = 0;
  base = 0;

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MMG2D_EOK(pt,color) || pt->ref < 0 ) continue;

      for (i=0; i<3; i++) {
        p0 = &mesh->point[pt->v[i]];
        if ( p0->flag == base || MG_SIN(p0->tag) || p0->tag & MG_NOM ) continue;

        ilist = MMG2D_boulet(mesh,k,i,list);

        if ( MG_EDG(p0->tag) ) {
          ier = MMG2D_movedgpt(mesh,met,ilist,list,improve);
          if ( ier ) ns++;
        }
        else {
          if ( met->size == 3 && met->m )
            ier = MMG2D_movintpt_ani(mesh,met,ilist,list,improve);
          else
            ier = MMG2D_movintpt(mesh,met,ilist,list,improve);
        }

        if ( ier ) {
          nm++;
          p0->flag = base;
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
  }
  while ( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm > 0 )
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

  return nnm;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 1 if success, 0 if strongly fail.
 *
 * Mesh adaptation -- new version of mmg2d1.c
 *
 **/
int MMG2D_mmg2d1n(MMG5_pMesh mesh,MMG5_pSol met) {

  /* Stage 0: mesh coloration with metis*/
   
  fprintf(stdout,"  -- Call Metis for partioning  \n");

  idx_t *part;
 
  int color;

  /** Allocate the table part */
  MMG5_SAFE_CALLOC(part,mesh->nt,idx_t,return 0);
  int status, i;


  int nelt= mesh->nt;

  MMG5_pTria pt;

  status=MMG_part_meshElts2metis( mesh, part, (idx_t)mesh->ncolors );

  for (i=0; i< mesh->nt; i++)
  {

	  pt= &mesh->tria[i+1];
          pt->color = part[i];
          pt->ref = pt->color;

          //fprintf(stdout,"--elt=%d, --color= %d \n",i, pt->color);


  }

  fprintf(stdout,"  --END Call Metis \n");

   /*save result*/
  if ( MMG2D_saveMesh(mesh, "metis.mesh") != 1 )
    exit(EXIT_FAILURE);
   
 /*************************************/
 
  /* Stage 1: creation of a geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");
  
  if ( !MMG2D_anatri(mesh,met,1) ) {
    fprintf(stderr,"  ## Unable to split mesh-> Exiting.\n");
    return 0;
  }

  /* Stage 2: creation of a computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  if ( !MMG2D_defsiz(mesh,met) ) {
    fprintf(stderr,"  ## Metric undefined. Exit program.\n");
    return 0;
  }

  MMG5_gradation_info(mesh);
  if ( mesh->info.hgrad > 0. ) {
    if (!MMG2D_gradsiz(mesh,met) ) {
      fprintf(stderr,"  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }

  if ( mesh->info.hgradreq > 0. ) {
    MMG2D_gradsizreq(mesh,met);
  }

  if ( !MMG2D_anatri(mesh,met,2)  ) {
    fprintf(stderr,"  ## Unable to proceed adaptation. Exit program.\n");
    return 0;
  }

  /* Stage 3: fine mesh improvements */
  if ( !MMG2D_adptri(mesh,met) ) {
    fprintf(stderr,"  ## Unable to make fine improvements. Exit program.\n");
    return 0;
  }


  return 1;
}
