// ============================================================================
// Copyright (c) 2013, FORTH-ICS / CARV 
//                     (Foundation for Research & Technology -- Hellas,
//                      Institute of Computer Science,
//                      Computer Architecture & VLSI Systems Laboratory)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 
// ============================================================================
// The code in this file is taken from the NAS benchmarks, modified to 
// run on the Formic architecture. The NAS benchmarks are available under
// the following copyright:
//
//  This benchmark is part of the NAS Parallel Benchmark 3.3 suite.
//  It is described in NAS Technical Report 95-020.
//
//  Permission to use, copy, distribute and modify this software
//  for any purpose with or without fee is hereby granted.  We
//  request, however, that all derived work reference the NAS
//  Parallel Benchmarks 3.3. This software is provided "as is"
//  without express or implied warranty.
// ============================================================================
// Code modified by Iakovos Mavroidis
// ============================================================================


#ifdef ARCH_MB
#include <fmpi.h>
#include <kernel_toolset.h>
#include <arch.h>
  
#define strlen kt_strlen
#define strcpy kt_strcpy
#define strdup kt_strdup
#define strstr kt_strstr
#define strcmp kt_strcmp
#define printf kt_printf
#define sprintf kt_sprintf
#define malloc kt_malloc
#define MPI_Wtime FMPI_Wtime
#define free   kt_free
#define memcpy kt_memcpy
  
#else

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

#endif
#define TYPE 2

#define CLASS 'S'
#define NUM_SAMPLES 1728
#define STD_DEVIATION 128
#define NUM_SOURCES 4

//#define CLASS 'W'
//#define NUM_SAMPLES 13824
//#define STD_DEVIATION 256
//#define NUM_SOURCES 8

//#define CLASS 'A'
//#define NUM_SAMPLES 110592
//#define STD_DEVIATION 512
//#define NUM_SOURCES 16

#define COMPILETIME "16 May 2012"
#define NPBVERSION "3.3.1"
#define MPICC "cc"
#define CFLAGS "-O3 "
#define CLINK "cc"
#define CLINKFLAGS "(none)"
#define CMPI_LIB "-lmpi"
#define CMPI_INC "(none)"


#ifndef CLASS
#define CLASS 'S'
#endif

#ifndef ARCH_MB
char *kt_strstr(const char *s1, const char *s2)
{
    char *cp = (char *) s1;
    char *str1, *str2;
  
    if (!*s2) return (char *) s1;
  
    while (*cp)
    {
      str1 = cp;
      str2 = (char *) s2;
  
      while (*str1 && *str2 && !(*str1 - *str2)) str1++, str2++;
      if (!*str2) return cp;
      cp++;
    }
  
    return NULL;
}

char *strdup_(const char *s)
{
  char *t;
  int len;

  if (!s) return NULL;
  len = strlen(s);
  t = (char *) malloc(len + 1);
  memcpy(t, s, len + 1);
  return t;
}

#define strdup strdup_

#endif


/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/


static void c_print_results( char   *name,
                      char   class,
                      int    n1, 
                      int    n2,
                      int    n3,
                      int    niter,
                      int    nprocs_compiled,
                      int    nprocs_total,
                      float t,
                      int mops,
		      char   *optype,
                      int    passed_verification,
                      char   *npbversion,
                      char   *compiletime,
                      char   *mpicc,
                      char   *clink,
                      char   *cmpi_lib,
                      char   *cmpi_inc,
                      char   *cflags,
                      char   *clinkflags )
{
//    char *evalue="1000";

    printf( "\r\n\r\n %s Benchmark Completed\r\n", name ); 

    printf( " Class           =                        %c\r\n", class );

    if( n3 == 0 ) {
        long nn = n1;
        if ( n2 != 0 ) nn *= n2;
        printf( " Size            =             %12ld\r\n", nn );   /* as in IS */
    }
    else
        printf( " Size            =              %3dx %3dx %3d\r\n", n1,n2,n3 );

    printf( " Iterations      =             %12d\r\n", niter );
 
//    printf( " Time in seconds =             %12.2f\r\n", t );
    printf( " Time in cycles =             %12d\r\n", (int)t );

    printf( " Total processes =             %12d\r\n", nprocs_total );

    if ( nprocs_compiled != 0 )
        printf( " Compiled procs  =             %12d\r\n", nprocs_compiled );

//    printf( " Mop/s total     =             %12.2f\r\n", mops );
    printf( " Op/Kcycles total     =             %12d\r\n", (int)mops );

//    printf( " Mop/s/process   =             %12.2f\r\n", mops/((float) nprocs_total) );
    printf( " Op/Kcycles/process   =             %12d/%12d\r\n", (int)mops, (int) nprocs_total);

    printf( " Operation type  = %24s\r\n", optype);

    if( passed_verification )
        printf( " Verification    =               SUCCESSFUL\r\n" );
    else
        printf( " Verification    =             UNSUCCESSFUL\r\n" );

    printf( " Version         =             %12s\r\n", npbversion );

    printf( " Compile date    =             %12s\r\n", compiletime );

    printf( "\r\n Compile options:\r\n" );

    printf( "    MPICC        = %s\r\n", mpicc );

    printf( "    CLINK        = %s\r\n", clink );

    printf( "    CMPI_LIB     = %s\r\n", cmpi_lib );

    printf( "    CMPI_INC     = %s\r\n", cmpi_inc );

    printf( "    CFLAGS       = %s\r\n", cflags );

    printf( "    CLINKFLAGS   = %s\r\n", clinkflags );
#ifdef SMP
    evalue = getenv("MP_SET_NUMTHREADS");
    printf( "   MULTICPUS = %s\r\n", evalue );
#endif

    printf( "\r\n\r\n" );
    printf( " Please send feedbacks and/or the results of this run to:\r\n\r\n" );
    printf( " NPB Development Team\r\n" );
    printf( " npb@nas.nasa.gov\r\n\r\n\r\n" );
/*    printf( " Please send the results of this run to:\n\n" );
    printf( " NPB Development Team\n" );
    printf( " Internet: npb@nas.nasa.gov\n \n" );
    printf( " If email is not available, send this to:\n\n" );
    printf( " MS T27A-1\n" );
    printf( " NASA Ames Research Center\n" );
    printf( " Moffett Field, CA  94035-1000\n\n" );
    printf( " Fax: 650-604-3957\n\n" );*/
}
 

int verify(char *bmname,float rnm2){
    float verify_value=0.0;
//    float epsilon=1.0E-8;
    float epsilon=1000;
    char cls=CLASS;
    int verified=-1;
    if (cls != 'U') {
       if(cls=='S') {
         if(strstr(bmname,"BH")){
//           verify_value=30892725.0;
           verify_value=30892744.0;
         }else if(strstr(bmname,"WH")){
           verify_value=67349758.0;
         }else if(strstr(bmname,"SH")){
           verify_value=58875767.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
         }
         verified = 0;
       }else if(cls=='W') {
         if(strstr(bmname,"BH")){
  	   verify_value = 4102461.0;
  	   verify_value = 4104716.0;
         }else if(strstr(bmname,"WH")){
  	   verify_value = 204280762.0;
         }else if(strstr(bmname,"SH")){
  	   verify_value = 186944764.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
         }
         verified = 0;
       }else if(cls=='A') {
         if(strstr(bmname,"BH")){
  	   verify_value = 17809491.0;
  	   verify_value = 17808000.0;
         }else if(strstr(bmname,"WH")){
  	   verify_value = 1289925229.0;
         }else if(strstr(bmname,"SH")){
  	   verify_value = 610856482.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
         }
  	 verified = 0;
       }else if(cls=='B') {
         if(strstr(bmname,"BH")){
  	   verify_value = 4317114.0;
         }else if(strstr(bmname,"WH")){
  	   verify_value = 7877279917.0;
         }else if(strstr(bmname,"SH")){
  	   verify_value = 1836863082.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
  	   verified = 0;
         }
       }else if(cls=='C') {
         if(strstr(bmname,"BH")){
  	   verify_value = 0.0;
         }else if(strstr(bmname,"WH")){
  	   verify_value = 0.0;
         }else if(strstr(bmname,"SH")){
  	   verify_value = 0.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
  	   verified = -1;
         }
       }else if(cls=='D') {
         if(strstr(bmname,"BH")){
  	   verify_value = 0.0;
         }else if(strstr(bmname,"WH")){
  	   verify_value = 0.0;
         }else if(strstr(bmname,"SH")){
  	   verify_value = 0.0;
         }else{
           printf("No such benchmark as %s.\n",bmname);
         }
         verified = -1;
       }else{
         printf("No such class as %c.\n",cls);
       }
       printf(" %s L2 Norm = %f\n",bmname,rnm2);
       if(verified==-1){
  	 printf(" No verification was performed.\n");
       }else if( rnm2 - verify_value < epsilon &&
                 rnm2 - verify_value > -epsilon) {  /* abs here does not work on ALTIX */
  	  verified = 1;
  	  printf(" Deviation = %f\n",(rnm2 - verify_value));
       }else{
  	 verified = 0;
  	 printf(" The correct verification value = %f\n",verify_value);
  	 printf(" Got value = %f\n",rnm2);
       }
    }else{
       verified = -1;
    }
    return  verified;  
  }

int ipowMod(int a,long long int n,int md){ 
  int seed=1,q=a,r=1;
  if(n<0){
    printf("ipowMod: exponent must be nonnegative exp=%lld\n",n);
    n=-n; /* temp fix */
/*    return 1; */
  }
  if(md<=0){
    printf("ipowMod: module must be positive mod=%d",md);
    return 1;
  }
  if(n==0) return 1;
  while(n>1){
    int n2 = n/2;
    if (n2*2==n){
       seed = (q*q)%md;
       q=seed;
       n = n2;
    }else{
       seed = (r*q)%md;
       r=seed;
       n = n-1;
    }
  }
  seed = (r*q)%md;
  return seed;
}

/////////////// DGraph ///////////

#ifndef _DGRAPH
#define _DGRAPH

#define BLOCK_SIZE  128
#define SMALL_BLOCK_SIZE 32

typedef struct{
  int id;
  void *tail,*head;
  int length,width,attribute,maxWidth;
}DGArc;

typedef struct{
  int maxInDegree,maxOutDegree;
  int inDegree,outDegree;
  int id;
  char *name;
  DGArc **inArc,**outArc;
  int depth,height,width;
  int color,attribute,address,verified;
  void *feat;
}DGNode;

typedef struct{
  int maxNodes,maxArcs;
  int id;
  char *name;
  int numNodes,numArcs;
  DGNode **node;
  DGArc **arc;
} DGraph;

DGArc *newArc(DGNode *tl,DGNode *hd);
void arcShow(DGArc *ar);
DGNode *newNode(char *nm);
void nodeShow(DGNode* nd);

DGraph* newDGraph(char *nm);
int AttachNode(DGraph *dg,DGNode *nd);
int AttachArc(DGraph *dg,DGArc* nar);
void graphShow(DGraph *dg,int DetailsLevel);

#endif

DGArc *newArc(DGNode *tl,DGNode *hd){
  DGArc *ar=(DGArc *)malloc(sizeof(DGArc));
  ar->tail=tl;
  ar->head=hd;
  return ar;
}
void arcShow(DGArc *ar){
  DGNode *tl=(DGNode *)ar->tail,
         *hd=(DGNode *)ar->head;
  printf("%d. |%s ->%s\n",ar->id,tl->name,hd->name);
}

DGNode *newNode(char *nm){
  DGNode *nd=(DGNode *)malloc(sizeof(DGNode));
  nd->attribute=0;
  nd->color=0;
  nd->inDegree=0;
  nd->outDegree=0;
  nd->maxInDegree=SMALL_BLOCK_SIZE;
  nd->maxOutDegree=SMALL_BLOCK_SIZE;
  nd->inArc=(DGArc **)malloc(nd->maxInDegree*sizeof(DGArc*));
  nd->outArc=(DGArc **)malloc(nd->maxOutDegree*sizeof(DGArc*));
  nd->name=strdup(nm);
  nd->feat=NULL;
  return nd;
}
void nodeShow(DGNode* nd){
  printf("%3d.%s: (%d,%d)\n",
	           nd->id,nd->name,nd->inDegree,nd->outDegree);
/*
  if(nd->verified==1) fprintf(stderr,"%ld.%s\t: usable.",nd->id,nd->name);
  else if(nd->verified==0)  fprintf(stderr,"%ld.%s\t: unusable.",nd->id,nd->name);
  else  fprintf(stderr,"%ld.%s\t: notverified.",nd->id,nd->name);   
*/
}

DGraph* newDGraph(char* nm){
  DGraph *dg=(DGraph *)malloc(sizeof(DGraph));
  dg->numNodes=0;
  dg->numArcs=0;
  dg->maxNodes=BLOCK_SIZE;
  dg->maxArcs=BLOCK_SIZE;
  dg->node=(DGNode **)malloc(dg->maxNodes*sizeof(DGNode*));
  dg->arc=(DGArc **)malloc(dg->maxArcs*sizeof(DGArc*));
  dg->name=strdup(nm);
  return dg;
}
int AttachNode(DGraph* dg, DGNode* nd) {
  int i=0,j,len=0;
  DGNode **nds =NULL, *tmpnd=NULL;
  DGArc **ar=NULL;

	if (dg->numNodes == dg->maxNodes-1 ) {
	  dg->maxNodes += BLOCK_SIZE;
//          nds =(DGNode **) calloc(dg->maxNodes,sizeof(DGNode*));
          nds =(DGNode **) malloc(dg->maxNodes*sizeof(DGNode*));
	  memcpy(nds,dg->node,(dg->maxNodes-BLOCK_SIZE)*sizeof(DGNode*));
	  free(dg->node);
	  dg->node=nds;
	}

        len = strlen( nd->name);
	for (i = 0; i < dg->numNodes; i++) {
	  tmpnd =dg->node[ i];
	  ar=NULL;
	  if ( strlen( tmpnd->name) != len ) continue;
//	  if ( strncmp( nd->name, tmpnd->name, len) ) continue;
	  if ( strcmp( nd->name, tmpnd->name) ) continue;
	  if ( nd->inDegree > 0 ) {
	    tmpnd->maxInDegree += nd->maxInDegree;
//            ar =(DGArc **) calloc(tmpnd->maxInDegree,sizeof(DGArc*));
            ar =(DGArc **) malloc(tmpnd->maxInDegree*sizeof(DGArc*));
	    memcpy(ar,tmpnd->inArc,(tmpnd->inDegree)*sizeof(DGArc*));
	    free(tmpnd->inArc);
	    tmpnd->inArc=ar;
	    for (j = 0; j < nd->inDegree; j++ ) {
	      nd->inArc[ j]->head = tmpnd;
	    }
	    memcpy( &(tmpnd->inArc[ tmpnd->inDegree]), nd->inArc, nd->inDegree*sizeof( DGArc *));
	    tmpnd->inDegree += nd->inDegree;
	  } 	
	  if ( nd->outDegree > 0 ) {
	    tmpnd->maxOutDegree += nd->maxOutDegree;
//            ar =(DGArc **) calloc(tmpnd->maxOutDegree,sizeof(DGArc*));
            ar =(DGArc **) malloc(tmpnd->maxOutDegree*sizeof(DGArc*));
	    memcpy(ar,tmpnd->outArc,(tmpnd->outDegree)*sizeof(DGArc*));
	    free(tmpnd->outArc);
	    tmpnd->outArc=ar;
	    for (j = 0; j < nd->outDegree; j++ ) {
	      nd->outArc[ j]->tail = tmpnd;
	    }			
	    memcpy( &(tmpnd->outArc[tmpnd->outDegree]),nd->outArc,nd->outDegree*sizeof( DGArc *));
	    tmpnd->outDegree += nd->outDegree;
	  } 
	  free(nd); 
	  return i;
	}
	nd->id = dg->numNodes;
	dg->node[dg->numNodes] = nd;
	dg->numNodes++;
return nd->id;
}
int AttachArc(DGraph *dg,DGArc* nar){
int	arcId = -1;
int i=0,newNumber=0;
DGNode	*head = nar->head,
	*tail = nar->tail; 
DGArc **ars=NULL,*probe=NULL;
/*fprintf(stderr,"AttachArc %ld\n",dg->numArcs); */
	if ( !tail || !head ) return arcId;
	if ( dg->numArcs == dg->maxArcs-1 ) {
	  dg->maxArcs += BLOCK_SIZE;
//          ars =(DGArc **) calloc(dg->maxArcs,sizeof(DGArc*));
          ars =(DGArc **) malloc(dg->maxArcs*sizeof(DGArc*));
	  memcpy(ars,dg->arc,(dg->maxArcs-BLOCK_SIZE)*sizeof(DGArc*));
	  free(dg->arc);
	  dg->arc=ars;
	}
	for(i = 0; i < tail->outDegree; i++ ) { /* parallel arc */
	  probe = tail->outArc[ i];
	  if(probe->head == head
	     &&
	     probe->length == nar->length
            ){
            free(nar);
	    return probe->id;   
	  }
	}
	
	nar->id = dg->numArcs;
	arcId=dg->numArcs;
	dg->arc[dg->numArcs] = nar;
	dg->numArcs++;
	
	head->inArc[ head->inDegree] = nar;
	head->inDegree++;
	if ( head->inDegree >= head->maxInDegree ) {
	  newNumber = head->maxInDegree + SMALL_BLOCK_SIZE;
//          ars =(DGArc **) calloc(newNumber,sizeof(DGArc*));
          ars =(DGArc **) malloc(newNumber*sizeof(DGArc*));
	  memcpy(ars,head->inArc,(head->inDegree)*sizeof(DGArc*));
	  free(head->inArc);
	  head->inArc=ars;
	  head->maxInDegree = newNumber;
	}
	tail->outArc[ tail->outDegree] = nar;
	tail->outDegree++;
	if(tail->outDegree >= tail->maxOutDegree ) {
	  newNumber = tail->maxOutDegree + SMALL_BLOCK_SIZE;
//          ars =(DGArc **) calloc(newNumber,sizeof(DGArc*));
          ars =(DGArc **) malloc(newNumber*sizeof(DGArc*));
	  memcpy(ars,tail->outArc,(tail->outDegree)*sizeof(DGArc*));
	  free(tail->outArc);
	  tail->outArc=ars;
	  tail->maxOutDegree = newNumber;
	}
/*fprintf(stderr,"AttachArc: head->in=%d tail->out=%ld\n",head->inDegree,tail->outDegree);*/
return arcId;
}
void graphShow(DGraph *dg,int DetailsLevel){
  int i=0,j=0;
  printf("%d.%s: (%d,%d)\n",dg->id,dg->name,dg->numNodes,dg->numArcs);
  if ( DetailsLevel < 1) return;
  for (i = 0; i < dg->numNodes; i++ ) {
    DGNode *focusNode = dg->node[ i];
    if(DetailsLevel >= 2) {
      for (j = 0; j < focusNode->inDegree; j++ ) {
	printf("\t ");
	nodeShow(focusNode->inArc[ j]->tail);
      }
    }
    nodeShow(focusNode);
    if ( DetailsLevel < 2) continue;
    for (j = 0; j < focusNode->outDegree; j++ ) {
      printf("\t ");
      nodeShow(focusNode->outArc[ j]->head);
    }	
    printf("---\n");
  }
  printf("----------------------------------------\n");
  if ( DetailsLevel < 3) return;
}

/////////////////////////////////

DGraph *buildSH(char cls){
/*
  Nodes of the graph must be topologically sorted
  to avoid MPI deadlock.
*/
  DGraph *dg;
  int numSources=NUM_SOURCES; /* must be power of 2 */
  int numOfLayers=0,tmpS=numSources>>1;
  int firstLayerNode=0;
  DGArc *ar=NULL;
  DGNode *nd=NULL;
  int mask=0x0,ndid=0,ndoff=0;
  int i=0,j=0;
  char nm[BLOCK_SIZE];
  
  sprintf(nm,"DT_SH.%c",cls);
  dg=newDGraph(nm);

  while(tmpS>1){
    numOfLayers++;
    tmpS>>=1;
  }
  for(i=0;i<numSources;i++){
    sprintf(nm,"Source.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  for(j=0;j<numOfLayers;j++){
    mask=0x00000001<<j;
    for(i=0;i<numSources;i++){
      sprintf(nm,"Comparator.%d",(i+j*firstLayerNode));
      nd=newNode(nm);
      AttachNode(dg,nd);
      ndoff=i&(~mask);
      ndid=firstLayerNode+ndoff;
      ar=newArc(dg->node[ndid],nd);     
      AttachArc(dg,ar);
      ndoff+=mask;
      ndid=firstLayerNode+ndoff;
      ar=newArc(dg->node[ndid],nd);     
      AttachArc(dg,ar);
    }
    firstLayerNode+=numSources;
  }
  mask=0x00000001<<numOfLayers;
  for(i=0;i<numSources;i++){
    sprintf(nm,"Sink.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
    ndoff=i&(~mask);
    ndid=firstLayerNode+ndoff;
    ar=newArc(dg->node[ndid],nd);     
    AttachArc(dg,ar);
    ndoff+=mask;
    ndid=firstLayerNode+ndoff;
    ar=newArc(dg->node[ndid],nd);     
    AttachArc(dg,ar);
  }
return dg;
}
DGraph *buildWH(char cls){
/*
  Nodes of the graph must be topologically sorted
  to avoid MPI deadlock.
*/
  int i=0,j=0;
  int numSources=NUM_SOURCES,maxInDeg=4;
  int numLayerNodes=numSources,firstLayerNode=0;
  int totComparators=0;
  int numPrevLayerNodes=numLayerNodes;
  int id=0,sid=0;
  DGraph *dg;
  DGNode *nd=NULL,*source=NULL,*tmp=NULL,*snd=NULL;
  DGArc *ar=NULL;
  char nm[BLOCK_SIZE];

  sprintf(nm,"DT_WH.%c",cls);
  dg=newDGraph(nm);

  for(i=0;i<numSources;i++){
    sprintf(nm,"Sink.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  totComparators=0;
  numPrevLayerNodes=numLayerNodes;
  while(numLayerNodes>maxInDeg){
    numLayerNodes=numLayerNodes/maxInDeg;
    if(numLayerNodes*maxInDeg<numPrevLayerNodes)numLayerNodes++;
    for(i=0;i<numLayerNodes;i++){
      sprintf(nm,"Comparator.%d",totComparators);
      totComparators++;
      nd=newNode(nm);
      id=AttachNode(dg,nd);
      for(j=0;j<maxInDeg;j++){
        sid=i*maxInDeg+j;
	if(sid>=numPrevLayerNodes) break;
        snd=dg->node[firstLayerNode+sid];
        ar=newArc(dg->node[id],snd);
        AttachArc(dg,ar);
      }
    }
    firstLayerNode+=numPrevLayerNodes;
    numPrevLayerNodes=numLayerNodes;
  }
  source=newNode("Source");
  AttachNode(dg,source);   
  for(i=0;i<numPrevLayerNodes;i++){
    nd=dg->node[firstLayerNode+i];
    ar=newArc(source,nd);
    AttachArc(dg,ar);
  }

  for(i=0;i<dg->numNodes/2;i++){  /* Topological sorting */
    tmp=dg->node[i];
    dg->node[i]=dg->node[dg->numNodes-1-i];
    dg->node[i]->id=i;
    dg->node[dg->numNodes-1-i]=tmp;
    dg->node[dg->numNodes-1-i]->id=dg->numNodes-1-i;
  }
return dg;
}
DGraph *buildBH(char cls){
/*
  Nodes of the graph must be topologically sorted
  to avoid MPI deadlock.
*/
  int i=0,j=0;
  int numSources=NUM_SOURCES,maxInDeg=4;
  int numLayerNodes=numSources,firstLayerNode=0;
  DGraph *dg;
  DGNode *nd=NULL, *snd=NULL, *sink=NULL;
  DGArc *ar=NULL;
  int totComparators=0;
  int numPrevLayerNodes=numLayerNodes;
  int id=0, sid=0;
  char nm[BLOCK_SIZE];

  sprintf(nm,"DT_BH.%c",cls);
  dg=newDGraph(nm);

  for(i=0;i<numSources;i++){
    sprintf(nm,"Source.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  while(numLayerNodes>maxInDeg){
    numLayerNodes=numLayerNodes/maxInDeg;
    if(numLayerNodes*maxInDeg<numPrevLayerNodes)numLayerNodes++;
    for(i=0;i<numLayerNodes;i++){
      sprintf(nm,"Comparator.%d",totComparators);
      totComparators++;
      nd=newNode(nm);
      id=AttachNode(dg,nd);
      for(j=0;j<maxInDeg;j++){
        sid=i*maxInDeg+j;
	if(sid>=numPrevLayerNodes) break;
        snd=dg->node[firstLayerNode+sid];
        ar=newArc(snd,dg->node[id]);
        AttachArc(dg,ar);
      }
    }
    firstLayerNode+=numPrevLayerNodes;
    numPrevLayerNodes=numLayerNodes;
  }
  sink=newNode("Sink");
  AttachNode(dg,sink);   
  for(i=0;i<numPrevLayerNodes;i++){
    nd=dg->node[firstLayerNode+i];
    ar=newArc(nd,sink);
    AttachArc(dg,ar);
  }
return dg;
}

typedef struct{
  int len;
  float* val;
} Arr;
Arr *newArr(int len){
  Arr *arr=(Arr *)malloc(sizeof(Arr));
  arr->len=len;
  arr->val=(float *)malloc(len*sizeof(float));
  return arr;
}
void arrShow(Arr* a){
  if(!a) printf("-- NULL array\n");
  else{
    printf("-- length=%d\n",a->len);
  }
}
float CheckVal(Arr *feat){
  float csum=0.0;
  int i=0;
  for(i=0;i<feat->len;i++){
    csum+=feat->val[i]*feat->val[i]/feat->len; /* The truncation does not work since 
                                                  result will be 0 for large len  */
  }
   return csum;
}
int GetFNumDPar(int* mean, int* stdev){
  *mean=NUM_SAMPLES;
  *stdev=STD_DEVIATION;
  return 0;
}

#ifndef ARCH_MB

//---------------------------------------------------------------------
//   This function is C verson of random number generator randdp.f 
//---------------------------------------------------------------------

double	randlc(X, A)
double *X;
double *A;
{
      static int        KS=0;
      static double	R23, R46, T23, T46;
      double		T1, T2, T3, T4;
      double		A1;
      double		A2;
      double		X1;
      double		X2;
      double		Z;
      int     		i, j;

      if (KS == 0) 
      {
        R23 = 1.0;
        R46 = 1.0;
        T23 = 1.0;
        T46 = 1.0;
    
        for (i=1; i<=23; i++)
        {
          R23 = 0.50 * R23;
          T23 = 2.0 * T23;
        }
        for (i=1; i<=46; i++)
        {
          R46 = 0.50 * R46;
          T46 = 2.0 * T46;
        }
        KS = 1;
      }

/*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */

      T1 = R23 * *A;
      j  = T1;
      A1 = j;
      A2 = *A - T23 * A1;

/*  Break X into two parts such that X = 2^23 * X1 + X2, compute
    Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    X = 2^23 * Z + A2 * X2  (mod 2^46).                            */

      T1 = R23 * *X;
      j  = T1;
      X1 = j;
      X2 = *X - T23 * X1;
      T1 = A1 * X2 + A2 * X1;
      
      j  = R23 * T1;
      T2 = j;
      Z = T1 - T23 * T2;
      T3 = T23 * Z + A2 * X2;
      j  = R46 * T3;
      T4 = j;
      *X = T3 - T46 * T4;
      return(R46 * *X);
} 

#endif

int GetFeatureNum(char *mbname,int id){
#ifndef ARCH_MB
  double tran=314159265.0;
  double A=2*id+1;
  double tmp=randlc(&tran,&A);
#endif

  float denom;
  switch (id) {
        case 0: denom=0.000004464471672; break;
        case 1: denom=0.000013393415017; break;
        case 2: denom=0.000022322358362; break;
        case 3: denom=0.000031251301706; break;
        case 4: denom=0.000040180245051; break;
        case 5: denom=0.000049109188396; break;
        case 6: denom=0.000058038131741; break;
        case 7: denom=0.000066967075085; break;
        case 8: denom=0.000075896018430; break;
        case 9: denom=0.000084824961775; break;
        case 10: denom=0.000093753905119; break;
        case 11: denom=0.000102682848464; break;
        case 12: denom=0.000111611791809; break;
        case 13: denom=0.000120540735153; break;
        case 14: denom=0.000129469678498; break;
        case 15: denom=0.000138398621843; break;
	default: 
#ifdef ARCH_MB
		printf("ERROR!!\n\r"); while(1);
#else
		printf("        case %d: denom=%.15f; break;\n", id, tmp); while(1);
#endif
  }
  char cval='S';
  int mean=NUM_SAMPLES,stdev=128;
  int rtfs=0,len=0;
  GetFNumDPar(&mean,&stdev);
  rtfs=ipowMod((int)(1/denom)*(int)cval,(long long int) (2*id+1),2*stdev);
  if(rtfs<0) rtfs=-rtfs;
  len=mean-stdev+rtfs;
  return len;
}
Arr* RandomFeatures(char *bmname,int fdim,int id){
  int len=GetFeatureNum(bmname,id)*fdim;
  Arr* feat=newArr(len);
  int nxg=2,nyg=2,nzg=2,nfg=5;
  int nx=421,ny=419,nz=1427,nf=3527;
  long long int expon=(len*(id+1))%3141592;
  int seedx=ipowMod(nxg,expon,nx),
      seedy=ipowMod(nyg,expon,ny),
      seedz=ipowMod(nzg,expon,nz),
      seedf=ipowMod(nfg,expon,nf);
  int i=0;
//  if(timer_on){
//    timer_clear(id+1);
//    timer_start(id+1);
//  }
  for(i=0;i<len;i+=fdim){
    seedx=(seedx*nxg)%nx;
    seedy=(seedy*nyg)%ny;
    seedz=(seedz*nzg)%nz;
    seedf=(seedf*nfg)%nf;
    feat->val[i]=seedx;
    feat->val[i+1]=seedy;
    feat->val[i+2]=seedz;
    feat->val[i+3]=seedf;
  }
//  if(timer_on){
//    timer_stop(id+1);
//    fprintf(stderr,"** RandomFeatures time in node %d = %f\n",id,timer_read(id+1));
//  }
  return feat;   
}
void Resample(Arr *a,int blen){
//    long long int i=0,j=0,jlo=0,jhi=0;
    long int i=0,j=0,jlo=0,jhi=0;
    float avval=0.0;
    float *nval=(float *)malloc(blen*sizeof(float));
//    Arr *tmp=newArr(10);
    for(i=0;i<blen;i++) nval[i]=0.0;
    for(i=1;i<a->len-1;i++){
      jlo=(int)(0.5*(2*i-1)*(blen/a->len)); 
      jhi=(int)(0.5*(2*i+1)*(blen/a->len));

      avval=a->val[i]/(jhi-jlo+1);    
      for(j=jlo;j<=jhi;j++){
        nval[j]+=avval;
      }
    }
    nval[0]=a->val[0];
    nval[blen-1]=a->val[a->len-1];
    free(a->val);
    a->val=nval;
    a->len=blen;
}
#define fielddim 4
Arr* WindowFilter(Arr *a, Arr* b,int w){
  int i=0,j=0,k=0;
  float rms0=0.0,rms1=0.0,rmsm1=0.0;
  float weight=((float) (w+1))/(w+2);
 
  w+=1;
//  if(timer_on){
//    timer_clear(w);
//    timer_start(w);
//  }
  if(a->len<b->len) Resample(a,b->len);
  if(a->len>b->len) Resample(b,a->len);
  for(i=fielddim;i<a->len-fielddim;i+=fielddim){
    rms0=(a->val[i]-b->val[i])*(a->val[i]-b->val[i])
	+(a->val[i+1]-b->val[i+1])*(a->val[i+1]-b->val[i+1])
	+(a->val[i+2]-b->val[i+2])*(a->val[i+2]-b->val[i+2])
	+(a->val[i+3]-b->val[i+3])*(a->val[i+3]-b->val[i+3]);
    j=i+fielddim;
    rms1=(a->val[j]-b->val[j])*(a->val[j]-b->val[j])
    	+(a->val[j+1]-b->val[j+1])*(a->val[j+1]-b->val[j+1])
    	+(a->val[j+2]-b->val[j+2])*(a->val[j+2]-b->val[j+2])
    	+(a->val[j+3]-b->val[j+3])*(a->val[j+3]-b->val[j+3]);
    j=i-fielddim;
    rmsm1=(a->val[j]-b->val[j])*(a->val[j]-b->val[j])
	 +(a->val[j+1]-b->val[j+1])*(a->val[j+1]-b->val[j+1])
	 +(a->val[j+2]-b->val[j+2])*(a->val[j+2]-b->val[j+2])
	 +(a->val[j+3]-b->val[j+3])*(a->val[j+3]-b->val[j+3]);
    k=0;
    if(rms1<rms0){
      k=1;
      rms0=rms1;
    }
    if(rmsm1<rms0) k=-1;
    if(k==0){
      j=i+fielddim;
      a->val[i]=weight*b->val[i];
      a->val[i+1]=weight*b->val[i+1];
      a->val[i+2]=weight*b->val[i+2];
      a->val[i+3]=weight*b->val[i+3];  
    }else if(k==1){
      j=i+fielddim;
      a->val[i]=weight*b->val[j];
      a->val[i+1]=weight*b->val[j+1];
      a->val[i+2]=weight*b->val[j+2];
      a->val[i+3]=weight*b->val[j+3];  
    }else { /*if(k==-1)*/
      j=i-fielddim;
      a->val[i]=weight*b->val[j];
      a->val[i+1]=weight*b->val[j+1];
      a->val[i+2]=weight*b->val[j+2];
      a->val[i+3]=weight*b->val[j+3];  
    }	   
  }
//  if(timer_on){
//    timer_stop(w);
//    fprintf(stderr,"** WindowFilter time in node %d = %f\n",(w-1),timer_read(w));
//  }
  return a;
}

int SendResults(DGraph *dg,DGNode *nd,Arr *feat){
  int *tmp;
  int i=0,tag=0;
  tmp = malloc(sizeof(int));
  DGArc *ar=NULL;
  DGNode *head=NULL;
  if(!feat) return 0;
  for(i=0;i<nd->outDegree;i++){
    ar=nd->outArc[i];
    if(ar->tail!=nd) continue;
    head=ar->head;
    tag=ar->id;
    if(head->address!=nd->address){
      *tmp = feat->len;
      MPI_Send(&feat->len,1,MPI_INT,head->address,tag,MPI_COMM_WORLD);
      MPI_Send(feat->val,feat->len,MPI_FLOAT,head->address,tag,MPI_COMM_WORLD);
    }
  }
  free(tmp);
  return 1;
}
Arr* CombineStreams(DGraph *dg,DGNode *nd){
  Arr *resfeat=newArr(NUM_SAMPLES*fielddim);
  int i=0,tag=0;
  int *len;
  DGArc *ar=NULL;
  DGNode *tail=NULL;
  MPI_Status status;
  Arr *feat=NULL,*featp=NULL;

  len = malloc(sizeof(int));
  *len = 0;
  if(nd->inDegree==0) return NULL;
  for(i=0;i<nd->inDegree;i++){
    ar=nd->inArc[i];
    if(ar->head!=nd) continue;
    tail=ar->tail;
    if(tail->address!=nd->address){
      *len=0;
      tag=ar->id;
      MPI_Recv(len,1,MPI_INT,tail->address,tag,MPI_COMM_WORLD,&status);
      feat=newArr(*len);
      MPI_Recv(feat->val,feat->len,MPI_FLOAT,tail->address,tag,MPI_COMM_WORLD,&status);
      resfeat=WindowFilter(resfeat,feat,nd->id);
      free(feat);
    }else{
      featp=(Arr *)tail->feat;
      feat=newArr(featp->len);
      memcpy(feat->val,featp->val,featp->len*sizeof(float));
      resfeat=WindowFilter(resfeat,feat,nd->id);  
      free(feat);
    }
  }
  for(i=0;i<resfeat->len;i++) resfeat->val[i]=((int)resfeat->val[i])/nd->inDegree;
  nd->feat=resfeat;
  free(len);
  return nd->feat;
}
float Reduce(Arr *a,int w){
  float retv=0.0;
//  if(timer_on){
//    timer_clear(w);
//    timer_start(w);
//  }
  retv=(int)(w*CheckVal(a));/* The casting needed for node  
                               and array dependent verifcation */
//  if(timer_on){
//    timer_stop(w);
//    fprintf(stderr,"** Reduce time in node %d = %f\n",(w-1),timer_read(w));
//  }
  return retv;
}

float ReduceStreams(DGraph *dg,DGNode *nd){
  float csum=0.0;
  int i=0,tag=0;
  int *len;
  DGArc *ar=NULL;
  DGNode *tail=NULL;
  Arr *feat=NULL;
  float retv=0.0;

  len = malloc(sizeof(int));
  *len = 0;
  for(i=0;i<nd->inDegree;i++){
    ar=nd->inArc[i];
    if(ar->head!=nd) continue;
    tail=ar->tail;
    if(tail->address!=nd->address){
      MPI_Status status;
      *len=0;
      tag=ar->id;
      MPI_Recv(len,1,MPI_INT,tail->address,tag,MPI_COMM_WORLD,&status);
      feat=newArr(*len);
      MPI_Recv(feat->val,feat->len,MPI_FLOAT,tail->address,tag,MPI_COMM_WORLD,&status);
      csum+=Reduce(feat,(nd->id+1));  
      free(feat);
    }else{
      csum+=Reduce(tail->feat,(nd->id+1));  
    }
  }
//  if(nd->inDegree>0)csum=(((long long int)csum)/nd->inDegree);
  if(nd->inDegree>0)csum=(((long int)csum)/nd->inDegree);
  retv=(nd->id+1)*csum;
  free(len);
  return retv;
}

int ProcessNodes(DGraph *dg,int me){
  float *chksum;
  Arr *feat=NULL;
  int i=0,verified=0,tag;
  DGNode *nd=NULL;
  float *rchksum;
  MPI_Status status;

  chksum = malloc(sizeof(float));
  rchksum = malloc(sizeof(float));
  *chksum = 0.0;
  *rchksum = 0.0;
  for(i=0;i<dg->numNodes;i++){
    nd=dg->node[i];
    if(nd->address!=me) continue;
    if(strstr(nd->name,"Source")){
      nd->feat=RandomFeatures(dg->name,fielddim,nd->id); 
      SendResults(dg,nd,nd->feat);
    }else if(strstr(nd->name,"Sink")){
      *chksum=ReduceStreams(dg,nd);
      tag=dg->numArcs+nd->id; /* make these to avoid clash with arc tags */
      MPI_Send(chksum,1,MPI_FLOAT,0,tag,MPI_COMM_WORLD);
    }else{
      feat=CombineStreams(dg,nd);
      SendResults(dg,nd,feat);
    }
  }
  if(me==0){ /* Report node */
    *rchksum=0.0;
    *chksum=0.0;
    for(i=0;i<dg->numNodes;i++){
      nd=dg->node[i];
      if(!strstr(nd->name,"Sink")) continue;
       tag=dg->numArcs+nd->id; /* make these to avoid clash with arc tags */
      MPI_Recv(rchksum,1,MPI_FLOAT,nd->address,tag,MPI_COMM_WORLD,&status);
      *chksum+=*rchksum;
    }
    verified=verify(dg->name,*chksum);
  }
free(chksum);
free(rchksum);
return verified;
}

#ifdef ARCH_MB
int nas_dt_mpi() {
#else
int main(int argc,char **argv ){
#endif
  int my_rank,comm_size;
  int i;
  DGraph *dg=NULL;
  int verified=0, featnum=0;
  float bytes_sent=2.0,tot_time=0.0;
  unsigned int start = 0;

    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

/*
     if(argc!=2||
                (  strncmp(argv[1],"BH",2)!=0
                 &&strncmp(argv[1],"WH",2)!=0
                 &&strncmp(argv[1],"SH",2)!=0
                )
      ){
      if(my_rank==0){
        printf("** Usage: mpirun -np N ../bin/dt.S GraphName\n");
        printf("** Where \n   - N is integer number of MPI processes\n");
        printf("   - S is the class S, W, or A \n");
        printf("   - GraphName is the communication graph name BH, WH, or SH.\n");
        printf("   - the number of MPI processes N should not be be less than \n");
        printf("     the number of nodes in the graph\n");
      }
      MPI_Finalize();
      exit(0);
    } 
   if(strncmp(argv[1],"BH",2)==0){
      dg=buildBH(CLASS);
    }else if(strncmp(argv[1],"WH",2)==0){
      dg=buildWH(CLASS);
    }else if(strncmp(argv[1],"SH",2)==0){
      dg=buildSH(CLASS);
    }
*/

   if(TYPE==0){
      dg=buildBH(CLASS);
    }else if(TYPE==1){
      dg=buildWH(CLASS);
    }else if(TYPE==2){
      dg=buildSH(CLASS);
    }
//    if(timer_on&&dg->numNodes+1>timers_tot){
//      timer_on=0;
//      if(my_rank==0)
//        fprintf(stderr,"Not enough timers. Node timeing is off. \n");
//    }
    if(dg->numNodes>comm_size){
      if(my_rank==0){
        printf("**  The number of MPI processes should not be less than \n");
        printf("**  the number of nodes in the graph\n");
        printf("**  Number of MPI processes = %d\n",comm_size);
        printf("**  Number nodes in the graph = %d\n",dg->numNodes);
      }
      MPI_Finalize();
//      exit(0);
	 while(1);
    }
    for(i=0;i<dg->numNodes;i++){ 
      dg->node[i]->address=i;
    }
    if( my_rank == 0 ){
      printf( "\n\n NAS Parallel Benchmarks 3.3 -- DT Benchmark\n\n" );
      graphShow(dg,0);
//      timer_clear(0);
//      timer_start(0);
      start = MPI_Wtime();
    }
    verified=ProcessNodes(dg,my_rank);
    
    featnum=NUM_SAMPLES*fielddim;
    bytes_sent=featnum*dg->numArcs;
    bytes_sent/=1048576;
    if(my_rank==0){
      //timer_stop(0);
      //tot_time=timer_read(0);
      tot_time = (float)(MPI_Wtime() - start);
      c_print_results( dg->name,
        	       CLASS,
        	       featnum,
        	       0,
        	       0,
        	       dg->numNodes,
        	       0,
        	       comm_size,
        	       tot_time,
        	       bytes_sent/tot_time,
        	       "bytes transmitted", 
        	       verified,
        	       NPBVERSION,
        	       COMPILETIME,
        	       MPICC,
        	       CLINK,
        	       CMPI_LIB,
        	       CMPI_INC,
        	       CFLAGS,
        	       CLINKFLAGS );
    }          
    MPI_Finalize();
  return 1;
}



