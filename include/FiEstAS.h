#ifndef FIESTAS_H
#define FIESTAS_H 1

//----------------------------------------------------------------------

#define DIM unsigned short int
#define NUM unsigned int

#ifndef DOUBLE
#define DOUBLE 1
#endif

#if (DOUBLE==0)

#define DATA float
#define yHUGE 1e+30
#define yTINY 1e-30

#else

#define DATA double
#define yHUGE 1e+300
#define yTINY 1e-300

#endif

//----------------------------------------------------------------------

#define ERROR(x,msg) if(x){ printf("\n ERROR on file %s, line %d: ",__FILE__,__LINE__); printf msg; printf("\n\n\a"); throw(-1); }

#if (VERBOSE>0)
#define yWARNING(x,msg) if(x){ printf(" WARNING: "); printf msg; printf("\n"); }
#define yINFO(msg) printf msg
#else
#define yWARNING(x,msg)
#define yINFO(msg)
#endif

#if (VERBOSE>1)
#define DEBUG(msg) printf msg
#define DEBUG_ERROR(x,msg) ERROR(x,msg)
#else
#define DEBUG(msg)
#define DEBUG_ERROR(x,msg)
#endif

#if (VERBOSE>2)
#define PANIC(stuff) { stuff fflush(stdout); }
#else
#define PANIC(msg)
#endif

//----------------------------------------------------------------------

#include"lib/ySnapshot.h"
#include"lib/yTree.h"
#include"lib/yASCII.h"
#include"lib/utils.h"

#endif
//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
