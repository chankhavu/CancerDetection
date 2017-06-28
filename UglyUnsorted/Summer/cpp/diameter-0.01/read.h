/*************************************************************************
 *  -
 *
 * $Id: read.h,v 1.1 2002/07/24 12:43:44 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon May 15 2000
 *
 *
 * ADDITIONS, CHANGES
 * 
 *
 */




#ifndef _read_h_
#define _read_h_


#include <stdio.h>
#include <string.h>

#include <alloc.h>
#include <rand.h>

extern void _VerboseInRead();
extern void _NoVerboseInRead();

extern void *_ReadPLYModel( char *name, int *nbpoints, int *dim  );
extern void *_ReadPTSModel( char *name, int *nbpoints, int *dim  );
extern void *_ReadModel( char *name, int *nbpoints, int *dim  );

#endif
