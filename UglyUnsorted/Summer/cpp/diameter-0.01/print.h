/*************************************************************************
 *  -
 *
 * $Id: print.h,v 1.1 2002/07/24 12:43:44 greg Exp $
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




#ifndef _print_h_
#define _print_h_

#include <stdio.h>
#include <math.h>

#include <alloc.h>
#include <pick.h>
#include <exact-diam.h>


extern void _PrintEnv( FILE *f,
		const long int init,
		const int nbpoints,
		const double diameterMax,
		const double diameterMin,
		const int psommet,
		enumDistribution typeDistribution,
		char *modelname,
		const int dim );

extern void _PrintSegment( FILE *fp, typeSegment *seg, char *desc, int dim );


extern void _PrintListOfPoints( double **theList,
			     const int first,
			     const int last,
			     const int dim, FILE *f );

extern int _PrintPointIndex( double **theList,
		      const int first,
		      const int last,
		      const double i1,
		      const double i2,
		      const int dim );

#endif
