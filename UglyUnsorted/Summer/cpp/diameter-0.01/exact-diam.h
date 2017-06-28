/*************************************************************************
 *  -
 *
 * $Id: exact-diam.h,v 1.1 2002/07/24 12:43:44 greg Exp $
 *
 * Copyright INRIA
 *
 * Author:
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

#ifndef _exact_diam_h_
#define _exact_diam_h_

#include <math.h>
#include <alloc.h>
#include <rand.h>
#include <util.h>

extern void _VerboseInExactDiam();
extern void _NoVerboseInExactDiam();

extern double _ExactDiameterInOneList( typeSegment *theDiam,
				    double **theList,
				    const int first,
				    const int last,
				    const int dim );

extern double _ExactDiameterInTwoLists( typeSegment *theDiam,
				 double **theList1,
				 const int first1,
				 const int last1,
				 double **theList2,
				 const int first2,
				 const int last2,
				 const int dim );

#endif
