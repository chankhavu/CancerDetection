/*************************************************************************
 *  -
 *
 * $Id: apprx-diam.h,v 1.1 2002/07/24 12:43:44 greg Exp $
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

#ifndef _apprx_diam_h_
#define _apprx_diam_h_

#include <exact-diam.h>

extern void _VerboseInApprxDiam();
extern void _NoVerboseInApprxDiam();

extern double _EstimeDiameterInOneList( typeSegment *theDiam,
					double **theList,
					const int first,
					const int last,
					const int dim,
					double _epsilon_ );

#endif
