/*************************************************************************
 *  -
 *
 * $Id: peled-like.h,v 1.1 2002/07/24 12:43:44 greg Exp $
 *
 * Copyright INRIA
 *
 * Author:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Sat Jun  9 15:22:00 MEST 2001
 *
 *
 * ADDITIONS, CHANGES
 * 
 *
 */

#ifndef _peled_like_h_
#define _peled_like_h_

#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <alloc.h>
#include <rand.h>
#include <util.h>

extern void _VerboseInPeledLike();
extern void _NoVerboseInPeledLike();


extern void set_evolution_peled_min_nb_points_to_grow();
extern void set_evolution_peled_min_nb_points_as_it_is();

extern void set_init_peled_diameter_to_largest_dim();
extern void set_init_peled_diameter_to_iterated_maximal_seg();
extern void set_init_peled_diameter_to_maximal_seg();

extern void set_update_peled_diameter_to_quadratic();
extern void set_update_peled_diameter_to_smart_max_exact();
extern void set_update_peled_diameter_to_max_exact();
extern void set_update_peled_diameter_to_max_exact_with_diameter();

extern void adapt_small_box_splitting();
extern void allow_to_split_small_box();
extern void do_not_allow_to_split_small_box();

extern void print_peled_env();


extern double _Peled_LikeDiameterInOneList( typeSegment *theDiam,
				       double **theList,
				       const int first,
				       const int last,
				       const int dim );


#endif
