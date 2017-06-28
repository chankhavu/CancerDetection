/*************************************************************************
 *  -
 *
 * $Id: exact-diam.c,v 1.1 2002/07/24 12:43:44 greg Exp $
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

#include <exact-diam.h>
#include <stdio.h>
#include <print.h>

static int _verbose_ = 0;
void _VerboseInExactDiam()
{
  _verbose_ = 1;
}
void _NoVerboseInExactDiam()
{
  _verbose_ = 0;
}










double _ExactDiameterInOneList( typeSegment *theDiam,
				double **theList,
				const int first,
				const int last,
				const int dim )
{
  int index;

  int f=first;
  int l=last;
  
  int newEstimateIsSmallerThanCurrentEstimate;
  typeSegment theSeg;

  typeListOfSegments theDoubleNormals;

  double newEstimate;

  int newlast;

  
  int verboseWhenReducing     = _GetVerboseWhenReducing();
  int _reduction_mode_in_iterative_ = _GetReductionModeInIterative();
  int tryToReduceQ            = _GetTryToReduceQ();
  int _reduction_mode_of_diameter_ = _GetReductionModeOfDiameter();
  int _reduction_mode_of_dbleNorm_ = _GetReductionModeOfDbleNorm();
  int _Q_scan_                     = _GetQscan();


  int i, j, k, n;

  int suspicion_of_convex_hull = 0;
  int fdn, ldn, idn;
  

  theDoubleNormals.n = 0;
  theDoubleNormals.nalloc = 0;
  theDoubleNormals.seg = NULL;




  theDiam->extremity1 = (double*)NULL;
  theDiam->extremity2 = (double*)NULL;
  theDiam->squareDiameter = 0.0;


  if ( first < 0 || last < 0 ) return( -1.0 );
  if ( first > last ) {
    l = first;
    f = last;
  }
  if ( f == l ) {
    theDiam->extremity1 = theList[f];
    theDiam->extremity2 = theList[l];
    return( 0.0 );
  }

  

  index = _GetRandomIntNb( first, last );


  do {

    /* end conditions
     */    
    newEstimateIsSmallerThanCurrentEstimate = 0;

    /* find a double normal
     */
    newEstimate = _MaximalSegmentInOneList( &theSeg, index, theList,
					    &f, &l, dim );
    
    /* if we get a better estimation
     */
    if ( newEstimate > theDiam->squareDiameter ) {
      
      /* update variables
       */
      *theDiam = theSeg;

      /* keep the maximal segment in list
       */
      if ( _AddSegmentToList( &theSeg, &theDoubleNormals ) != 1 ) {
	if ( theDoubleNormals.nalloc > 0 ) free( theDoubleNormals.seg );
	return( -1.0 );
      }


      /* find the farthest point outside the sphere
       */
      newlast = l;
      index = _FarthestPointFromSphere( &theSeg, theList, 
					f, &newlast, dim,
					_reduction_mode_in_iterative_ );
      if ( _reduction_mode_in_iterative_ == 1 ) {
	if ( verboseWhenReducing )
	  fprintf( stdout, "...processing frth: remove %d points\n", l-newlast );
	if ( newlast == l ) {
	  suspicion_of_convex_hull = 1;
	  _reduction_mode_of_diameter_ = 0;
	  _reduction_mode_of_dbleNorm_ = 0;
	}
	l = newlast;
      }




      /* stopping condition
	 no point outside the sphere 
      */
      if ( index < f ) {
	if ( theDoubleNormals.nalloc > 0 ) free( theDoubleNormals.seg );
	return( theDiam->squareDiameter );
      }

    } else {

      newEstimateIsSmallerThanCurrentEstimate = 1;

      /*  I add the found segment to the list
	  there is no evidence that it is a maximal segment for
	  the initial set P (but it is for the current one
	  ie the initial minus the points which have been removed),
	  it may somehow help in case of reduction of the set
	  of potential extremities
	  
	  
	  The list of double normals is sorted (from 0 to n)
	  by increasing square diameter:
	  - the best diameter is added (again) at the end
	  - we search the right place for this new element
	  
      */

      if ( _AddSegmentToList( theDiam, &theDoubleNormals ) != 1 ) {
	if ( theDoubleNormals.nalloc > 0 ) free( theDoubleNormals.seg );
	return( -1.0 );
      }
	
      for ( n = theDoubleNormals.n-2; n >= 0 ; n-- ) {
	if ( n == 0 ) {
	  theDoubleNormals.seg[ n ] = theSeg;
	} else {
	  if ( theSeg.squareDiameter <= theDoubleNormals.seg[ n ].squareDiameter &&
	       theSeg.squareDiameter >  theDoubleNormals.seg[ n-1 ].squareDiameter ) {
	    theDoubleNormals.seg[ n ] = theSeg;
	    n = -1;
	  } else {
	    theDoubleNormals.seg[ n ] = theDoubleNormals.seg[ n-1 ];
	  }
	}
      }

    }
    
  } while ( newEstimateIsSmallerThanCurrentEstimate == 0 );










  /* last processing with the found diameter
     - points inside the smallest sphere of the
       diameter may have been already removed
   */
  if ( _reduction_mode_in_iterative_ > 0 && _reduction_mode_of_diameter_ == 1 )
    _reduction_mode_of_diameter_ = 0;



  newlast = l;
  index = _LastPointOutsideSphereWithDiameter( theDiam, theDiam->squareDiameter,
					       theList, f, &newlast, dim,
					       _reduction_mode_of_diameter_ );
  if ( _reduction_mode_of_diameter_ == 1 ||
       _reduction_mode_of_diameter_ == 2 ) {
    if ( verboseWhenReducing )
      fprintf( stdout, "...processing diam: remove %d points\n", l-newlast );
    if ( newlast == l ) {
      suspicion_of_convex_hull = 1;
      _reduction_mode_of_dbleNorm_ = 0;
    }
    l = newlast;
  }
  
  
  /* in some (rare) case, the remaining points outside the largest 
     sphere are removed while searching for a better diameter
     thus it is still an avantageous case.
     if any, we have 
     #f       -> #index : points outside the sphere
     #index+1 -> #l     : points inside the sphere
  */
  if ( index < f ) {
    if ( theDoubleNormals.nalloc > 0 ) free( theDoubleNormals.seg );
    return( theDiam->squareDiameter );
  }







  /* to get some information on 
     the points 
  */
  if ( 0 ) {
    for ( n = theDoubleNormals.n-1; n >= 0; n -- ) {
      _CountPointsInSpheres( &theDoubleNormals.seg[ n ], theDiam->squareDiameter,
			     theList, f, l, dim );
    }
  }







  /* here we will reduce the set of potential extremities
     for the diameter,
     ie the set of points which are to be compared against all 
     the other points

     right now, we have
     #f       -> #index : points outside the sphere
     #index+1 -> #l     : points inside the sphere

     we have a set of maximal segments
     theDoubleNormals.seg[ #i ] for #i from 0 to theDoubleNormals.n-1
     with 
     theDoubleNormals.seg[ theDoubleNormals.n-1 ] == theDiam
     
  */

  if ( tryToReduceQ && theDoubleNormals.n > 1 ) {



    for ( k = 0; k < theDoubleNormals.n; k ++ ) 
      theDoubleNormals.seg[k].reduction_mode = _reduction_mode_of_dbleNorm_;



    switch ( _Q_scan_ ) {
    default :
    case 0 :
      /* backward
       */
      ldn = 0;   fdn = theDoubleNormals.n-2;   idn = -1;
      break;
    case 1 :
      /* forward
       */
      fdn = 0;   ldn = theDoubleNormals.n-2;   idn = +1;
      break;
    }



    for ( n = fdn; n != (ldn+idn) && index >= f ; n += idn ) {
      
      /* in [ #f #index ] find the points outside the sphere
	 theDoubleNormals.seg[ n ]
	 
	 as a result
	 #f   -> #i     are to be compared with all other points
	 #i+1 -> #index are to be compared with a subset
	                if this subset is empty, continue

      */
      i = _LastPointOutsideSphereWithDiameter( &theDoubleNormals.seg[ n ], 
					       theDiam->squareDiameter,
					       theList, f, &index, dim, 0 );
      if ( i >= index ) continue;

      /* in [ #index+1 #l ] find the points outside the sphere
	 theDoubleNormals.seg[ n ]

	 as a result
	 #index+1 -> #j   are to be compared with the previous subset
      */
      newlast = l;
      j = _LastPointOutsideSphereWithDiameter( &theDoubleNormals.seg[ n ], 
					       theDiam->squareDiameter,
					       theList, index+1, &newlast, dim,
					       theDoubleNormals.seg[ n ].reduction_mode );

      if ( theDoubleNormals.seg[ n ].reduction_mode == 1 ||
	   theDoubleNormals.seg[ n ].reduction_mode == 2 ) {

	if ( verboseWhenReducing )
	  fprintf( stdout, "...processing dbNR: remove %d points\n", l-newlast );
	if ( newlast == l ) {
	  suspicion_of_convex_hull = 1;
	  for ( k = 0; k < theDoubleNormals.n; k ++ ) 
	    theDoubleNormals.seg[k].reduction_mode = 0;
	}
	l = newlast;
      }
      
      if ( j <= index ) {
	index = i;
	continue;
      }

      /* right now
	 #f       -> #i     : points to be compared with all other points
	 #i+1     -> #index : points to be compared with the below set
	 #index+1 -> #j     : points to be compared with the above set
	 #j+1     -> #l     : remaining points

	 
      */

      theSeg.extremity1 = (double*)NULL;
      theSeg.extremity2 = (double*)NULL;
      theSeg.squareDiameter = 0.0;

      newEstimate = _QuadraticDiameterInTwoLists( &theSeg, NULL, NULL,
						  theList, i+1, index,
						  theList, index+1, j,
						  dim );
      index = i;
      if ( newEstimate > theDiam->squareDiameter ) {
	/* update variables
	 */
	*theDiam = theSeg;
	/* we find a better estimate
	   it is perhaps not the diameter, according that one
	   diameter extremity can be in [ #f #i ]
	   The question are : 
	   1. have we to look for a maximal segment with these two points 
	      or not ? 
	      -> Seems not necessary ...
	   2. have we to consider this segment with the others ?
	      -> yes
	         but we can not reduce the whole set with it !!!
	      
	*/
	theDoubleNormals.seg[ n ] = theSeg;
	theDoubleNormals.seg[ n ].reduction_mode = 0;
	n -= idn;
      }
      
    }

  }
















  if ( theDoubleNormals.nalloc > 0 ) free( theDoubleNormals.seg );

  /* exhautive search
     
     comparison of points from #f to #index
     against all others points
  */

  if ( dim == 2 ) {
    for ( i=f;   i<=index; i++ )
    for ( j=i+1; j<=l;     j++ ) {
      newEstimate = _SquareDistance2D( theList[i], theList[j] );
      if ( newEstimate > theDiam->squareDiameter ) {
	theDiam->extremity1 = theList[i];
	theDiam->extremity2 = theList[j];
	theDiam->squareDiameter = newEstimate;
      }
    }
    return( theDiam->squareDiameter );
  }


  if ( dim == 3 ) {
    for ( i=f;   i<=index; i++ )
    for ( j=i+1; j<=l;     j++ ) {
      newEstimate = _SquareDistance3D( theList[i], theList[j] );
      if ( newEstimate > theDiam->squareDiameter ) {
	theDiam->extremity1 = theList[i];
	theDiam->extremity2 = theList[j];
	theDiam->squareDiameter = newEstimate;
      }
    }
    return( theDiam->squareDiameter );
  }

  
  for ( i=f;   i<=index; i++ )
  for ( j=i+1; j<=l;     j++ ) {
    newEstimate = _SquareDistance( theList[i], theList[j], dim );
    if ( newEstimate > theDiam->squareDiameter ) {
      theDiam->extremity1 = theList[i];
      theDiam->extremity2 = theList[j];
      theDiam->squareDiameter = newEstimate;
    }
  }
  return( theDiam->squareDiameter );
}



















double _ExactDiameterInTwoLists( typeSegment *theDiam,
				 double **theList1,
				 const int first1,
				 const int last1,
				 double **theList2,
				 const int first2,
				 const int last2,
				 const int dim )
{

  int min_nb_of_points = 5;
  
  int index1;
  int index2;
  int f1 = first1;
  int l1 = last1;
  int f2 = first2;
  int l2 = last2;

  double newEstimate;
  int newEstimateIsSmallerThanCurrentEstimate;
  typeSegment theSeg;

  double dmax1, dmax2;
  
  /*
  theDiam->extremity1 = (double*)NULL;
  theDiam->extremity2 = (double*)NULL;
  theDiam->squareDiameter = 0.0;
  */
  
  if ( first1 < 0 || last1 < 0 ) return( -1.0 );
  if ( first1 > last1 ) {
    l1 = first1;
    f1 = last1;
  }
  if ( first2 < 0 || last2 < 0 ) return( -1.0 );
  if ( first2 > last2 ) {
    l2 = first2;
    f2 = last2;
  }

  
  if ( last1 - first1 <= min_nb_of_points ||
       last2 - first2 <= min_nb_of_points ) {
    return( _QuadraticDiameterInTwoLists( theDiam, NULL, NULL,
					  theList1, first1, last1,
					  theList2, first2, last2,
					  dim ) );
  }
  
  

  index1 = _GetRandomIntNb( first1, last1 );
  index2 = first2 - 1;

  /* at the beginning:
     index1 in [first1, last1]
     
     after only one index_i is in [first_i, last_i]
     it indicates the first point
  */

  do {

    /* end conditions
     */    
    newEstimateIsSmallerThanCurrentEstimate = 0;

    /* find a double normal between the two list
       case #1 : the first point is in list #1
       case #2 : the first point is in list #2
     */
    if ( index2 < f2 ) {
      if ( f2 <= l2 ) {
	newEstimate = _MaximalSegmentInTwoLists( &theSeg, index1,
						 theList1, &f1, &l1,
						 theList2, &f2, &l2, 
						 dim );
	if ( newEstimate > theDiam->squareDiameter ) {
	  theDiam->extremity1     = theSeg.extremity1;
	  theDiam->extremity2     = theSeg.extremity2;
	  theDiam->squareDiameter = theSeg.squareDiameter;
	  newEstimateIsSmallerThanCurrentEstimate = 1;
	}
      }
    } else {
      if ( f1 <= l1 ) {
	newEstimate = _MaximalSegmentInTwoLists( &theSeg, index2,
						 theList2, &f2, &l2, 
						 theList1, &f1, &l1,
						 dim );
	if ( newEstimate > theDiam->squareDiameter ) {
	  theDiam->extremity1     = theSeg.extremity2;
	  theDiam->extremity2     = theSeg.extremity1;
	  theDiam->squareDiameter = theSeg.squareDiameter;
	  newEstimateIsSmallerThanCurrentEstimate = 1;
	}
      }
    }

    /* if we get a better estimation
     */
    if ( newEstimateIsSmallerThanCurrentEstimate == 1 ) {
      
      /* find points outside the sphere of diameter theDiam
       */
      index1 = _FarthestPointFromSphere( theDiam, theList1, f1, &l1, dim, 0 );
      index2 = _FarthestPointFromSphere( theDiam, theList2, f2, &l2, dim, 0 );
      
      /* if there are no points in both list,
	 or in one list
	 it's done
      */
      if ( (index1 <  f1 && index2 < f2) ||
	   (index1 >= f1 &&     l2 < f2) ||
	   (index2 >= f2 &&     l1 < f1) )
	return( theDiam->squareDiameter );
      
      

      if ( index1 >= f1 && index2 >= f2 ) {
	dmax1 = _ScalarProduct( theList1[index1], theDiam->extremity1,
				theList1[index1], theDiam->extremity2, dim ) +
	  0.25 * theDiam->squareDiameter;
	dmax2 = _ScalarProduct( theList2[index2], theDiam->extremity1,
				theList2[index2], theDiam->extremity2, dim ) +
	  0.25 * theDiam->squareDiameter;
	if ( dmax1 > dmax2 ) {
	  index2 = f2 - 1;
	} else {
	  index1 = f1 - 1;
	}
      }
      
    }


  } while ( newEstimateIsSmallerThanCurrentEstimate == 1 );
    

  
  index1 = _LastPointOutsideSphereWithDiameter( theDiam, theDiam->squareDiameter,
						theList1, f1, &l1, dim, 0 );
  index2 = _LastPointOutsideSphereWithDiameter( theDiam, theDiam->squareDiameter,
						theList2, f2, &l2, dim, 0 );

  /* in some (rare) case, the remaining points outside the largest 
     sphere are removed while searching for a better diameter
     thus it is still an avantageous case.
     if any, we have 
     #f       -> #index : points outside the sphere
     #index+1 -> #l     : points inside the sphere
  */

  if ( (index1 <  f1 && index2 < f2) ||
       (index1 >= f1 &&     l2 < f2) ||
       (index2 >= f2 &&     l1 < f1) )
    return( theDiam->squareDiameter );


  if ( index1 >= f1 ) {
    _QuadraticDiameterInTwoLists( theDiam, NULL, NULL,
				  theList1, f1, index1, 
				  theList2, f2, l2,
				  dim );
    if ( index2 >= f2 ) {
      _QuadraticDiameterInTwoLists( theDiam, NULL, NULL,
				    theList1, index1+1, l1,
				    theList2, f2, index2,
				    dim );
    }
    return( theDiam->squareDiameter );
  }

  _QuadraticDiameterInTwoLists( theDiam, NULL, NULL,
				theList1, f1, l1,
				theList2, f2, index2,
				dim );
  return( theDiam->squareDiameter );
}






