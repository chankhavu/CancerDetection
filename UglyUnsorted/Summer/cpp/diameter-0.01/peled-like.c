/*************************************************************************
 *  -
 *
 * $Id: peled-like.c,v 1.1 2002/07/24 12:43:44 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Sat Jun  9 15:23:25 MEST 2001
 *
 *
 * ADDITIONS, CHANGES
 * 
 *
 */



#include <stdio.h>
#include <string.h>

#include <peled-like.h>

#include <print.h>
#include <util.h>
#include <exact-diam.h>


/*
  % make -f Makefile OPT=-g test-diameter 
  % test-diameter -meth peled -dim 2 -p 10   
 */

typedef struct {

  int first;
  int last;

  int n_endpoints;

  double *min;
  double *max;

  int    i_maxlength;
  double l_maxlength;

  /* int diameter; 
   */

  int left;
  int right;

} typeBox;



typedef struct {
  int n_box;
  int a_box;
  typeBox *theBox;
  double  *dbltmp;
} typeListOfBoxes;



#define _N_BOX_  10000



static int add_box_to_list( typeListOfBoxes *t,
			    double **theList,
			    int f, int l, int _dim_ );
static void init_list_of_boxes( typeListOfBoxes *t );
static void free_list_of_boxes( typeListOfBoxes *t );
static void print_box( typeBox *b, int i );
static void update_box( typeBox *b, double **theList, int dim );
static void init_diam_from_box( typeSegment *theDiam,
				typeBox *b,
				double **theList, 
				int dim );











typedef struct {
  
  int box_one;
  int box_two;

  double upperDiameter;

} typePair;


typedef struct {
  int n_pair;
  int a_pair;
  typePair *thePair;
} typeListOfPairs;



#define _N_PAIR_      10000
#define _EXTRA_PAIR_  10000



static int cut_largest_box_in_largest_pair( typeSegment *theDiam, 
					    typeListOfPairs *list_of_pairs,
					    typeListOfBoxes *b,
					    double **oneList, double **twoList,
					    int _dim_,
					    double *current_estimate );


static int add_pair_in_queue( typeSegment *theDiam,
			      typeListOfPairs *t, int b1, int b2, 
			      typeListOfBoxes *b, double **theList1, double **theList2, int dim,
				     double *current_estimate );

static double randomly_update_diameter_from_pair( typeSegment *theDiam,
						  typePair *p, typeListOfBoxes *b, 
						  double **theList, int _dim_ );
static void update_pair( typePair *p, typeListOfBoxes *b, int dim );



static void init_list_of_pairs( typeListOfPairs *t );
static void free_list_of_pairs( typeListOfPairs *t );
static void print_list_of_pairs( typeListOfPairs *t, typeListOfBoxes *b );
static void print_pair( typePair *p, typeListOfBoxes *b, char *s ); 


static void update_diameter_with_a_pair( typeSegment *theDiam,
					 typeBox *b1,
					 typeBox *b2,
					 double **theList1, 
					 double **theList2, 
					 int _dim_,
					 double *current_estimate );











static void _search_pair( typeListOfPairs *p,
			   typeListOfBoxes *b,
			   double **theList,
			   int f,
			   int l,
			   double x1,
			   double y1,
			   double x2,
			   double y2 );
static void _search_point( typeListOfPairs *p,
			   typeListOfBoxes *b,
			   double **theList,
			   int f,
			   int l,
			   double x,
			   double y );


static int _FarthestPointFromSphereAndCount( typeSegment *theSeg,
					     double **theList,
					     const int first,
					     int *last,
					     const int dim,
					     int *nb,
					     const int _reduction_mode_ );









/* 
#define _TRACE_
 */

static int _verbose_ = 1;
void _VerboseInPeledLike()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}
void _NoVerboseInPeledLike()
{
  _verbose_ = 0;
}








/* */
typedef enum {
  _ALLOW_TO_GROW_,
  _REMAINS_AS_IT_IS_
} enumEvolutionMinNbPtsInBox;

static enumEvolutionMinNbPtsInBox _evolution_min_nb_of_pts_ = _ALLOW_TO_GROW_;
#define _MIN_NB_POINTS_IN_A_BOX_ 40
static int min_nb_points_in_a_box = _MIN_NB_POINTS_IN_A_BOX_;

void set_evolution_peled_min_nb_points_to_grow()
{
  _evolution_min_nb_of_pts_ = _ALLOW_TO_GROW_;
}
void set_evolution_peled_min_nb_points_as_it_is()
{
  _evolution_min_nb_of_pts_ = _REMAINS_AS_IT_IS_;
}




/* initialisation du diametre
 */

typedef enum {
  _LARGEST_DIM_,
  _ITERATED_MAXIMAL_SEG_,
  _MAXIMAL_SEG_
} enumInitDiameter;

static enumInitDiameter _init_diameter_ = _LARGEST_DIM_;

void set_init_peled_diameter_to_largest_dim()
{
  _init_diameter_ = _LARGEST_DIM_;
}
void set_init_peled_diameter_to_iterated_maximal_seg()
{
  _init_diameter_ = _ITERATED_MAXIMAL_SEG_;
}
void set_init_peled_diameter_to_maximal_seg()
{
  _init_diameter_ = _MAXIMAL_SEG_;
}






/* calcul du diametre entre 2 boites
 */

typedef enum {
  _QUADRATIC_,
  _SMART_MAX_EXACT_,
  _MAX_EXACT_,
  _MAX_EXACT_WITH_DIAMETER_
} enumUpdateWithPair;

static enumUpdateWithPair _update_with_pair_ = _QUADRATIC_;


void set_update_peled_diameter_to_quadratic()
{
  _update_with_pair_ = _QUADRATIC_;
}
void set_update_peled_diameter_to_smart_max_exact()
{
  _update_with_pair_ = _SMART_MAX_EXACT_;
}
void set_update_peled_diameter_to_max_exact()
{
  _update_with_pair_ = _MAX_EXACT_;
}
void set_update_peled_diameter_to_max_exact_with_diameter()
{
  _update_with_pair_ = _MAX_EXACT_WITH_DIAMETER_;
}












/*
 */
typedef enum {
  _SPLIT_SMALL_BOX_,
  _DO_NOT_SPLIT_SMALL_BOX_,
  _ADAPTED_SPLITTING_
} enumSplitBoxMode;

static enumSplitBoxMode _split_box_mode_ = _SPLIT_SMALL_BOX_;

void adapt_small_box_splitting()
{
  _split_box_mode_ = _ADAPTED_SPLITTING_;
}
void allow_to_split_small_box()
{
  _split_box_mode_ = _SPLIT_SMALL_BOX_;
}
void do_not_allow_to_split_small_box()
{
  _split_box_mode_ = _DO_NOT_SPLIT_SMALL_BOX_;
}











void print_peled_env()
{

  printf( "  " );


  switch ( _evolution_min_nb_of_pts_ ) {
  case _ALLOW_TO_GROW_ :
    printf( "nb_pts=can-grow " ); break;
  case _REMAINS_AS_IT_IS_ :
    printf( "nb_pts=stable   " ); break;
  }

  switch ( _init_diameter_ ) {
  case _LARGEST_DIM_ :
    printf( "diam-init=large-dim    " ); break;
  case _MAXIMAL_SEG_ :
    printf( "diam-init=max-seg      " ); break;
  case _ITERATED_MAXIMAL_SEG_ :
    printf( "diam-init=iter-max-seg " ); break;
  }

  switch ( _update_with_pair_ ) {
  case _QUADRATIC_ :
    printf( "update=quadratic       " ); break;
  case _SMART_MAX_EXACT_ :
    printf( "update=smart-max-exact " ); break;
  case _MAX_EXACT_ :
    printf( "update=max-exact       " ); break;
  case _MAX_EXACT_WITH_DIAMETER_ :
    printf( "update=max-exact-diam  " ); break;
  }
 
  switch( _split_box_mode_ ) {
  case _SPLIT_SMALL_BOX_ :
    printf( "split-small-box=yes " ); break;
  case _DO_NOT_SPLIT_SMALL_BOX_ :
    printf( "split-small-box=no  " ); break;
  case _ADAPTED_SPLITTING_ :
    printf( "adapted-splitting   " ); break;
  }
  
  printf( "\n" );

}





















static int nb_pair;


double _Peled_LikeDiameterInOneList( typeSegment *theDiam,
				     double **theList,
				     const int first,
				     const int last,
				     const int dim )
{
  char *proc = "_Peled_LikeDiameterInOneList";
  typeSegment exact_pair1;
  double exact_diameter;

  typeListOfBoxes list_of_boxes;
  typeListOfPairs list_of_pairs;

  int f = first;
  int l = last;
  int index = 0;
  
  double init_estimate  = 0.0;
  double current_estimate = 0.0;

  int i = 0;

  double **oneList = NULL;
  double **twoList = NULL;

  int _reduction_mode_ = 1;
  typeSegment minBall, theSeg;
  int outside, minOutside = last-first+1;
  int newEstimateIsSmallerThanCurrentEstimate;
  double newEstimate;
  int newlast;

  
  /* global
     initialisation
  */
  nb_pair = 0;
  min_nb_points_in_a_box = _MIN_NB_POINTS_IN_A_BOX_;




  if ( 0 ) {
    exact_diameter = _ExactDiameterInOneList( &exact_pair1, theList,
					      f, l, dim );
    _PrintSegment( stdout, &exact_pair1, "Points realizing the diameter (exact)", dim );
    printf( "\n" );
  }


  theDiam->extremity1 = (double*)NULL;
  theDiam->extremity2 = (double*)NULL;
  theDiam->squareDiameter = 0.0;






  init_list_of_boxes( &list_of_boxes );
  init_list_of_pairs( &list_of_pairs );


  /* oneList is for the potential ends of diameter
     twoList is for the other points
  */



  switch( _init_diameter_ ) {
  default :
  case _LARGEST_DIM_ :

    oneList = twoList = theList;
    if ( add_box_to_list( &list_of_boxes, theList, f, l, dim ) < 0 ) {
      fprintf( stderr, "%s: error in inserting box in list\n", proc );
      return( 0.0 );
    }
    list_of_boxes.theBox[0].n_endpoints = l-f+1;
    init_diam_from_box( theDiam, &(list_of_boxes.theBox[0]), theList, dim );
    init_estimate = current_estimate = sqrt( theDiam->squareDiameter );
    (void)add_pair_in_queue( theDiam, &list_of_pairs, 0, 0, 
			     &list_of_boxes, oneList, twoList, dim, &current_estimate );
    break;


  case _ITERATED_MAXIMAL_SEG_ :

    index = _GetRandomIntNb( first, last );
    i = 0;

    do {

      newEstimateIsSmallerThanCurrentEstimate = 0;
      newEstimate = _MaximalSegmentInOneList( &theSeg, index, theList,
					      &f, &l, dim );

      if ( newEstimate > theDiam->squareDiameter ) {

	*theDiam = theSeg;
	newlast = l;
	index = _FarthestPointFromSphereAndCount( &theSeg, theList, 
						  f, &newlast, dim, &outside,
						  _reduction_mode_ );
	if ( minOutside > outside ) {
	  minBall = theSeg;
	  minOutside = outside;
	}
	if ( _reduction_mode_ == 1 ) {
	  if ( newlast == l ) _reduction_mode_ = 0;
	  l = newlast;
	}
	if ( index < f ) {
	  fprintf( stderr, "%s: the maximal segment is the diameter\n", proc );
	  fprintf( stderr, "     no further processing needed\n" );
	  return( theDiam->squareDiameter );
	}

      } else {

	newEstimateIsSmallerThanCurrentEstimate = 1;

      }

    } while ( newEstimateIsSmallerThanCurrentEstimate == 0 );

    init_estimate = current_estimate = sqrt( theDiam->squareDiameter );
    index = _LastPointOutsideSphereWithDiameter( &minBall, theDiam->squareDiameter,
						 theList, f, &l, dim,
						 _reduction_mode_ );

  case _MAXIMAL_SEG_ :

    if ( _init_diameter_ == _MAXIMAL_SEG_ ) {

      index = _GetRandomIntNb( f, l );
      init_estimate = current_estimate 
	= sqrt( _MaximalSegmentInOneList( theDiam, index, theList,
					  &f, &l, dim ) );
      index = _LastPointOutsideSphereWithDiameter( theDiam, theDiam->squareDiameter,
						   theList, f, &l, dim,
						   _reduction_mode_ );
    }

    /* from f       to index => outside
       form index+1 to l     => inside 
    */
    if ( index < f ) {
      fprintf( stderr, "%s: the maximal segment is the diameter\n", proc );
      fprintf( stderr, "     no further processing needed\n" );
      return( theDiam->squareDiameter );
    }
    oneList = (double**)malloc( (index - f + 1) * sizeof( double* ) );
    if ( oneList == NULL ) {
      fprintf( stderr, "%s: error when allocating list\n", proc );
      return( 0.0 );
    }
    (void)memcpy( oneList, &(theList[f]), (index - f + 1) * sizeof( double* ) );
    twoList = theList;
    if ( add_box_to_list( &list_of_boxes, oneList, 0, index-f, dim ) < 0 ) {
      fprintf( stderr, "%s: error in inserting box in list\n", proc );
      free( oneList );
      return( 0.0 );
    }
    list_of_boxes.theBox[0].n_endpoints = (index - f + 1);
    if ( add_box_to_list( &list_of_boxes, twoList, f, l, dim ) < 0 ) {
      fprintf( stderr, "%s: error in inserting box in list\n", proc );
      free_list_of_boxes( &list_of_boxes );
      free( oneList );
      return( 0.0 );
    }
    list_of_boxes.theBox[1].n_endpoints = 0;
    (void)add_pair_in_queue( theDiam, &list_of_pairs, 0, 1, 
			      &list_of_boxes, oneList, twoList, dim, &current_estimate );

    break;
  }






  


  if ( _verbose_ ) {
    fprintf( stdout, "%s:\n", proc );
    fprintf( stdout, "     init diameter = %f\n",
	     init_estimate );
    fprintf( stdout, "     list = [%d %d] X [%d %d]\n",
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_one ].first,
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_one ].last,
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_two ].first,
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_two ].last );
    fprintf( stdout, "     pot. ends = {%d %d}\n", 
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_one ].n_endpoints,
	     list_of_boxes.theBox[ list_of_pairs.thePair[0].box_two ].n_endpoints );
  }











  while ( list_of_pairs.n_pair > 0 ) {


    
#ifdef _TRACE_ 
      i++;
      printf( "#%2d estimate = %f\n", i, current_estimate );
#endif





    /* coupe la boite et nouvelle queue
     */
    (void) cut_largest_box_in_largest_pair( theDiam, &list_of_pairs, 
					    &list_of_boxes, 
					    oneList, twoList,
					    dim, &current_estimate );
    
  }


#ifdef _TRACE_ 
    printf( "\n" );
    printf( "estimate = %f\n", current_estimate );
    print_list_of_pairs( queue, &list_of_boxes );


    _PrintSegment( stdout, &exact_pair1, "Points realizing the diameter (exact)", dim );
    _PrintSegment( stdout, theDiam, "Points realizing the diameter (peled like)", dim );

    printf( "\n" );
#endif


  if ( _verbose_ ) {
    fprintf( stdout, "%s: use %d boxes, %d pairs\n", 
	     proc, list_of_boxes.n_box, nb_pair );
  }



  free_list_of_pairs( &list_of_pairs );
  free_list_of_boxes( &list_of_boxes );

  switch( _init_diameter_ ) {
  default :
  case _LARGEST_DIM_ :
    break;
  case _MAXIMAL_SEG_ :
    free( oneList );
  }

  return( theDiam->squareDiameter );
}
















static void remove_first_pair_from_queue( typeListOfPairs *list_of_pairs ) 
{
  int  ind, left, right, max_ind;
  typePair pair;
  
  list_of_pairs->n_pair --;
  list_of_pairs->thePair[ 0 ] = list_of_pairs->thePair[ list_of_pairs->n_pair ];
  ind = 0;
    
  while  ( ind < list_of_pairs->n_pair ) {
    left = 2 * ind + 1;
    right = 2 * ind + 2;
    if  ( left >= list_of_pairs->n_pair ) 
      break;
    if  ( right >= list_of_pairs->n_pair ) 
      right = left;
    
    if ( list_of_pairs->thePair[ left ].upperDiameter <
	 list_of_pairs->thePair[ right ].upperDiameter )
      max_ind = right;
    else
      max_ind = left;
    if ( list_of_pairs->thePair[ ind ].upperDiameter >
	 list_of_pairs->thePair[ max_ind ].upperDiameter )
      break;
    
    pair = list_of_pairs->thePair[ind];
    list_of_pairs->thePair[ind] = list_of_pairs->thePair[max_ind];
    list_of_pairs->thePair[max_ind] = pair;
    
    ind = max_ind;
  }
}




static int cut_box( int i_box,
		     typeListOfBoxes *list_of_boxes, double **theList,
		     int _dim_ )
{
  char *proc = "cut_box";
  int i_largest_dim;
  int left, right;
  double threshold;

  
  /* already cut
   */
  if ( list_of_boxes->theBox[ i_box ].left  >= 0 || 
       list_of_boxes->theBox[ i_box ].right >= 0 )
    return( 1 );
  
  /* contains only one point
   */
  if ( list_of_boxes->theBox[ i_box ].last - list_of_boxes->theBox[ i_box ].first <= 0 ) {
    list_of_boxes->theBox[ i_box ].left  = i_box;
    list_of_boxes->theBox[ i_box ].right = -1;
    return( 1 );
  }

  i_largest_dim = list_of_boxes->theBox[ i_box ].i_maxlength;
  
  /* tri de la liste de points
   */
  left  = list_of_boxes->theBox[ i_box ].first;
  right = list_of_boxes->theBox[ i_box ].last;
  
  threshold = ( list_of_boxes->theBox[ i_box ].min[ i_largest_dim ] + 
		list_of_boxes->theBox[ i_box ].max[ i_largest_dim ] ) / 2.0;

    
  while ( left < right ) {
    if ( theList[left ][i_largest_dim] < threshold ) {
      left ++;
    }
    else if ( theList[right][i_largest_dim] >= threshold ) {
      right--; 
    }
    else {
      _SwapPoints( theList, right, left );
    }
  }

  if ( theList[left ][i_largest_dim] >= threshold ) {
    left --;
  } 
  else if ( theList[right][i_largest_dim] < threshold ) {
    right ++;
  }





  if ( add_box_to_list( list_of_boxes,
			theList, list_of_boxes->theBox[ i_box ].first, left, _dim_ ) < 0  ) {
    fprintf( stderr, "%s: error in inserting first box\n", proc );
    return( -1 );
  }
  list_of_boxes->theBox[ i_box ].left  = list_of_boxes->n_box - 1;


  if ( add_box_to_list( list_of_boxes, 
			theList, right, list_of_boxes->theBox[ i_box ].last, _dim_ ) < 0 ) {
    fprintf( stderr, "%s: error in inserting second box\n", proc );
    return( -1 );
  }

  list_of_boxes->theBox[ i_box ].right = list_of_boxes->n_box - 1;

  
  return( 1 );
}






/* diviser une boite
   trouver ses occurences
   pour chacune d'entre 
      mettre la paire dans la liste 
      ou calculer un diametre exact
 */










static int cut_largest_box_in_largest_pair( typeSegment *theDiam, 
					    typeListOfPairs *list_of_pairs,
					    typeListOfBoxes *list_of_boxes,
					    double **oneList,
					    double **twoList,
					    int _dim_,
					    double *current_estimate )
{
  char *proc = "cut_largest_box_in_largest_pair";
  int one_left, one_right;
  int two_left, two_right;
  int i_one, i_two;



  
  if ( list_of_pairs->n_pair <= 0 ) return( 1 );





  if (  list_of_pairs->thePair[0].upperDiameter <= *current_estimate ) {
    remove_first_pair_from_queue( list_of_pairs );
    return( 1 );
  }




  i_one = list_of_pairs->thePair[0].box_one;
  i_two = list_of_pairs->thePair[0].box_two;

  if (  list_of_boxes->theBox[ i_one ].last - list_of_boxes->theBox[ i_one ].first
	<= min_nb_points_in_a_box &&
	list_of_boxes->theBox[ i_two ].last - list_of_boxes->theBox[ i_two ].first
	<= min_nb_points_in_a_box ) {

    update_diameter_with_a_pair( theDiam,
				 &(list_of_boxes->theBox[ i_one ]),
				 &(list_of_boxes->theBox[ i_two ]),
				 oneList, twoList, _dim_,
				 current_estimate );

    remove_first_pair_from_queue( list_of_pairs );
    return( 1 );

  }




  switch( _split_box_mode_ ) {
  case _ADAPTED_SPLITTING_ :

    if ( (list_of_boxes->theBox[ i_two ].last - list_of_boxes->theBox[ i_two ].first > min_nb_points_in_a_box
	  && (list_of_boxes->theBox[ i_one ].last - list_of_boxes->theBox[ i_one ].first > min_nb_points_in_a_box 
	      || list_of_boxes->theBox[ i_one ].l_maxlength >= list_of_boxes->theBox[ i_two ].l_maxlength)) ||
	 (list_of_boxes->theBox[ i_one ].last - list_of_boxes->theBox[ i_one ].first > min_nb_points_in_a_box
	  && (list_of_boxes->theBox[ i_two ].last - list_of_boxes->theBox[ i_two ].first > min_nb_points_in_a_box 
	      || list_of_boxes->theBox[ i_two ].l_maxlength >= list_of_boxes->theBox[ i_one ].l_maxlength)) ) {
	 
      
      if ( cut_box( i_two, list_of_boxes, twoList, _dim_ ) < 0 ) {
	fprintf( stderr, "%s: error in spliting second box\n", proc );
	return( -1 );
      }
      two_left  = list_of_boxes->theBox[ i_two ].left;
      two_right = list_of_boxes->theBox[ i_two ].right;
      
      
      if ( cut_box( i_one, list_of_boxes, oneList, _dim_ ) < 0 ) {
	fprintf( stderr, "%s: error in spliting first box\n", proc );
	return( -1 );
      }
      one_left  = list_of_boxes->theBox[ i_one ].left;
      one_right = list_of_boxes->theBox[ i_one ].right;
      
      remove_first_pair_from_queue( list_of_pairs );
      
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_left, two_left, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_left, two_right, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_right, two_left, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_right, two_right, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      
      return( 1 );
    }


    if ( list_of_boxes->theBox[ i_one ].last - list_of_boxes->theBox[ i_one ].first > min_nb_points_in_a_box ) {

      if ( cut_box( i_one, list_of_boxes, oneList, _dim_ ) < 0 ) {
	fprintf( stderr, "%s: error in spliting first box\n", proc );
	return( -1 );
      }
      one_left  = list_of_boxes->theBox[ i_one ].left;
      one_right = list_of_boxes->theBox[ i_one ].right;
      
      remove_first_pair_from_queue( list_of_pairs );
      
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_left, i_two, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      (void)add_pair_in_queue( theDiam, list_of_pairs, 
			       one_right, i_two, list_of_boxes, 
			       oneList, twoList, _dim_, current_estimate );
      
      return( 1 );
    }

    

    if ( cut_box( i_two, list_of_boxes, twoList, _dim_ ) < 0 ) {
      fprintf( stderr, "%s: error in spliting second box\n", proc );
      return( -1 );
    }
    two_left  = list_of_boxes->theBox[ i_two ].left;
    two_right = list_of_boxes->theBox[ i_two ].right;
      
      
    remove_first_pair_from_queue( list_of_pairs );
      
    (void)add_pair_in_queue( theDiam, list_of_pairs, 
			     i_one, two_left, list_of_boxes, 
			     oneList, twoList, _dim_, current_estimate );
      
    (void)add_pair_in_queue( theDiam, list_of_pairs, 
			     i_one, two_right, list_of_boxes, 
			     oneList, twoList, _dim_, current_estimate );
    
    return( 1 );
    


    break;
  case _DO_NOT_SPLIT_SMALL_BOX_ :
    
    /* cette solution semble
       - creer moins de boites 
       - creer a peu pres autant de paires
       - etre plus rapide (avec les procedures rusees de calcul entre boites
       % test-diameter -mode comp -p 20000 -dim 10 -init 993057538
    */

    if (  list_of_boxes->theBox[ i_one ].last - list_of_boxes->theBox[ i_one ].first >
	  min_nb_points_in_a_box ) {
      
      if ( cut_box( i_one, list_of_boxes, oneList, _dim_ ) < 0 ) {
	fprintf( stderr, "%s: error in spliting first box\n", proc );
	return( -1 );
      }
      one_left  = list_of_boxes->theBox[ i_one ].left;
      one_right = list_of_boxes->theBox[ i_one ].right;
      
    } else {
      one_left  = i_one;
      one_right = -1;
    }
    
    if (  list_of_boxes->theBox[ i_two ].last - list_of_boxes->theBox[ i_two ].first
	  > min_nb_points_in_a_box ) {
      
      if ( cut_box( i_two, list_of_boxes, twoList, _dim_ ) < 0 ) {
	fprintf( stderr, "%s: error in spliting second box\n", proc );
	return( -1 );
      }
      two_left  = list_of_boxes->theBox[ i_two ].left;
      two_right = list_of_boxes->theBox[ i_two ].right;
      
    }  else {
      two_left  = i_two;
      two_right = -1;
    }
    
    break;

  case _SPLIT_SMALL_BOX_ :
  default :
    
    /* cette solution semble
       - creer plus de boites 
       - creer a peu pres autant de paires
       - etre plus lente (avec les procedures rusees de calcul entre boites
       % test-diameter -mode comp -p 20000 -dim 10 -init 993057538
    */

    if ( cut_box( i_one, list_of_boxes, oneList, _dim_ ) < 0 ) {
      fprintf( stderr, "%s: error in spliting first box\n", proc );
      return( -1 );
    }
    one_left  = list_of_boxes->theBox[ i_one ].left;
    one_right = list_of_boxes->theBox[ i_one ].right;
    
    if ( cut_box( i_two, list_of_boxes, twoList, _dim_ ) < 0 ) {
      fprintf( stderr, "%s: error in spliting second box\n", proc );
      return( -1 );
    }
    two_left  = list_of_boxes->theBox[ i_two ].left;
    two_right = list_of_boxes->theBox[ i_two ].right;
    
  }
  



  remove_first_pair_from_queue( list_of_pairs );






  (void)add_pair_in_queue( theDiam, list_of_pairs, 
			   one_left, two_left, list_of_boxes, 
			   oneList, twoList, _dim_, current_estimate );

  (void)add_pair_in_queue( theDiam, list_of_pairs, 
			   one_left, two_right, list_of_boxes, 
			    oneList, twoList, _dim_, current_estimate );
  if ( one_right != two_right )
    (void)add_pair_in_queue( theDiam, list_of_pairs, 
			     one_right, two_left, list_of_boxes, 
			     oneList, twoList, _dim_, current_estimate );
  (void)add_pair_in_queue( theDiam, list_of_pairs, 
			   one_right, two_right, list_of_boxes, 
			   oneList, twoList, _dim_, current_estimate );

  return( 1 );
}















static int add_pair_in_queue( typeSegment *theDiam,
			      typeListOfPairs *t, int b1, int b2, 
			      typeListOfBoxes *b, double **oneList, double **twoList, int dim,
			      double *current_estimate )
{
  char *proc = "add_pair_in_queue";
  typePair test;
  typePair pair, *tmpPair = NULL;
  double d, d2;
  int ind, n;


  if ( b1 < 0 || b2 < 0 ) return( 1 );

  if ( b->theBox[b1].n_endpoints == 0 && b->theBox[b2].n_endpoints == 0 )
    return( 1 );


  if ( b->theBox[b1].last - b->theBox[b1].first <= min_nb_points_in_a_box &&
       b->theBox[b2].last - b->theBox[b2].first <= min_nb_points_in_a_box ) {

    update_diameter_with_a_pair( theDiam,
				 &(b->theBox[b1]),
				 &(b->theBox[b2]),
				 oneList, twoList, 
				 dim,
				 current_estimate );
    return( 1 );
  }

  nb_pair ++;

  
  if ( t->n_pair >= t->a_pair ) {
    
    if ( t->a_pair <= 0 ) {
      n = _N_PAIR_;
    } else {
      n = t->a_pair + _EXTRA_PAIR_;
      if ( _evolution_min_nb_of_pts_ == _ALLOW_TO_GROW_ ) {
	min_nb_points_in_a_box *= 2;
	if ( _verbose_ )
	  fprintf( stdout, "minimal nb of points in box set to %d\n", min_nb_points_in_a_box );
      }
    }
    
    tmpPair = (typePair*)calloc( n, sizeof( typePair ) );
    if ( tmpPair == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate list of pair\n", proc );
      return( -1 );
    }
    
    if ( t->a_pair > 0 ) {
      (void)memcpy( tmpPair, t->thePair, t->a_pair  * sizeof( typePair ) );
      if ( t->thePair != NULL ) free( t->thePair );
    }
    
    t->thePair = tmpPair;
    t->a_pair  = n;

  }




  test.box_one = b1;
  test.box_two = b2;
  update_pair( &(test), b, dim );
  
  if ( test.upperDiameter <= *current_estimate ) return( 1 );


  
  /*
    d = randomly_update_diameter_from_pair( theDiam, &test, b, theList, dim );
  */
  d2 = _SquareDistance( oneList[ b->theBox[ test.box_one ].first ],
			twoList[ b->theBox[ test.box_two ].first ], dim ) ;
  d = sqrt( d2 );
  if ( *current_estimate < d )  {
    theDiam->extremity1 = oneList[ b->theBox[ test.box_one ].first ];
    theDiam->extremity2 = twoList[ b->theBox[ test.box_two ].first ];
    theDiam->squareDiameter = d2;
    *current_estimate = d;
  }




  ind = t->n_pair;
  t->n_pair++;

  t->thePair[ind] = test;
  while  ( ind > 0 ) {
    n = ( ind - 1 ) / 2;
    if ( t->thePair[ind].upperDiameter < t->thePair[n].upperDiameter )
      break;
    pair = t->thePair[ind];
    t->thePair[ind] = t->thePair[n];
    t->thePair[n] = pair;
    ind = n;
  }

  
  return( 1 );
}











/* pick a pair of points from two and computes the length
 */

static double randomly_update_diameter_from_pair( typeSegment *theDiam,
						  typePair *p, typeListOfBoxes *b, 
						  double **theList, int _dim_ )
{

  typeBox *b_one, *b_two;
  int i1, i2;
  double d;

  b_one = &(b->theBox[ p->box_one ]);
  b_two = &(b->theBox[ p->box_two ]);

  i1 = ( b_one->last > b_one->first ) ? _GetRandomIntNb( b_one->first, b_one->last ) : b_one->first;
  i2 = ( b_two->last > b_two->first ) ? _GetRandomIntNb( b_two->first, b_two->last ) : b_two->first;

  d = _SquareDistance( theList[i1], theList[i2], _dim_  );
  if ( d > theDiam->squareDiameter ) {
      theDiam->extremity1 = theList[i1];
      theDiam->extremity2 = theList[i2];
      theDiam->squareDiameter = d;
  }
  
  return( sqrt( d ) );

}








static void update_pair( typePair *p, typeListOfBoxes *b, int dim )
{
  typeBox *b_one, *b_two;
  int i;
  double max, min, d, u;

  b_one = &(b->theBox[ p->box_one ]);
  b_two = &(b->theBox[ p->box_two ]);

  for ( u=0, i=0; i<dim; i++ ) {
    max = ( b_one->max[i] > b_two->max[i] ) ? b_one->max[i] : b_two->max[i];
    min = ( b_one->min[i] < b_two->min[i] ) ? b_one->min[i] : b_two->min[i];
    d = max - min;
    u += d*d;
  }
  p->upperDiameter = sqrt( u );
  
}



static void init_list_of_pairs( typeListOfPairs *t )
{
  t->n_pair = 0;
  t->a_pair = 0;
  t->thePair = NULL;
}

static void free_list_of_pairs( typeListOfPairs *t )
{
  if ( t->thePair != NULL ) free( t->thePair );
  init_list_of_pairs( t );
}




static void print_list_of_pairs( typeListOfPairs *t, typeListOfBoxes *b )
{

  char str[10];
  int i = 0;

  for (i=0; i < t->n_pair; i ++ ) {
    sprintf( str, "%d", i );
    print_pair( &(t->thePair[i]), b, str );
  }
}



static void print_pair( typePair *p, typeListOfBoxes *b, char *s ) 
{
  printf( "PAIR" );
  if ( s != NULL ) printf( "[%s]", s );
  printf( " = (%d %d)", p->box_one, p->box_two );
  printf( " , D < %f", p->upperDiameter );
  printf( "\n" );
  /*
  printf( " %f < D < %f , %f < D < %f \n",
	  p->my_lowerDiameter, p->my_upperDiameter,
	  p->randomDiameter, p->upperDiameter );
  */
  
  if ( b != NULL ) {
    print_box( &(b->theBox[ p->box_one ]), p->box_one );
    print_box( &(b->theBox[ p->box_two ]), p->box_two );
  }
}
























static int add_box_to_list( typeListOfBoxes *t,
			    double **theList,
			    int f, int l, int _dim_ )
{
  char *proc = "add_box_to_list";
  typeBox *tmpBox = NULL;
  double *tmpDbl = NULL;
  int i, n;

  if ( t->n_box >= t->a_box ) {
    for ( n=t->a_box ; n <= t->n_box ; n += _N_BOX_ )
      ;
    
    
    tmpBox = (typeBox *)calloc( n, sizeof( typeBox ) );
    if ( tmpBox == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate list of box\n", proc );
      return( -1 );
    }
    tmpDbl = (double*)calloc( 2 * n * _dim_,  sizeof( double ) );
    if ( tmpDbl == NULL ) {
    free( tmpBox );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate list of pointers\n", proc );
      return( -1 );
    }

    if ( t->a_box > 0 ) {
      (void)memcpy( tmpBox, t->theBox,     t->a_box * sizeof( typeBox ) );
      (void)memcpy( tmpDbl, t->dbltmp, 2 * t->a_box * _dim_ * sizeof( double ) );
      if ( t->theBox != NULL ) free( t->theBox );
      if ( t->dbltmp != NULL ) free( t->dbltmp );
    }
    for ( i=0; i<n; i ++ ) {
      tmpBox[i].min = &(tmpDbl[(    i*2) * _dim_]);
      tmpBox[i].max = &(tmpDbl[(1 + i*2) * _dim_]);
    }
    t->theBox = tmpBox;
    t->dbltmp = tmpDbl;
    t->a_box = n;
  }
  
  t->theBox[t->n_box].first = f;
  t->theBox[t->n_box].last  = l;

  t->theBox[t->n_box].n_endpoints = -1;

  t->theBox[t->n_box].left = -1;
  t->theBox[t->n_box].right = -1;

  update_box( &(t->theBox[t->n_box]), theList, _dim_ );
  t->n_box ++;

  return( t->n_box - 1 );
}



static void init_list_of_boxes( typeListOfBoxes *t )
{
  t->n_box = 0;
  t->a_box = 0;
  t->theBox = NULL;
  t->theBox = NULL;
}

static void free_list_of_boxes( typeListOfBoxes *t )
{
  if ( t->theBox != NULL ) free( t->theBox );
  if ( t->dbltmp != NULL ) free( t->dbltmp );
  init_list_of_boxes( t );
}




static void print_box( typeBox *b, int s )
{
  printf( "    BOX" );
  printf( "[%3d]", s );
  printf( "= {#%d->#%d} , ",
	  b->first, b->last );
  /*
    printf( "l_max=%f along i_max=%d , R=%f\n",
	  b->l_maxlength, b->i_maxlength, b->diameter / 2.0 );
  */
  printf( "l_max=%f along i_max=%d\n",
	  b->l_maxlength, b->i_maxlength );
}







/* met a jour les parametres de la boite
   sachant que les indices des premiers
   et derniers points sont donnes
*/
static void update_box( typeBox *b, double **theList, int dim )
{
  int i, j;
  /* double r; */
  double d;



  /* coins inferieurs et superieurs
   */
  memcpy( b->min, theList[b->first], dim * sizeof( double ) );
  memcpy( b->max, theList[b->first], dim * sizeof( double ) );



  for ( i=b->first+1; i<=b->last; i++ ) {
    for ( j =0; j <dim ; j ++ ) {
      if ( b->min[j] > theList[i][j] ) {
	b->min[j] = theList[i][j];
      } else if ( b->max[j] < theList[i][j] ) {
	b->max[j] = theList[i][j];
      }
    }
  }


  /* milieu
     longueur max
     rayon
  */

  b->i_maxlength = 0;
  b->l_maxlength = b->max[0] - b->min[0];

  d = b->max[0] - b->min[0];
  
  /* r = d*d; */

  for ( j=1; j <dim ; j ++ ) {
    d = b->max[j] - b->min[j];
    /* r += d*d; */
    if ( b->l_maxlength < d ) {
      b->l_maxlength = d;
      b->i_maxlength = j;
    }
  }
  /* b->diameter  =  sqrt( r ); */
}



static void init_diam_from_box( typeSegment *theDiam,
				typeBox *b,
				double **theList, 
				int dim )
{
  int j_max=0, i, j;
  double d, d_max;
  double *pt_f, *pt_l;

  d_max = b->max[0] - b->min[0];
  for ( j=1; j <dim ; j ++ ) {
    d = b->max[j] - b->min[j];
    if ( d > d_max ) {
      d_max = d;
      j_max = j;
    }
  }
  
  pt_f = pt_l = theList[b->first];
  for ( i=b->first+1; i<=b->last; i++ ) {
    if ( pt_f[j_max] > theList[i][j_max] ) {
      pt_f = theList[i];
    } else if ( pt_l[j_max] < theList[i][j_max] ) {
      pt_l = theList[i];
    }
  }
  
  theDiam->extremity1 = pt_f;
  theDiam->extremity2 = pt_l;
  theDiam->squareDiameter = _SquareDistance( pt_f, pt_l, dim );
}






















static void update_diameter_with_a_pair( typeSegment *theDiam,
					 typeBox *b1,
					 typeBox *b2,
					 double **oneList, 
					 double **twoList, 
					 int _dim_,
					 double *current_estimate )
{
  int min_nb_of_points = 5;
  
  typeSegment theSeg;
  double d;
  enumUpdateWithPair update = _update_with_pair_;

  double **theList1 = NULL;
  int f1 = b1->first;
  int l1 = b1->last;

  double **theList2 = NULL;
  int f2 = b2->first;
  int l2 = b2->last;
  
  int i, n;

  theSeg.extremity1 = NULL;
  theSeg.extremity2 = NULL;
  theSeg.squareDiameter = 0;



  /* si les boites contiennent peu de points
   */
  
  if ( b1->last - b1->first <= min_nb_of_points ||
       b2->last - b2->first <= min_nb_of_points ) {
    d =_QuadraticDiameterInTwoLists( &theSeg, NULL, NULL,
				     oneList, b1->first, b1->last,
				     twoList, b2->first, b2->last,
				     _dim_ );
    update = _QUADRATIC_;
  }

  
  
  /* si les boites ont deja ete coupees
   */
  if ( update == _QUADRATIC_ ) {
    theList1 = oneList;
    theList2 = twoList;
  } else {

    if ( b1->left >= 0 && b1->right >= 0 ) {
      theList1 = (double**)malloc( (b1->last - b1->first + 1)*sizeof( double* ) );
      (void)memcpy( theList1, &(oneList[b1->first]), (b1->last - b1->first + 1)*sizeof( double* ) );
      f1 = 0;
      l1 = b1->last - b1->first;
    } else {
      theList1 = oneList;
    }

    if ( b2->left >= 0 && b2->right >= 0 ) {
      theList2 = (double**)malloc( (b2->last - b2->first + 1)*sizeof( double* ) );
      (void)memcpy( theList2, &(twoList[b2->first]), (b2->last - b2->first + 1)*sizeof( double* ) );
      f2 = 0;
      l2 = b2->last - b2->first;
    } else {
      theList2 = twoList;
    }

  }




  switch ( update ) {
  default :
  case _QUADRATIC_ :
    
    d = _QuadraticDiameterInTwoLists( &theSeg, NULL, NULL,
				      oneList, b1->first, b1->last,
				      twoList, b2->first, b2->last,
				      _dim_ );

    
    break;


  case _SMART_MAX_EXACT_ :
    
    /* we count le number of points outside the sphere
       if they are too many points (we count only in the first
       list, which could be the list of potential diameter
       extremity), compute a "new" diameter
       else, initialise with the current diameter
    */

    switch ( _dim_ ) {
    case 3 :
      for ( n=0, i=b1->first; i <= b1->last; i++ ) {
	if ( _ScalarProduct3D( oneList[i], theDiam->extremity1,
			       oneList[i], theDiam->extremity2  ) > 0 )
	  n ++;
      }
      for ( i=b2->first; i <= b2->last; i++ ) {
	if ( _ScalarProduct3D( twoList[i], theDiam->extremity1,
			       twoList[i], theDiam->extremity2  ) > 0 )
	  n ++;
      }
      break;
    default :
      for ( n=0, i=b1->first; i <= b1->last; i++ ) {
	if ( _ScalarProduct( oneList[i], theDiam->extremity1,
			     oneList[i], theDiam->extremity2, _dim_  ) > 0 )
	  n ++;
      }
      for ( i=b2->first; i <= b2->last; i++ ) {
	if ( _ScalarProduct( twoList[i], theDiam->extremity1,
			     twoList[i], theDiam->extremity2, _dim_  ) > 0 )
	  n ++;
      }
    }
  
    
    if ( 8*n > (b1->last - b1->first + 1) + (b2->last - b2->first + 1) ) {
      d = _ExactDiameterInTwoLists( &theSeg, 
				    theList1, f1, l1, 
				    theList2, f2, l2,
				    _dim_ );
      break;
    }

    
   case _MAX_EXACT_WITH_DIAMETER_ :

    *current_estimate = sqrt( _ExactDiameterInTwoLists( theDiam,
							theList1, f1, l1, 
							theList2, f2, l2,
							_dim_ ) );
    
    if ( theList1 != oneList ) free( theList1 );
    if ( theList2 != twoList ) free( theList2 );
    return;

  case _MAX_EXACT_ :

    d = _ExactDiameterInTwoLists( &theSeg, 
				  theList1, f1, l1, 
				  theList2, f2, l2,
				  _dim_ );
    
    break;

  }
  
  if ( theDiam->squareDiameter < theSeg.squareDiameter ) {
    theDiam->squareDiameter = theSeg.squareDiameter;
    theDiam->extremity1 = theSeg.extremity1;
    theDiam->extremity2 = theSeg.extremity2;
    *current_estimate = sqrt( theSeg.squareDiameter );
  }

  if ( theList1 != oneList ) free( theList1 );
  if ( theList2 != twoList ) free( theList2 );

}


















static void _search_pair( typeListOfPairs *p,
			   typeListOfBoxes *b,
			   double **theList,
			   int f,
			   int l,
			   double x1,
			   double y1,
			   double x2,
			   double y2 )
{
  int i, j, n1, n2;
  double e = 0.00001;

  for ( n1=f-1, i=f; i<=l && n1<f; i++ ) {
    if ( x1-e < theList[i][0] && 
	 x1+e > theList[i][0] &&
	 y1-e < theList[i][1] && 
	 y1+e > theList[i][1] ) {
      n1 = i;
    }
  }
  fprintf( stdout, " POINT[%4d] = [%f %f %f]\n", n1, theList[n1][0], theList[n1][1], theList[n1][2] );
  
  for ( n2=f-1, i=f; i<=l && n2<f; i++ ) {
    if ( x2-e < theList[i][0] && 
	 x2+e > theList[i][0] &&
	 y2-e < theList[i][1] && 
	 y2+e > theList[i][1] ) {
      n2 = i;
    }
  }
  fprintf( stdout, " POINT[%4d] = [%f %f %f]\n",
	   n2, theList[n2][0], theList[n2][1], theList[n2][2] );

  for ( j=0; j<p->n_pair; j++ ) {
    if ( (b->theBox[ p->thePair[j].box_one ].first <= n1 &&
	  b->theBox[ p->thePair[j].box_one ].last  >= n1 &&
	  b->theBox[ p->thePair[j].box_two ].first <= n2 &&
	  b->theBox[ p->thePair[j].box_two ].last  >= n2) ||
	 (b->theBox[ p->thePair[j].box_one ].first <= n2 &&
	  b->theBox[ p->thePair[j].box_one ].last  >= n2 &&
	  b->theBox[ p->thePair[j].box_two ].first <= n1 &&
	  b->theBox[ p->thePair[j].box_two ].last  >= n1) ) {
      if ( 0 ) {
	fprintf( stdout, " PAIR[%2d] = {%2d %2d}",
		 j, p->thePair[j].box_one, p->thePair[j].box_two );
	fprintf( stdout, " = { [%3d %3d] [%3d %3d] }\n",
		 b->theBox[ p->thePair[j].box_one ].first,
		 b->theBox[ p->thePair[j].box_one ].last,
		 b->theBox[ p->thePair[j].box_two ].first,
		 b->theBox[ p->thePair[j].box_two ].last );
      }
      fprintf( stdout, " PAIR[%2d] = ", j );
      print_pair( &(p->thePair[j]), b, NULL );
    }
    
  }
  
}


static void _search_point( typeListOfPairs *p,
			   typeListOfBoxes *b,
			   double **theList,
			   int f,
			   int l,
			   double x,
			   double y )
{
  int i, j, pa, n, bo;
  double e = 0.00001;
  for ( n=f-1, i=f; i<=l && n<f; i++ ) {
    if ( x-e < theList[i][0] && 
	 x+e > theList[i][0] &&
	 y-e < theList[i][1] && 
	 y+e > theList[i][1] ) {
      n = i;
    }
  }
  if ( n == f-1 ) return;
  fprintf( stdout, " POINT[%4d] = [%f %f %f]\n",
	   n, theList[n][0], theList[n][1], theList[n][2] );
  
  for ( bo=-1 , i=0; i<b->n_box; i++ ) {
    if ( n >= b->theBox[i].first &&
	 n <= b->theBox[i].last ) {
      fprintf( stdout, "   in BOX[%4d] = [%3d %3d]\n", 
	       i, b->theBox[i].first, b->theBox[i].last  );
      for ( pa=-1 , j=0; j<p->n_pair; j++ ) {
	if ( p->thePair[j].box_one == i ||
	     p->thePair[j].box_two == i )
	  fprintf( stdout, "      in PAIR[%4d] = {%3d %3d}\n",
		   j, p->thePair[j].box_one, p->thePair[j].box_two );
      }
    }
  }
    


}














static int _FarthestPointFromSphereAndCount( typeSegment *theSeg,
					     double **theList,
					     const int first,
					     int *last,
					     const int dim,
					     int *nb,
					     const int _reduction_mode_ ) 
{
  int i, l = (*last);
  int index = first - 1;
  double diff, maxdiff = 0.0;
  /* threshold = 1.5 - sqrt(3)
   */
  double threshold = -0.23205080756887729352 * theSeg->squareDiameter;


  if ( l < first ) return( index );
  *nb = 0;

  switch ( _reduction_mode_ ) {

  case 0 :

    /* NO REDUCTION CASE
     */

    if ( dim == 2 ) {
      for ( i=first; i<=l; i++ ) {
	diff = _ScalarProduct2D( theList[i], theSeg->extremity1,
				 theList[i], theSeg->extremity2 );
	if ( diff > 0 ) (*nb) ++;
	if ( maxdiff < diff ) {
	  index   = i;
	  maxdiff = diff;
	}
      }
      return( index );
    }
    
    if ( dim == 3 ) {
      for ( i=first; i<=l; i++ ) {
	diff = _ScalarProduct3D( theList[i], theSeg->extremity1,
				 theList[i], theSeg->extremity2 );
	if ( diff > 0 ) (*nb) ++;
	if ( maxdiff < diff ) {
	  index   = i;
	  maxdiff = diff;
	}
      }
    return( index );
    }

    for ( i=first; i<=l; i++ ) {
      diff = _ScalarProduct( theList[i], theSeg->extremity1,
			     theList[i], theSeg->extremity2, dim );
      if ( diff > 0 ) (*nb) ++;
      if ( maxdiff < diff ) {
	index   = i;
	maxdiff = diff;
      }
    }
    return( index );

    /* END
       NO REDUCTION CASE
    */
    break;
    
  default :
  case 1 :

    /* REDUCTION INSIDE THE SMALLEST SPHERE
     */

    /* AB = diameter extremities
       MA.MB = (MC+CA).(MC+CB) = MC^2 + CA.CB + MC ( CB+CA )
       = MC^2 - R^2   + 0
    */
    if ( dim == 2 ) {
      for ( i=first; i<=l; i++ ) {
	diff = _ScalarProduct2D( theList[i], theSeg->extremity1,
				 theList[i], theSeg->extremity2 );
	if ( diff > 0 ) (*nb) ++;
	if ( diff > maxdiff ) {
	  index   = i;
	  maxdiff = diff;
	} 
	else if ( diff <= threshold ) {
	  _SwapPoints( theList, i, l );
	  i --;   l --;
	  continue;
	}
      }
      *last = l;
      return( index );
    }
    
    
    if ( dim == 3 ) {
      for ( i=first; i<=l; i++ ) {
	diff = _ScalarProduct3D( theList[i], theSeg->extremity1,
				 theList[i], theSeg->extremity2 );
	if ( diff > 0 ) (*nb) ++;
	if ( maxdiff < diff ) {
	  index   = i;
	  maxdiff = diff;
	}
	else if ( diff <= threshold ) {
	  _SwapPoints( theList, i, l );
	  i --;   l --;
	  continue;
	}
      }
      *last = l;
      return( index );
    }
    
    
    for ( i=first; i<=l; i++ ) {
      diff = _ScalarProduct( theList[i], theSeg->extremity1,
			     theList[i], theSeg->extremity2, dim );
      if ( diff > 0 ) (*nb) ++;
      if ( maxdiff < diff ) {
	index   = i;
	maxdiff = diff;
      }
      else if ( diff <= threshold ) {
	_SwapPoints( theList, i, l );
	i --;   l --;
	continue;
      }
    }
    *last = l;
    return( index );
    /* END
       REDUCTION INSIDE THE SMALLEST SPHERE
     */
  }
  
  return( -1 );  
}
