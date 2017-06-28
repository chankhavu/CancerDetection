/*************************************************************************
 *  -
 *
 * $Id: print.c,v 1.1 2002/07/24 12:43:44 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue May 16 2000
 *
 *
 * ADDITIONS, CHANGES
 * 
 *
 */

#include <print.h>






void _PrintEnv( FILE *f,
		const long int init,
		const int nbpoints,
		const double diameterMax,
		const double diameterMin,
		const int psommet,
		enumDistribution typeDistribution,
		char *modelname,
		const int dim )
{
  fprintf( f, "\n" );
  fprintf( f, "%%   random seed      = %ld\n", init );
  fprintf( f, "%%   number of points = %d\n", nbpoints );
  fprintf( f, "%%   dimension        = %d\n", dim );
  fprintf( f, "%%   try to reduce P in iterative = %d\n", _GetReductionModeInIterative() );
  fprintf( f, "%%   try to reduce P in diameter  = %d\n", _GetReductionModeOfDiameter() );
  fprintf( f, "%%   try to reduce Q     = %d\n", _GetTryToReduceQ() );
  if ( _GetTryToReduceQ() ) {
    fprintf( f, "%%   try to reduce P in dble norm = %d\n", _GetReductionModeOfDbleNorm() );
    switch( _GetQscan() ) {
    default :
    case 0 :
      fprintf( f, "%%   backward scan\n" );
      break;
    case 1 :
      fprintf( f, "%%   forward scan\n" );
    }
  }
  fprintf( f, "%%   try to get tight bounds = %d\n", _GetTightBounds() );
  fprintf( f, "\n" );

  fprintf( f, "%%   distribution      = " );
  switch( typeDistribution ) {
  case IN_CUBE :          fprintf( f, "in a cube" ); break;
  case ON_SPHERE :        fprintf( f, "on a ball" ); break;
  case IN_SPHERE :        fprintf( f, "in a ball" ); break;
  case ON_ELLIPSOID :     fprintf( f, "on an ellipsoid" ); break;
  case ON_REG_ELLIPSOID : fprintf( f, "on a regular ellipsoid" ); break;
  case IN_CST_DIAMETER :  fprintf( f, "in a Reuleaux polygon (vertices=%d)", 2*psommet+1 ); break;
  case IN_PLY_MODEL :         fprintf( f, "from a ply model '%s'", modelname ); break;
  case IN_MODEL :            fprintf( f, "from a     model '%s'", modelname ); break;
  case IN_PTS_MODEL :         fprintf( f, "from a pts model '%s'", modelname ); break;
  default :
    fprintf( f, "unknown" );
  }
  fprintf( f, "\n" );

  switch( typeDistribution ) {
  case ON_ELLIPSOID :  
  case ON_REG_ELLIPSOID :
    fprintf( f, "%%   small diameter    = %g\n", diameterMin );
  case IN_CUBE :
  case ON_SPHERE :
  case IN_SPHERE :
  case IN_CST_DIAMETER :
    fprintf( f, "%%   diameter          = %g\n", diameterMax );
    break;
  default :
    break;
  }
  fprintf( f, "\n" );
 
}






void _PrintSegment( FILE *fp, typeSegment *seg, char *desc, int dim )
{
  double *ext1 = seg->extremity1;
  double *ext2 = seg->extremity2;

  if ( desc != NULL ) fprintf( fp, "%s (length = %g)\n", 
			       desc, sqrt( seg->squareDiameter ) );

  if ( seg->extremity1 == NULL ||seg->extremity2 == NULL ) {

    fprintf( fp, "   (NULL)\n" );
    return;
  }
    
  if ( ext1[0] > ext2[0] ) {
    ext1 = seg->extremity2;
    ext2 = seg->extremity1;
  }
  
  fprintf( fp, "  " );


  if ( dim >=2 ) fprintf( fp, "(%g,%g", ext1[0], ext1[1] );
  if ( dim == 2 ) {
    fprintf( fp, ")" );
  } else if (dim == 3) {
    fprintf( fp, ",%g)", ext1[2] );
  } else {
    fprintf( fp, ",%g,...)", ext1[2] );
  }


  fprintf( fp, " - " );
    

  if ( dim >=2 ) fprintf( fp, "(%g,%g", ext2[0], ext2[1] );
  if ( dim == 2 ) {
    fprintf( fp, ")" );
  } else if (dim == 3) {
    fprintf( fp, ",%g)", ext2[2] );
  } else {
    fprintf( fp, ",%g,...)", ext2[2] );
  }

  fprintf( fp, "\n" );
}
		    






void _PrintListOfPoints( double **theList,
			 const int first,
			 const int last,
			 const int dim,
			 FILE *f ) 
{
  int i, j;
  FILE *fp = f;
  if ( fp == NULL ) fp = stdout;

  if ( theList == NULL || first > last ) return;

  fprintf( fp, "%d", last-first+1 );
  fprintf( fp, " %d", dim );
  fprintf( fp, "\n" );

  for (i=first; i<=last; i++ ) {
    for ( j=0; j<dim; j++ ) {
      fprintf( fp, "%11g", theList[i][j] );
      if ( j < dim-1 ) fprintf( fp, " " );
    }
    fprintf( fp, "\n" );
  }
}



int _PrintPointIndex( double **theList,
		      const int first,
		      const int last,
		      const double i1,
		      const double i2,
		      const int dim ) 
{
  int i, k =-1;
  double e = 0.0001;
  for (i=first; i<=last; i++ ) {
    if ( theList[i][0] - e < i1 && i1 < theList[i][0] + e &&
	 theList[i][1] - e < i2 && i2 < theList[i][1] + e ) {
      k = i;
      fprintf( stderr, "%d : (%g,%g", i, theList[i][0], theList[i][1] );
      if ( dim == 2 ) {
	fprintf( stderr, ")\n" );
      } 
      else if ( dim == 3 ) {
	fprintf( stderr, ",%g)\n", theList[i][2] );
      }
      else {
	fprintf( stderr, ",%g,...)\n", theList[i][2] );
      }
    }
  }
  return( k );
}
