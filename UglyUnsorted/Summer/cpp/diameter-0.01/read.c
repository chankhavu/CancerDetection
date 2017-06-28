/*************************************************************************
 *  -
 *
 * $Id: read.c,v 1.1 2002/07/24 12:43:44 greg Exp $
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




#include <read.h>

static int _verbose_ = 0;

void _VerboseInRead()
{
  _verbose_ = 1;
}

void _NoVerboseInRead()
{
  _verbose_ = 0;
}






void *_ReadPLYModel( char *name, int *nbpoints, int *dim  )
{
  char *proc = "_ReadPLYModel";
  FILE *f;
  int stillread = 1;
  char string[100];
  int l=0, n=0;
  
  double **listOfPoints = (double**)NULL;


  *dim = 3;
  *nbpoints = 0;

  f = fopen( name, "r" );
  if ( f == NULL ) return( NULL );

  do {

    l++;

    if ( fgets( string, 99, f ) == NULL ) {
      stillread = 0;
      continue;
    }
    
    if ( strncmp( string, "element vertex ", 15 ) == 0 ) {
      sscanf( &(string[15]), "%d", nbpoints );
      
      listOfPoints = (double**)_AllocateListOfPoints( *nbpoints, *dim );
      if ( listOfPoints == NULL ) return( NULL );
      continue;
    }

    if ( string[0] >= 'a' && string[0] <= 'z' ) {
      continue;
    }

    if ( listOfPoints != NULL ) {
      if ( sscanf( string, "%lf %lf %lf", &(listOfPoints[n][0]), 
		   &(listOfPoints[n][1]), &(listOfPoints[n][2]) ) == 3 )
	n++;
    }

    if ( n == *nbpoints ) stillread = 0;

  } while ( stillread == 1 );

  if ( _verbose_ )
    fprintf( stderr, "%s: has read %d points (dim=%d) in '%s'\n",
	     proc, *nbpoints, *dim, name );

  fclose(f);

  return( (void*)listOfPoints );
}





#define STR_LENGTH 1024

void *_ReadPTSModel( char *name, int *nbpoints, int *dim  )
{
  char *proc = "_ReadPTSModel";
  FILE *f;
  int stillread = 1;
  int i, d, n=0;
  char *string;

  double **listOfPoints = (double**)NULL;


  *dim = 3;
  *nbpoints = 0;

  f = fopen( name, "r" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s'\n", proc, name );
    return( NULL );
  }

  string = (char*)malloc( STR_LENGTH );
  if ( string == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate string while reading '%s'\n", 
	     proc, name );
    fclose( f );
    return( NULL );
  }


  /* skip comments
   */
  do {
    if ( fgets( string, STR_LENGTH, f ) == NULL ) {
      if ( _verbose_ )
      fprintf( stderr, "%s: unable to read string from '%s'\n", 
	       proc, name );
      stillread = 0;
      continue;
    }
    for ( i=0; i<STR_LENGTH 
	    && string[i] == ' ' 
	    && string[i] == '\t' ; i++ )
      ;
    if ( string[i] == '#' ) 
      continue;
    stillread = 0;
  } while ( stillread == 1 );
  
  /* number of points
     and dim
  */
  for ( i=0; i<STR_LENGTH 
	  && string[i] == ' ' 
	  && string[i] == '\t' ; i++ )
    ;
  
  (void)sscanf( &(string[i]), "%d", nbpoints );
  while ( '0' <= string[i] && string[i] <= '9' ) i++;
  while ( ' ' == string[i] || string[i] == '\t' ) i++;
  if ( '0' <= string[i] && string[i] <= '9' ) {
    (void)sscanf( &(string[i]), "%d", dim );
  }
  
  listOfPoints = (double**)_AllocateListOfPoints( *nbpoints, *dim );
  if ( listOfPoints == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate points while reading '%s'\n", 
	     proc, name );
    free( string );
    fclose( f );
    return( NULL );
  }
  

  stillread = 1;
  n = 0;
  d = 0;
  do {
    if ( fscanf( f, "%lg", &(listOfPoints[n][d]) ) == 1 ) {
      d ++;
      if ( d == *dim ) { 
	d = 0;
	n ++;
	if ( n == *nbpoints ) stillread = 0;
      }
    } 
    else {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: '%s' seems truncated\n", 
		 proc, name );
	fprintf( stderr, "   stop reading at coordinate #%d of point #%d,",
		 d, n );
	fprintf( stderr, " should have read %d points (dim=%d) \n", 
		 *nbpoints, *dim );
      }
      free( listOfPoints );
      free( string );
      fclose( f );
      return( NULL );
    } 
  } while ( stillread == 1 );
  

  if ( _verbose_ )
    fprintf( stderr, "%s: has read %d points (dim=%d) in '%s'\n",
	   proc, *nbpoints, *dim, name );

  free( string );
  fclose(f);

  return( (void*)listOfPoints );

}






void *_ReadModel( char *name, int *nbpoints, int *dim  )
{
  FILE *f;
  int stillread = 1;
  char *str, string[100];
  int i, l=0, n=0;
  long int seed, call;
  
  double **listOfPoints = (double**)NULL;


  *dim = 0;
  *nbpoints = 0;

  f = fopen( name, "r" );
  if ( f == NULL ) return( NULL );

  do {

    l++;

    if ( fgets( string, 99, f ) == NULL ) {
      stillread = 0;
      continue;
    }
    
    if ( string[0] == '#' ) continue;

    if ( strncmp( string, "points ", 7 ) == 0 ) {
      sscanf( string, "points = %d", nbpoints );
      
      if ( (*nbpoints) > 0 && (*dim) > 0 && listOfPoints == NULL ) {
	listOfPoints = (double**)_AllocateListOfPoints( *nbpoints, *dim );
	if ( listOfPoints == NULL ) return( NULL );
	for ( n=0; n<(*nbpoints); n++ )
	  for (i=0; i<(*dim); i++ ) 
	    listOfPoints[n][i] = 0.0;
	n = 0;
      }
      continue;
    }

    if ( strncmp( string, "dim ", 4 ) == 0 ) {
      sscanf( string, "dim = %d", dim );
      
      if ( (*nbpoints) > 0 && (*dim) > 0 && listOfPoints == NULL ) {
	listOfPoints = (double**)_AllocateListOfPoints( *nbpoints, *dim );
	if ( listOfPoints == NULL ) return( NULL );
	for ( n=0; n<(*nbpoints); n++ )
	  for (i=0; i<(*dim); i++ ) 
	    listOfPoints[n][i] = 0.0;
	n = 0;
      }
      continue;
    }

    if ( strncmp( string, "random seed ", 12 ) == 0 ) {
      sscanf( string, "random seed = %20ld", &seed );
      _SetRandomSeed( seed );
      continue;
    }

    if ( strncmp( string, "random calls ", 13 ) == 0 ) {
      sscanf( string, "random calls = %20ld", &call );
      if ( 0 ) 
	for (i=0; i<call; i++ )
	  (void)_GetRandomDoubleNb();
      continue;
    }

    if ( string[0] >= 'a' && string[0] <= 'z' ) {
      continue;
    }
    
    if ( listOfPoints != NULL ) {
      if ( sscanf( string, "%lf", &(listOfPoints[n][0]) ) == 1 ) {
	str = string;
	while( (str[0] >= '0' && str[0] <= '9') || str[0] == '-' || str[0] =='.' ) {
	  str++;
	}
	for (i=1; i<(*dim); i++ ) {
	  (void)sscanf( str, "%lf", &(listOfPoints[n][i]) );
	  str ++;
	  while( (str[0] >= '0' && str[0] <= '9') || str[0] == '-' || str[0] =='.' ) {
	    str++;
	  }
	}
	n++;
      }
    }

    if ( n == *nbpoints ) stillread = 0;

  } while ( stillread == 1 );
  fclose(f);

  return( (void*)listOfPoints );
}
