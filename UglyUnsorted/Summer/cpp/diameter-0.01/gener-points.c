/*************************************************************************
 *  -
 *
 * $Id: gener-points.c,v 1.3 2004/06/10 09:23:33 greg Exp $
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

#include <stdio.h>
#include <string.h>
#include <time.h>

#include <pick.h>
#include <rand.h>
#include <read.h>
#include <util.h>

#include <print.h>

#define _STRLENGTH_ 256


typedef enum {
  _TEST_,
  _EVAL_,
  _TIME_
} enumMode;

typedef enum {
  _EXACT_,
  _QUADR_,
  _APPRX_,
  _PELEG_
} enumMethod;


static char *program;
static char *usage = "\
 [-dim %d]\n\
 [-distrib|-dist  \n\tc[ube]|e[llipse]|r[ellipse]|g[ellipse]|s[phere]|b[all]|d[iamcst]]\n\
 [-help | --help | -h]\
 [-init | -seed %ld]\
 [-model %s]\n\
 [-no-verbose | -nv]\
 [-out %s]\
 [-ply %s]\
 [-points | -p %d]\
 [-pts %s]\n\
 [-random-rotation | -rr]\
 [-verbose | -v]\
 [-vertices %d]\n\
 [%s]\
\n";

static char *detail = "\
-dim %d # dimension of the embedding space\n\
-distrib | -dist %s # allows to choose a random distribution\n\
                    # %s should be among\n\
     c[ube]     # distribution inside a cube\n\
     e[ellipse] # distribution on an ellipsoid\n\
     r[ellipse] # distribution on a regular/gentle ellipsoid\n\
     g[ellipse] # id.\n\
     s[phere]   # distribution on a sphere\n\
     b[all]     # distribution inside a sphere\n\
     d[iamcst]  # distribution inside a set on constant width (must be 2-D)\n\
-help | --help | -h #\n\
-init | -seed %ld # set the initial seed for the sequence of\n\
     # pseudo-random integers\n\
-model %s # obsolete\n\
-no-verbose | -nv #\n\
%s | -out %s # writes distribution in file '%s' (else on stdout)\n\
-ply %s # read points from an ascii ply file, as the ones available at\n\
     # http://www.cs.gatech.edu/projects/large_models/\n\
     # home made routine, not robust at all\n\
-points | -p %d # number of generated points for a given random distribution\n\
-pts %s # read points from a file as the ones generated from 'gener-points'\n\
-random-rotation | -rr # apply a random rotation to the points of a random\n\
     # distribution. Works only in 2-D and 3-D.\n\
-verbose | -v #\n\
-vertices %d # gives the number of vertices of Reuleaux polygon\n\
     # (only with '-dim 2')\n\
     # '-vertices 1' -> triangle (2*1+1) = 3\n\
     # '-vertices 2' -> pentagon (2*2+1) = 5\n\
\n";

static char *explain="\nGenerates list of points\n";



static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s %s\n",program, usage, explain);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s\n",str);
  exit(0);
}






int main( int argc, char *argv[] )
{
  char modlname[_STRLENGTH_];
  char *outname = NULL;
  FILE *fout = NULL;

  enumDistribution typeDistribution = IN_CUBE;
#ifdef WIN32
  unsigned int seedRandom = time(0);
#else
  long int seedRandom = time(0);
#endif


  int _dim_ = 3;  
  double _max_diameter_ = 1.0;
  double _min_diameter_ = 0.2;
  int    _psommet_ = 1;
  int    _nbpoints_ = 10;

  int i, status;

  double **listOfPoints = (double**)NULL;


  program = argv[0];

  for ( i=1; i<argc; i++ ) {

    
    if ( argv[i][0] == '-' ) {
      if ( strcmp ( argv[i], "-help" ) == 0 ||
	   strcmp ( argv[i], "--help" ) == 0 ||
	   strcmp ( argv[i], "-h" ) == 0 ) {
	_ErrorParse( "help message\n", 1 );
      }
      
      
      /* distributions
       */
      if ( (strcmp ( argv[i], "-distrib" ) == 0) ||
	   (strcmp ( argv[i], "-dist" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -distrib...\n", 0 );
	if ( (strcmp ( argv[i], "cube" ) == 0) ||
	     (strcmp ( argv[i], "c" ) == 0) ) {
	  typeDistribution = IN_CUBE;
	}
	else if ( (strcmp ( argv[i], "ellipse" ) == 0) ||
		  (strcmp ( argv[i], "e" ) == 0) ) {
	  typeDistribution = ON_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "gellipse" ) == 0) ||
		  (strcmp ( argv[i], "g" ) == 0) ) {
	  typeDistribution = ON_REG_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "rellipse" ) == 0) ||
		  (strcmp ( argv[i], "r" ) == 0) ) {
	  typeDistribution = ON_REG_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "sphere" ) == 0) ||
		  (strcmp ( argv[i], "s" ) == 0) ) {
	  typeDistribution = ON_SPHERE;
	}
	 else if ( (strcmp ( argv[i], "ball" ) == 0) ||
		   (strcmp ( argv[i], "b" ) == 0) ) {
	   typeDistribution = IN_SPHERE;
	 } 
	else if ( (strcmp ( argv[i], "diamcst" ) == 0) ||
		  (strcmp ( argv[i], "d" ) == 0) ) {
	  typeDistribution = IN_CST_DIAMETER;
	} 
	else  {
	  sprintf( modlname, "unknown distribution = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      } 
    
      
      else if ( strcmp ( argv[i], "-out" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -out...\n", 0 );
	outname = argv[i];
      }
      else if ( strcmp ( argv[i], "-model" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -model...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_MODEL;
      }
      else if ( strcmp ( argv[i], "-ply" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -ply...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_PLY_MODEL;
      }
      else if ( strcmp ( argv[i], "-pts" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -pts...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_PTS_MODEL;
      }
      

      /* number of points
       */
      else if ( (strcmp ( argv[i], "-p" ) == 0) ||
		(strcmp ( argv[i], "-points" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -points...\n", 0 );
	status = sscanf( argv[i],"%d",&_nbpoints_ );
	if ( status <= 0 ) _ErrorParse( "parsing -points...\n", 0 );
      }

      /* dimension
       */
      else if ( strcmp ( argv[i], "-dim" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -dim...\n", 0 );
	status = sscanf( argv[i],"%d",&_dim_ );
	if ( status <= 0 ) _ErrorParse( "parsing -dim...\n", 0 );
      }

      /* seed for random distribution
       */
      else if ( (strcmp ( argv[i], "-init" ) == 0) ||
		(strcmp ( argv[i], "-seed" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -seed...\n", 0 );
#ifdef WIN32
	status = sscanf( argv[i],"%u",&seedRandom );
#else
	status = sscanf( argv[i],"%ld",&seedRandom );
#endif
	if ( status <= 0 ) _ErrorParse( "parsing -seed...\n", 0 );
      }

      /* number of vertices for Reuleaux polygons = 2*p+1
	 _psommet_ = 1 => triangle.
       */
      else if ( strcmp ( argv[i], "-vertices" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -vertices...\n", 0 );
	status = sscanf( argv[i],"%d",&_psommet_ );
	if ( status <= 0 ) _ErrorParse( "parsing -vertices...\n", 0 );
      }
      

      else if ( (strcmp ( argv[i], "-random-rotation" ) == 0) ||
		(strcmp ( argv[i], "-rr" ) == 0) ) {
	_VerboseInPick();
	_DoApplyARandomRotation();
      }



      else if ( (strcmp ( argv[i], "-verbose" ) == 0) ||
		(strcmp ( argv[i], "-v" ) == 0) ) {
	_VerboseInPick(); 
	_VerboseInRead();
	
      }
      else if ( (strcmp ( argv[i], "-no-verbose" ) == 0) ||
		(strcmp ( argv[i], "-nv" ) == 0) ) {
	_NoVerboseInPick(); 
	_NoVerboseInRead();
      }


      /*--- option inconnue ---*/
      else {
	sprintf( modlname, "unknown option %s\n", argv[i] );
	_ErrorParse( modlname, 0);
      }
      
      
    } 

    /* arg without '-'
     */
    else if ( argv[i][0] != 0 ) {
      if ( outname == NULL ) 
	outname = argv[i];
      else {
	sprintf( modlname, "too many file names at '%s'\n", argv[i] );
	_ErrorParse( modlname, 0);
      }
    }
  }



  _SetRandomSeed( seedRandom );

  switch( typeDistribution ) {
  default :
    listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					  _min_diameter_, _psommet_, 
					  typeDistribution, _dim_ );
    break;
  case IN_PLY_MODEL:
    listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
    break;
  case IN_MODEL:
    listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
    break;
  case IN_PTS_MODEL:
    listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
    break;
  }

  if ( listOfPoints == NULL ) {
    _ErrorParse( "unable to generate/read list of points\n", 0 );
  }


  _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );



  if ( outname != NULL ) {
    fout = fopen( outname, "w" );
  }
  if ( fout != NULL ) {
    fprintf( fout, "#\n" );
    fprintf( fout, "#" );
    for ( i=0; i<argc; i++ ) fprintf( fout, " %s", argv[i] );;
    fprintf( fout, "\n" );
    fprintf( fout, "#\n" );
  }
  if ( fout != NULL ) {
#ifdef WIN32
    fprintf( fout, "# Random seed = %u\n", seedRandom );
#else
    fprintf( fout, "# Random seed = %ld\n", seedRandom );
#endif
    fprintf( fout, "#\n" );
  }
  else {
    if ( 0 )
#ifdef WIN32
      fprintf( stderr, "Random seed = %u\n", seedRandom );
#else
      fprintf( stderr, "Random seed = %ld\n", seedRandom );
#endif
  }

  _PrintListOfPoints( listOfPoints, 0, _nbpoints_-1, _dim_, fout );

  if ( fout != NULL ) fclose ( fout );
  
  free( listOfPoints );

  return( 0 );
}
