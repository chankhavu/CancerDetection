#include <stdio.h>
#include <string.h>
#include <time.h>

#include <pick.h>
#include <rand.h>
#include <read.h>
#include <util.h>
#include <exact-diam.h>
#include <apprx-diam.h>

#include <peled-like.h>

#include <print.h>

#define _STRLENGTH_ 256

typedef enum {
  _EXACT_,
  _QUADR_,
  _APPRX_,
  _PELED_
} enumMethod;

static char program[_STRLENGTH_];

static char *usage = "\
 [-dim %d]\n\
 [-distrib|-dist  \n\tc[ube]|e[llipse]|r[ellipse]|g[ellipse]|s[phere]|b[all]|d[iamcst]]\n\
 [-epsilon | -e %lf]\
 [-help | --help | -h]\
 [-init | -seed %ld]\n\
 [-method|-meth exact|appr[o]x|quadr[atic]|[har-]peled|hybrid1|hybrid2]\n\
 [-model %s]\
 [-no-verbose | -nv]\
 [-ply %s]\
 [-points | -p %d]\
 [%s | -pts %s]\n\
 [-Q-reduction]\
 [-random-rotation | -rr]\
 [-verbose | -v]\
 [-vertices %d]\
 [%s]\
\n";

static char *explain="\nComputes the diameter (or an approximation) of a point set\n";

static char *detail1 = "\
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
-epsilon | -e %lf # epsilon for the epsilon-approximation of the diameter\n\
-help | --help | -h #\n\
-init | -seed %ld # set the initial seed for the sequence of\n\
     # pseudo-random integers\n\
-method | -meth %s # allows to choose the method of computation\n\
     exact       # computation of the exact diameter, Malandain & Boissonnat\n\
     approx      # computation of a epsilon-approximation,\n\
                 # Malandain & Boissonnat\n\
     brute-force # computation of the exact diameter, brute-force method\n\
     brute       # id.\n\
     quadratic   # id.\n\
     quadr       # id.\n\
     har-peled   # computation of the exact diameter, Har-Peled\n\
     peled       # id.\n\
     hybrid1     # first hybrid method (see related publications)\n\
     hybrid2     # second hybrid method (see related publications)\n\
-model %s # obsolete\n\
-no-verbose | -nv #\n";
static char *detail2 = "\
-ply %s # read points from an ascii ply file, as the ones available at\n\
     # http://www.cs.gatech.edu/projects/large_models/\n\
     # home made routine, not robust at all\n\
-points | -p %d # number of generated points for a given random distribution\n\
%s | -pts %s # read points from a file as the ones generated from 'gener-points'\n\
-Q-reduction # reduces the last set of potential diameter endpoints\n\
     # (see related publications)\n\
-random-rotation | -rr # apply a random rotation to the points of a random\n\
     # distribution. Works only in 2-D and 3-D.\n\
-verbose | -v #\n\
-vertices %d # gives the number of vertices of Reuleaux polygon\n\
     # (only with '-dim 2')\n\
     # '-vertices 1' -> triangle (2*1+1) = 3\n\
     # '-vertices 2' -> pentagon (2*2+1) = 5\n\
\n";



static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s %s\n",program, usage, explain);
  if ( flag == 1 ) (void)fprintf(stderr,"%s%s",detail1,detail2);
  (void)fprintf(stderr,"Erreur : %s\n",str);
  exit(0);
}






int main( int argc, char *argv[] )
{
  char modlname[_STRLENGTH_];

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
  int    _nbpoints_ = 1000;
  double _epsilon_  = 0.01;
  enumMethod typeMethod = _EXACT_;

  double **listOfPoints = (double**)NULL;
  
  typeSegment pair1;
  int c0, c1;
  int status, i;
  double upper = 0.0;
  
  /* some default tuning
   */
  _SetReductionModeInIterative( 0 );
  _SetReductionModeOfDiameter( 0 );
  _SetReductionModeOfDbleNorm( 0 );
  _DoNotTryToReduceQ();



  /* parse args
   */
  for ( i=1; i<argc; i++ ) {

    if ( argv[i][0] == '-' ) {
      if ( strcmp ( argv[i], "-help" ) == 0 ||
	   strcmp ( argv[i], "--help" ) == 0 ||
	   strcmp ( argv[i], "-h" ) == 0 ) {
	_ErrorParse( "help message\n", 1 );
      }

      

      /* file 
       */
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


      /* random distributions
       */
      else if ( (strcmp ( argv[i], "-distrib" ) == 0) ||
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
	else if ( (strcmp ( argv[i], "rellipse" ) == 0) ||
		  (strcmp ( argv[i], "r" ) == 0) ) {
	  typeDistribution = ON_REG_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "gellipse" ) == 0) ||
		  (strcmp ( argv[i], "g" ) == 0) ) {
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



      /* number of points
       */
      else if ( (strcmp ( argv[i], "-p" ) == 0) ||
		(strcmp ( argv[i], "-points" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -points...\n", 0 );
	status = sscanf( argv[i],"%d",&_nbpoints_ );
	if ( status <= 0 ) _ErrorParse( "parsing -points...\n", 0 );
      }

      /* epsilon for approximation
       */
      else if ( (strcmp ( argv[i], "-e" ) == 0) ||
		(strcmp ( argv[i], "-epsilon" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -epsilon...\n", 0 );
	status = sscanf( argv[i],"%lf",&_epsilon_ );
	if ( status <= 0 ) _ErrorParse( "parsing -epsilon...\n", 0 );
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
      
      else if ( (strcmp ( argv[i], "-method" ) == 0) ||
		(strcmp ( argv[i], "-meth" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -method...\n", 0 );  
	if ( strcmp ( argv[i], "exact" ) == 0 ) {
	  typeMethod = _EXACT_;
	}
	else if ( strcmp ( argv[i], "quadr" ) == 0 
		  || strcmp ( argv[i], "quadratic" ) == 0 
		  || strcmp ( argv[i], "brute" ) == 0  
		  || strcmp ( argv[i], "brute-force" ) == 0 ) {
	  typeMethod = _QUADR_;
	} 
	else if ( strcmp ( argv[i], "apprx" ) == 0 ) {
	  typeMethod = _APPRX_;
	}
	else if ( strcmp ( argv[i], "approx" ) == 0 ) {
	  typeMethod = _APPRX_;
	}
	else if ( strcmp ( argv[i], "har-peled" ) == 0
		  || strcmp ( argv[i], "peled" ) == 0 ) {
	  set_init_peled_diameter_to_largest_dim();
	  set_update_peled_diameter_to_quadratic();
	  typeMethod = _PELED_;
	}
	else if ( strcmp ( argv[i], "hybrid1" ) == 0 ) {
	  set_init_peled_diameter_to_maximal_seg();
	  set_update_peled_diameter_to_quadratic();
	  typeMethod = _PELED_;
	}
	else if ( strcmp ( argv[i], "hybrid2" ) == 0 ) {
	  set_init_peled_diameter_to_maximal_seg();
	  set_update_peled_diameter_to_max_exact_with_diameter();
	  typeMethod = _PELED_;
	} 
	else  {
	  sprintf( modlname, "unknown method = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      }

      else if ( strcmp ( argv[i], "-Q-reduction" ) == 0 ) {
	_DoTryToReduceQ();
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
      sprintf( modlname, "%s", argv[i] );
      typeDistribution = IN_PTS_MODEL;
      if ( 0 ) {
	sprintf( modlname, "unknown option %s\n", argv[i] );
	_ErrorParse( modlname, 0);
      }
    }
  }


  
  _SetRandomSeed( seedRandom );

  fprintf( stderr, "\n" );
  fprintf( stderr, "generating / reading points ... " );

  switch( typeDistribution ) {
  default :
    listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					  _min_diameter_, _psommet_, 
					  typeDistribution, _dim_ );
    break;
  case IN_PTS_MODEL:
    listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
    break;
  case IN_PLY_MODEL:
    listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
    break;
  case IN_MODEL:
    listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
    break;
  }

  fprintf( stderr, "done\n" );

  if ( listOfPoints == NULL ) {
    _ErrorParse( "unable to generate/read list of points\n", 0 );
  }

  _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );

  if ( _nbpoints_ <= 0 ) {
    fprintf( stderr, "no points ?\n" );
    exit( 0 );
  }


#ifdef WIN32
  fprintf( stdout, "random seed = %u\n", seedRandom );
#else
  fprintf( stdout, "random seed = %ld\n", seedRandom );
#endif
  fprintf( stdout, "set = %d points of dimension %d from",
	   _nbpoints_, _dim_ );
  switch( typeDistribution ) {
  default :
  case IN_CUBE : 
    fprintf( stdout, " distribution in a cube\n" );   break;
  case ON_ELLIPSOID :
    fprintf( stdout, " distribution on an ellipsoid\n" );   break;
  case ON_REG_ELLIPSOID :
    fprintf( stdout, " distribution on a regular/gentle ellipsoid\n" );   break;
  case ON_SPHERE :
    fprintf( stdout, " distribution on a sphere\n" );   break;
  case IN_SPHERE :
    fprintf( stdout, " distribution in a sphere\n" );   break;
  case IN_CST_DIAMETER :
    fprintf( stdout, " distribution in a set of cst diameter\n" );
    fprintf( stdout, "(generated from a polygon with %d vertices)\n", 2*_psommet_+1 );
    
    break;
  case IN_PTS_MODEL:
  case IN_PLY_MODEL:
  case IN_MODEL:
    fprintf( stdout, " file '%s'\n", modlname );
  }


  c0 = clock();

  switch ( typeMethod ) {
  default :
  case _EXACT_ :
    (void)_ExactDiameterInOneList( &pair1, listOfPoints, 
				   0, _nbpoints_-1, _dim_ );
    break;
  case _PELED_ :
    (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					0, _nbpoints_-1, _dim_ );
    break;
  case _APPRX_ :
    upper = _EstimeDiameterInOneList( &pair1, listOfPoints, 
				      0, _nbpoints_-1, _dim_, _epsilon_ );
    break;
  case _QUADR_ :
    (void)_QuadraticDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
  }
  
  c1 = clock();
  
  fprintf( stdout, "\n" );
  switch ( typeMethod ) {
  default :
    _PrintSegment( stdout, &pair1, "Points realizing the diameter", _dim_ );
    break;
  case _APPRX_ : 
    _PrintSegment( stdout, &pair1, "Points approximating the diameter", _dim_ );
    break;
  }
  fprintf( stdout, "computation time = %g sec.\n",
	  (double)(c1-c0)/(double)CLOCKS_PER_SEC );
  fprintf( stdout, "\n" );

  free( listOfPoints );

 return( 0 );
}
