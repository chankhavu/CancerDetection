/*************************************************************************
 *  -
 *
 * $Id: pick.c,v 1.1 2002/07/24 12:43:44 greg Exp $
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




#include <pick.h>


static int _do_apply_a_random_rotation_ = 0;
static int _verbose_ = 0;

static void RotationMatrixFromRotationVector( double *mat,
					      double *rot );


void _VerboseInPick()
{
  _verbose_ = 1;
}

void _NoVerboseInPick()
{
  _verbose_ = 0;
}


void _DoApplyARandomRotation()
{
  _do_apply_a_random_rotation_ = 1;
}

void _DoNotApplyARandomRotation()
{
  _do_apply_a_random_rotation_ = 0;
}







void *_PickPoints( const int nbpoints,
		   const double diameterMax,
		   const double diameterMin,
		   const int psommet,
		   enumDistribution typeDistribution,
		   const int dim )
{
  double **listOfPoints;

  listOfPoints = (double**)_AllocateListOfPoints( nbpoints, dim );
  if ( listOfPoints == NULL ) return( NULL );

  switch( typeDistribution ) {
  default :
  case IN_CUBE :
    _PickPointsInUnitCube( listOfPoints, dim, nbpoints, diameterMax );
    break;
  case ON_SPHERE :
    _PickPointsOnUnitSphere( listOfPoints, dim, nbpoints, diameterMax );
    break;
  case IN_SPHERE :
    _PickPointsInUnitSphere( listOfPoints, dim, nbpoints, diameterMax );
    break;
  case ON_ELLIPSOID :
    _PickPointsOnEllipsoid( listOfPoints, dim, nbpoints, 
			    diameterMax, diameterMin );
    break;
  case ON_REG_ELLIPSOID :
    _PickPointsOnRegularEllipsoid( listOfPoints, dim, nbpoints, 
				   diameterMax, diameterMin );
    break;
  case IN_CST_DIAMETER :
    _Pick2PointsInCstDiameter( listOfPoints, dim, nbpoints, diameterMax, psommet );
  }

  

  return( listOfPoints );
}


void _ApplyARandomRotation( double **listOfPoints,
			    const int nbpoints, const int dim )
{
  double mat[9], v[3], n, c, s, angle;
  int i;

  if ( _do_apply_a_random_rotation_  && ( dim == 2 || dim == 3 ) ) {

    angle = (2.0 * _GetRandomDoubleNb() - 1.0 ) *  3.1415926536;

    if ( _verbose_ ) {
      fprintf( stderr, " angle = %f rad = %f deg \n", angle, angle * 180.0 /  3.1415926536 );
    }

    if ( dim == 2 ) {
      c = cos( angle );
      s = sin( angle );
      
      for ( i=0; i<nbpoints; i++ ) {
	v[0] = c * listOfPoints[i][0] - s * listOfPoints[i][1];
	v[1] = s * listOfPoints[i][0] + c * listOfPoints[i][1];
	listOfPoints[i][0] = v[0];
	listOfPoints[i][1] = v[1];
      }
    }

    if ( dim == 3 ) {
      do { 
	v[0] = (2.0 * _GetRandomDoubleNb() - 1.0 );
	v[1] = (2.0 * _GetRandomDoubleNb() - 1.0 );
	v[2] = (2.0 * _GetRandomDoubleNb() - 1.0 );
	n = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
      } while ( n < 0.0001 );

      if ( _verbose_ ) {
	fprintf( stderr, " axe   = [ %f %f %f ]\n", v[0]/n, v[1]/n, v[2]/n );
      }

      v[0] *= angle / n;
      v[1] *= angle / n;
      v[2] *= angle / n;
      RotationMatrixFromRotationVector( mat, v );
      for ( i=0; i<nbpoints; i++ ) {
	v[0] = mat[0] * listOfPoints[i][0] + mat[1] * listOfPoints[i][1] + mat[2] * listOfPoints[i][2];
	v[1] = mat[3] * listOfPoints[i][0] + mat[4] * listOfPoints[i][1] + mat[5] * listOfPoints[i][2];
	v[2] = mat[6] * listOfPoints[i][0] + mat[7] * listOfPoints[i][1] + mat[8] * listOfPoints[i][2];
	listOfPoints[i][0] = v[0];
	listOfPoints[i][1] = v[1];
	listOfPoints[i][2] = v[2];
      }
    }

  }
}




void _PickPointsInUnitCube( double **listOfPoints,
			    const int dim,
			    const int nbpoints,
			    const double diameter )
{
  int i, j;
  double arete = diameter / sqrt( (double)dim );
  
  for ( i=0; i<nbpoints; i++ ) {
    for ( j=0; j<dim; j++ )
      listOfPoints[i][j] = arete * (2.0 * _GetRandomDoubleNb() - 1.0) / 2.0;
  }
}





void _PickPointsOnUnitSphere( double **listOfPoints,
			      const int dim,
			      const int nbpoints,
			      const double diameter )
{
  int i, j;
  double t;
  double *pt;
  double rayon = diameter / 2.0;
  
  for ( i=0; i<nbpoints; i++ ) {
    pt = listOfPoints[i];
    do {
      
      j = 0;
      t = 0;
      do {
	pt[j] = 2.0 * _GetRandomDoubleNb() - 1.0;
	t += pt[j] * pt[j];
	if ( t > 1 ) {
	  j = 0;
	  t = 0;
	}
	else         j ++;
      } while( j < dim );
      t = sqrt( t );
    } while( t > 1 || t < 1e-10 );

    for ( j=0; j<dim; j++ ) 
      listOfPoints[i][j] *= rayon / t;
  }
}







void _PickPointsInUnitSphere( double **listOfPoints,
			      const int dim,
			      const int nbpoints,
			      const double diameter )
{
  int i, j;
  double t;
  double *pt;
  double rayon = diameter / 2.0;
  
  for ( i=0; i<nbpoints; i++ ) {
    pt = listOfPoints[i];
    do {
      
      j = 0;
      t = 0;
      do {
	pt[j] = 2.0 * _GetRandomDoubleNb() - 1.0;
	t += pt[j] * pt[j];
	if ( t > 1 ) {
	  j = 0;
	  t = 0;
	}
	else         j ++;
      } while( j < dim );
      t = sqrt( t );
    } while( t > 1 );

    for ( j=0; j<dim; j++ ) 
      listOfPoints[i][j] *= rayon;
  }
}



void _PickPointsOnEllipsoidWithRadii( double **listOfPoints,
					     double *radii,
					     const int dim,
					     const int nbpoints )
{
  int i, j;
  double t;
  double *pt;
  
  for ( i=0; i<nbpoints; i++ ) {
    pt = listOfPoints[i];
    do {
      
      j = 0;
      t = 0;
      do {
	pt[j] = 2.0 * _GetRandomDoubleNb() - 1.0;
	t += pt[j] * pt[j];
	if ( t > 1 ) {
	  j = 0;
	  t = 0;
	}
	else         j ++;
      } while( j < dim );
      t = sqrt( t );
    } while( t > 1 || t < 1e-10 );

    t = 0;
    for ( j=0; j<dim; j++ ) 
      t += pt[j] * pt[j] / ( radii[j] * radii[j] );
    t = sqrt( t );
    
    for ( j=0; j<dim; j++ ) 
      pt[j] /= t;
  }
}






void _PickPointsOnEllipsoid( double **listOfPoints,
			     const int dim,
			     const int nbpoints,
			     const double diameterMax,
			     const double diameterMin )
{
  char *proc = "_PickPointsOnEllipsoid";
  double *radii = (double*)NULL;
  double rmax = diameterMax / 2.0;
  double rmin = diameterMin / 2.0;
  int j;
  
  radii = (double*)malloc( dim * sizeof(double) );
  if ( radii == (double*)NULL ) {
    if ( _verbose_ )
    fprintf( stderr, " %s: can not allocate radii.\n", proc );
    return;
  }

  radii[0] = rmax;
  for ( j=1; j<dim; j++ )
    radii[j] = rmin + (rmax - rmin) * _GetRandomDoubleNb();
  
  _PickPointsOnEllipsoidWithRadii( listOfPoints, radii, dim, nbpoints );
  
  free( radii );
  
}








void _PickPointsOnRegularEllipsoid( double **listOfPoints,
				    const int dim,
				    const int nbpoints,
				    const double diameterMax,
				    const double diameterMin )
{
  char *proc = "_PickPointsOnRegularEllipsoid";
  double *radii = (double*)NULL;
  double rmax = diameterMax / 2.0;
  double rmin = diameterMin / 2.0;
  double intervalLength = (rmax - rmin)/((double)(dim - 1));
  int j;
  
  radii = (double*)malloc( dim * sizeof(double) );
  if ( radii == (double*)NULL ) {
    if ( _verbose_ )
    fprintf( stderr, " %s: can not allocate radii.\n", proc );
    return;
  }
  
  radii[0] = rmax;
  radii[dim-1] = rmin;
  for ( j=1; j<dim-1; j++ ) {
    radii[j] = rmax - j*intervalLength;
    radii[j] += (2.0 * _GetRandomDoubleNb() - 1.0) * intervalLength/2.0;
  }

  _PickPointsOnEllipsoidWithRadii( listOfPoints, radii, dim, nbpoints );
  
  free( radii );
}











static double _pi_ = 3.141592653589793238462643383279;


/* distribution dans un ensemble a diametre constant
   bati sur un polygone a (2 * psommet + 1) sommets
*/
void _Pick2PointsInCstDiameter( double **listOfPoints,
				const int dim,
				const int nbpoints,
				const double diameter,
				const int psommet )
{
  char *proc = "_Pick2PointsInCstDiameter";
  int i;
  double r;
  double x, y, t;
  double xo, yo;
  double alpha;
  int j;
  int p = psommet;


  if ( dim != 2 ) {
    if ( _verbose_ )
    fprintf( stderr, " %s: can not deal with 'dim' != 2 \n", proc );
    exit( 0 );
  }

  if ( p < 1 ) p = 1;


  /* rayon englobant la structure
     1ere ligne => diametre = 1;
   */
  r = sqrt( 1.0 / (2.0 * (1.0 - cos( (2.0 * _pi_ * (double)p )/(2.0 * (double)p + 1.0) ))) );

  for ( i=0; i<nbpoints; i++ ) {
    do {
      x = r * (2.0 * _GetRandomDoubleNb() - 1.0);
      y = r * (2.0 * _GetRandomDoubleNb() - 1.0);
      t = sqrt( x*x + y*y );
      alpha = acos( x/t );
      if ( asin( y/t ) < 0 ) alpha = 2.0 * _pi_ - alpha;
      /* le point est
	 entre j*(2p+1)/(2*PI) et (j+1)*(2p+1)/(2*PI)
	 -> le point 'oppose' est donc j - p modulo (2p+1)
      */
      
      j = (int)( (2.0 * (double)p + 1.0)*alpha / (2.0 * _pi_) );
      j += p+1;
      
      xo = r * cos( (j * 2.0 * _pi_) / (2.0 * (double)p + 1.0) );
      yo = r * sin( (j * 2.0 * _pi_) / (2.0 * (double)p + 1.0) );
      t = (x-xo)*(x-xo) + (y-yo)*(y-yo);
    } while ( t > 1 );
    /* jusqu'ici on a suppose un diametre egal a 1
     */
    listOfPoints[i][0] = x * diameter;
    listOfPoints[i][1] = y * diameter;
  }
}






















int _DivideListOfPoints( double **listOfPoints,
			 const int nbpoints, const int dim )
{
  char *proc = "_DivideListOfPoints";
  double t;
  int i, j;
  int l;
  double *v = NULL;
  double *tmp;
  
  /* on tire un vecteur unitaire au hasard
   */
  v = (double*)malloc( dim * sizeof( double ) );
  if ( v == (double*)NULL ) {
    if ( _verbose_ )
    fprintf( stderr, " %s: can not allocate vector.\n", proc );
    return( -1 );
  }

  do {
    for ( t=0.0, j=0; j<dim; j++ ) {
      v[j] = 2.0 * _GetRandomDoubleNb() - 1.0;
      t += v[j] * v[j];
    }
    t = sqrt( t );
  } while( t > 1 || t < 1e-6 );
  
  for ( j=0; j<dim; j++ ) 
    v[j] /= t;
  

  
  /* on separe la liste en deux
     selon le produit scalaire
  */
  l = nbpoints-1;
  for ( i=0; i<=l; i++ ) {
  
    for ( t=0.0, j=0; j<dim; j++ ) {
      t += v[j] * listOfPoints[i][j];
    }

    if ( t <= 0 ) {
      /* on met en fin de liste
	 => l
      */
      tmp             = listOfPoints[i];
      listOfPoints[i] = listOfPoints[l];
      listOfPoints[l] = tmp;
      l --;
      i --;
    }
  }
  
  /* de 0 a l => t >  0
     de l+1 a nbpoints-1   => t <= 0

     l est l'indice du dernier point tel que t>0
  */
  
  free( v );
  
  return( l );
}












/* calcul la matrice de rotation a partir
   d'un vecteur de rotation r

   R = I + f(theta) X(r) + g(theta) [X(r) * X(r)]

   avec theta qui est l'angle de la rotation (norme du vecteur rotation), 
        g(theta) = (1 - cos(theta)) / (theta * theta)
        f(theta) = sin(theta) / theta
   et X(r) qui est la matrice du produit vectoriel par r.

   Si la rotation est donnee par un vecteur unitaire v et un angle
   theta, on a alors
   
   R = I + sin(theta) X(v) + (1 - cos(theta)) X(v)*X(v)

   D'ou on deduit que 
   trace(R) = 3 - 2*(1 - cos(theta)) = 1 + 2*cos(theta)


*/
void RotationMatrixFromRotationVector( double *mat,
				       double *rot )
{
  double f, g, theta, t2;
  
  t2 = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
  theta = sqrt( t2 );
  
  if ( theta > 1e-8 ) {
    f = sin( theta ) / theta;
    g = ( 1.0 - cos( theta ) ) / ( t2 );
    
    mat[0] = 1.0 - g * (rot[1]*rot[1] + rot[2]*rot[2]);
    mat[4] = 1.0 - g * (rot[2]*rot[2] + rot[0]*rot[0]);
    mat[8] = 1.0 - g * (rot[0]*rot[0] + rot[1]*rot[1]);
    
    mat[3] = mat[1] = g * rot[0] * rot[1];
    mat[6] = mat[2] = g * rot[0] * rot[2];
    mat[7] = mat[5] = g * rot[2] * rot[1];
    
    mat[1] -= f * rot[2];

    mat[2] += f * rot[1];
    mat[5] -= f * rot[0];
    
    mat[3] += f * rot[2];
    mat[6] -= f * rot[1];
    mat[7] += f * rot[0];
  }
  else {
    mat[0] = mat[4] = mat[8] = 1.0;
    mat[1] = mat[2] = mat[3] = 0.0;
    mat[5] = mat[6] = mat[7] = 0.0;
  }
}



