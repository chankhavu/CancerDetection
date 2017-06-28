/*************************************************************************
 *  -
 *
 * $Id: pick.h,v 1.1 2002/07/24 12:43:44 greg Exp $
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




#ifndef _pick_h_
#define _pick_h_







#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <alloc.h>
#include <rand.h>


typedef enum {
  IN_CUBE=0,
  ON_ELLIPSOID=1,
  ON_REG_ELLIPSOID=2,
  IN_SPHERE=3,
  IN_CST_DIAMETER=4,
  ON_SPHERE=5,
  IN_PLY_MODEL=6,
  IN_MODEL=7,
  IN_PTS_MODEL=8
} enumDistribution;  


extern void _VerboseInPick();
extern void _NoVerboseInPick();
extern void _DoApplyARandomRotation();
extern void _DoNotApplyARandomRotation();

extern void _ApplyARandomRotation( double **listOfPoints,
			    const int nbpoints, const int dim );


extern void * _PickPoints( const int nbpoints,
			 const double diameterMax,
			 const double diameterMin,
			 const int psommet,
			 enumDistribution typeDistribution,
			 const int dim );




extern void _PickPointsInUnitCube( double **listOfPoints,
				   const int dim,
				   const int nbpoints,
				   const double diameter );

extern void _PickPointsOnUnitSphere( double **listOfPoints,
				     const int dim,
				     const int nbpoints,
				     const double diameter );

extern void _PickPointsInUnitSphere( double **listOfPoints,
				     const int dim,
				     const int nbpoints,
				     const double diameter );


extern void _PickPointsOnEllipsoidWithRadii( double **listOfPoints,
					     double *radii,
					     const int dim,
					     const int nbpoints );

extern void _PickPointsOnEllipsoid( double **listOfPoints,
				    const int dim,
				    const int nbpoints,
				    const double diameterMax,
				    const double diameterMin );

extern void _PickPointsOnRegularEllipsoid( double **listOfPoints,
					   const int dim,
					   const int nbpoints,
					   const double diameterMax,
					   const double diameterMin );


extern void _Pick2PointsInCstDiameter( double **listOfPoints,
				       const int dim,
				       const int nbpoints,
				       const double diameter,
				       const int psommet );



extern int _DivideListOfPoints( double **listOfPoints,
				const int nbpoints, const int dim );






#endif
