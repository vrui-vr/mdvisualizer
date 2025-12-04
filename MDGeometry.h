/***********************************************************************
MDGeometry - Definition of the Euclidean space underlying the visualized
molecular dynamics trace data.
Copyright (c) 2005 Oliver Kreylos
***********************************************************************/

#ifndef MD_GEOMETRY_INCLUDED
#define MD_GEOMETRY_INCLUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Box.h>

namespace MD {

typedef float Scalar; // Scalar type of MD space
typedef Geometry::Point<Scalar,3> Point; // Point type of MD space
typedef Geometry::Vector<Scalar,3> Vector; // Vector type of MD space
typedef Geometry::Box<Scalar,3> Box; // Type for axis-aligned boxes in MD space

}

#endif
