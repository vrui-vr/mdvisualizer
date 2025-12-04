/***********************************************************************
MDGeometry - Definition of the Euclidean space underlying the visualized
molecular dynamics trace data.
Copyright (c) 2005-2025 Oliver Kreylos

This file is part of the MD Visualizer (MDVisualizer).

The MD Visualizer is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

The MD Visualizer is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with the MD Visualizer; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
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
