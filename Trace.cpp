/***********************************************************************
Trace - Base class for molecular dynamics simulation traces, i.e., time
series of atom positions.
Copyright (c) 2019-2025 Oliver Kreylos

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

#include "Trace.h"

namespace MD {

/**********************
Methods of class Trace:
**********************/

Trace::Trace(void)
	:atomTypeMap(151)
	{
	}

Trace::~Trace(void)
	{
	}

Scalar Trace::calcMaxHalfBondLength(void) const
	{
	/* Find the maximum half bond length among all defined atom types: */
	Scalar result(0);
	for(Atom::TypeMap::ConstIterator atmIt=atomTypeMap.begin();!atmIt.isFinished();++atmIt)
		if(result<atmIt->getDest().halfBondLength)
			result=atmIt->getDest().halfBondLength;
	
	return result;
	}

}
