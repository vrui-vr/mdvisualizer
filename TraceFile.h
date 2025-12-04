/***********************************************************************
TraceFile - Base class to encapsulate operations on MD trace files.
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

#ifndef MD_TRACEFILE_INCLUDED
#define MD_TRACEFILE_INCLUDED

#include "MDGeometry.h"
#include "Atom.h"

namespace MD {

class TraceFile
	{
	/* Elements: */
	protected:
	unsigned int numAtoms; // Number of atoms in the trace file
	
	/* Constructors and destructors: */
	public:
	TraceFile(void);
	virtual ~TraceFile(void);
	
	/* Methods: */
	unsigned int getNumAtoms(void) const // Returns the number of atoms in the trace file
		{
		return numAtoms;
		};
	virtual Box getBoundingBox(void) =0; // Returns the bounding box of all atoms in the trace file
	virtual void readAtoms(Atom atoms[]) =0; // Reads all atoms in the trace file into an array
	};

}

#endif
