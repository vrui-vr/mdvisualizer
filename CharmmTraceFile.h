/***********************************************************************
CharmmTraceFile - Class to encapsulate operations on MD trace files in
CHARMM format.
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

#ifndef MD_CHARMMTRACEFILE_INCLUDED
#define MD_CHARMMTRACEFILE_INCLUDED

#include <Misc/File.h>

#include "TraceFile.h"

namespace MD {

class CharmmTraceFile:public TraceFile
	{
	/* Elements: */
	private:
	Misc::File file; // The underlying file
	Misc::File::Offset atomsStart; // Offset of the first atom's line in the CHARMM file
	Atom* atomCache; // Temporary buffer for all atoms in the file; since we have to read the entire file for the bounding box, we might as well cache the atoms
	
	/* Private methods: */
	void readAtomCache(void); // Reads atoms from the file into the cache
	
	/* Constructors and destructors: */
	public:
	CharmmTraceFile(const char* fileName); // Opens the CHARMM file of the given name
	virtual ~CharmmTraceFile(void);
	
	/* Methods: */
	virtual Box getBoundingBox(void);
	virtual void readAtoms(Atom atoms[]);
	};

}

#endif
