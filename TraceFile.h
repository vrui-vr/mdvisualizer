/***********************************************************************
TraceFile - Base class to encapsulate operations on MD trace files.
Copyright (c) 2005-2019 Oliver Kreylos
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
