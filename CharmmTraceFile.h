/***********************************************************************
CharmmTraceFile - Class to encapsulate operations on MD trace files in
CHARMM format.
Copyright (c) 2005 Oliver Kreylos
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
