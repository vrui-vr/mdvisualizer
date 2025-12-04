/***********************************************************************
BinaryTraceFile - Class to encapsulate operations on MD trace files in
compact binary format.
Copyright (c) 2005 Oliver Kreylos
***********************************************************************/

#ifndef MD_BINARYRACEFILE_INCLUDED
#define MD_BINARYRACEFILE_INCLUDED

#include <Misc/File.h>

#include "TraceFile.h"

namespace MD {

class BinaryTraceFile:public TraceFile
	{
	/* Elements: */
	private:
	Misc::File file; // The underlying file
	Box boundingBox; // Bounding box of all atoms in the trace file
	Misc::File::Offset atomsStart; // Offset of the first atom's data in the binary file
	
	/* Constructors and destructors: */
	public:
	BinaryTraceFile(const char* fileName); // Opens the binary trace file of the given name
	virtual ~BinaryTraceFile(void);
	
	/* Methods: */
	virtual Box getBoundingBox(void);
	virtual void readAtoms(Atom atoms[]);
	};

}

#endif
