/***********************************************************************
CharmmTrace - Class for molecular dynamics simulation traces in CHARMM
format.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#ifndef MD_CHARMMTRACE_INCLUDED
#define MD_CHARMMTRACE_INCLUDED

#include <string>
#include <vector>
#include <IO/Directory.h>

#include "Trace.h"

namespace MD {

class CharmmTrace:public Trace
	{
	/* Elements: */
	private:
	IO::DirectoryPtr baseDirectory; // Directory containing all trace files
	std::string traceFileNameTemplate; // Template for trace file names relative to base directory
	std::vector<unsigned int> timeStepNumbers; // Sorted list of trace file time step numbers
	size_t numAtoms; // Number of atoms in each time step; taken from first time step
	Box bbox; // Bounding box containing all atoms; estimated from first time step
	
	/* Constructors and destructors: */
	public:
	CharmmTrace(const char* sTraceFileNameTemplate); // Reads all trace files fitting the given template
	virtual ~CharmmTrace(void);
	
	/* Methods from class Trace: */
	virtual size_t getNumTimeSteps(void) const;
	virtual double getTimePoint(size_t timeStepIndex) const;
	virtual size_t getNumAtoms(size_t timeStepIndex) const;
	virtual Box calcBoundingBox(void) const;
	virtual void loadAtoms(size_t timeStepIndex,Atom atoms[]);
	};

}

#endif
