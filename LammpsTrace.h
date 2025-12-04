/***********************************************************************
LammpsTrace - Class for molecular dynamics simulation traces produced by
the LAMMPS MD simulation code.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#ifndef MD_LAMMPSTRACE_INCLUDED
#define MD_LAMMPSTRACE_INCLUDED

#include <vector>
#include <IO/SeekableFile.h>

#include "MDGeometry.h"
#include "Trace.h"

namespace MD {

class LammpsTrace:public Trace
	{
	/* Embedded classes: */
	private:
	struct TimeStep // Structure representing a single time step in a trace file
		{
		/* Elements: */
		public:
		double timePoint; // Time step's time point in some arbitrary unit
		size_t numAtoms; // Time step's number of atoms
		Box bbox; // Time step's bounding box
		IO::SeekableFile::Offset atomsFilePos; // Position in trace file where time step's atom positions begin
		};
	
	/* Elements: */
	private:
	IO::SeekableFilePtr traceFile; // Handle to the trace file
	std::vector<TimeStep> timeSteps; // List of time steps in the trace file
	Box bbox; // Overall bounding box of all time steps in the trace file
	
	/* Constructors and destructors: */
	public:
	LammpsTrace(const char* metaDataFileName,const char* traceFileName); // Opens the Lammp trace file of the given name using the meta data file of the given name
	
	/* Methods from class Trace: */
	virtual size_t getNumTimeSteps(void) const;
	virtual double getTimePoint(size_t timeStepIndex) const;
	virtual size_t getNumAtoms(size_t timeStepIndex) const;
	virtual Box calcBoundingBox(void) const;
	virtual void loadAtoms(size_t timeStepIndex,Atom atoms[]);
	};

}

#endif
