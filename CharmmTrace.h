/***********************************************************************
CharmmTrace - Class for molecular dynamics simulation traces in CHARMM
format.
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
