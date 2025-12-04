/***********************************************************************
LammpsTrace - Class for molecular dynamics simulation traces produced by
the LAMMPS MD simulation code.
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

#include "LammpsTrace.h"

#include <string>
#include <Misc/StdError.h>
#include <Misc/MessageLogger.h>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>

namespace MD {

namespace {

/****************
Helper functions:
****************/

void checkItem(IO::ValueSource& file,const char* item) // Checks an ITEM: line in a Lammp file
	{
	if(!file.isLiteral("ITEM:")||!file.isString(item)||!file.isLiteral('\n'))
		throw Misc::makeStdErr(0,"Missing %s item",item);
	}

void checkEndline(IO::ValueSource& file) // Checks whether the next literal is a line end
	{
	if(!file.isLiteral('\n'))
		throw std::runtime_error("Bad line");
	}

}

/****************************
Methods of class LammpsTrace:
****************************/

LammpsTrace::LammpsTrace(const char* metaDataFileName,const char* traceFileName)
	:traceFile(IO::openSeekableFile(traceFileName)),
	 bbox(Box::empty)
	{
	/* Parse the meta data file: */
	IO::ValueSource metaData(IO::openFile(metaDataFileName));
	metaData.setPunctuation('\n',true);
	metaData.skipWs();
	int state=0;
	Scalar scale(1);
	while(!metaData.eof())
		{
		if(metaData.peekc()=='\n'||metaData.peekc()=='#')
			{
			/* Skip an empty or comment line: */
			metaData.skipLine();
			metaData.skipWs();
			}
		else
			{
			if(state==0) // Read unit conversion factor
				{
				scale=Scalar(metaData.readNumber());
				state=1; // Read atom type definitions next
				}
			else // Read atom type definition
				{
				/* Read an atom type identifer: */
				Atom::TypeID typeId(metaData.readUnsignedInteger());
				
				/* Parse an atom type definition: */
				Atom::Type atomType;
				atomType.element=Atom::parseElement(metaData.readString().c_str());
				atomType.halfBondLength=Atom::getHalfBondLength(atomType.element);
				atomType.radius=Atom::getVanDerWaalsRadius(atomType.element);
				
				/* Check if there is a custom color definition: */
				if(metaData.peekc()!='\n')
					{
					/* Parse an RGB color: */
					for(int i=0;i<3;++i)
						atomType.color[i]=float(metaData.readNumber());
					}
				else
					{
					/* Use the element's standard color: */
					atomType.color=Atom::getColor(atomType.element);
					}
				
				/* Store the type definition: */
				setAtomType(typeId,atomType);
				}
			
			/* Skip the end-of-line: */
			if(!metaData.isLiteral('\n'))
				throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Missing line end in meta data file %s",metaDataFileName);
			}
		}
	
	/* Scan the trace file to find the starting positions of all included time steps: */
	IO::ValueSource trace(traceFile);
	trace.setPunctuation('\n',true);
	trace.skipWs();
	try
		{
		while(!trace.eof())
			{
			/*******************
			Parse the time step:
			*******************/
			
			TimeStep timeStep;
			
			/* Read the time step's time point: */
			checkItem(trace,"TIMESTEP");
			timeStep.timePoint=trace.readNumber();
			checkEndline(trace);
			
			/* Read the number of atoms: */
			checkItem(trace,"NUMBER OF ATOMS");
			timeStep.numAtoms=trace.readUnsignedInteger();
			checkEndline(trace);
			
			/* Read the bounding box: */
			checkItem(trace,"BOX BOUNDS pp pp pp");
			for(int i=0;i<3;++i)
				{
				timeStep.bbox.min[i]=Scalar(trace.readNumber())*scale;
				timeStep.bbox.max[i]=Scalar(trace.readNumber())*scale;
				checkEndline(trace);
				}
			
			/* Read the atoms: */
			checkItem(trace,"ATOMS id type xs ys zs");
			
			/* Save the file position of the time step's atom positions: */
			timeStep.atomsFilePos=traceFile->getReadPos()-1; // Subtract one to account for IO::ValueSource's one-character read-ahead
			
			for(size_t atomIndex=0;atomIndex<timeStep.numAtoms;++atomIndex)
				{
				/* Skip atom ID and type: */
				trace.readUnsignedInteger();
				trace.readUnsignedInteger();
				
				/* Skip atom position: */
				for(int i=0;i<3;++i)
					trace.readNumber();
				
				/* Check end-of-line: */
				checkEndline(trace);
				}
			
			/* Save the time step: */
			timeSteps.push_back(timeStep);
			
			/* Add the time step's bounding box to the overall box: */
			bbox.addBox(timeStep.bbox);
			}
		}
	catch(const std::runtime_error& err)
		{
		Misc::sourcedLogWarning(__PRETTY_FUNCTION__,"Stopped reading trace file %s after %u time steps due to exception %s",traceFileName,(unsigned int)(timeSteps.size()),err.what());
		}
	}

size_t LammpsTrace::getNumTimeSteps(void) const
	{
	return timeSteps.size();
	}

double LammpsTrace::getTimePoint(size_t timeStepIndex) const
	{
	return timeSteps[timeStepIndex].timePoint;
	}

size_t LammpsTrace::getNumAtoms(size_t timeStepIndex) const
	{
	return timeSteps[timeStepIndex].numAtoms;
	}

Box LammpsTrace::calcBoundingBox(void) const
	{
	return bbox;
	}

void LammpsTrace::loadAtoms(size_t timeStepIndex,Atom atoms[])
	{
	const TimeStep& timeStep=timeSteps[timeStepIndex];
	
	/* Seek to the beginning of the time step's atom positions: */
	traceFile->setReadPosAbs(timeStep.atomsFilePos);
	
	/* Create a value source to read data: */
	IO::ValueSource trace(traceFile);
	trace.setPunctuation('\n',true);
	trace.skipWs();
	
	/* Set up coordinate transformation; coordinates are normalized to [0, 1]: */
	Scalar scale[3],offset[3];
	for(int i=0;i<3;++i)
		{
		scale[i]=timeStep.bbox.max[i]-timeStep.bbox.min[i];
		offset[i]=timeStep.bbox.min[i];
		}
	
	/* Read the time step's atoms: */
	Atom* atomPtr=atoms;
	for(size_t atomIndex=0;atomIndex<timeStep.numAtoms;++atomIndex,++atomPtr)
		{
		/* Skip atom ID for now: */
		trace.readUnsignedInteger();
		
		/* Read atom type ID: */
		atomPtr->setTypeId(Atom::TypeID(trace.readUnsignedInteger()));
		
		/* Read atom position: */
		Point pos;
		for(int i=0;i<3;++i)
			pos[i]=Scalar(trace.readNumber())*scale[i]+offset[i];
		atomPtr->setPosition(pos);
		
		/* Skip end-of-line: */
		trace.skipString();
		}
	}

}
