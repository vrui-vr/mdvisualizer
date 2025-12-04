/***********************************************************************
EspressoTrace - Class for molecular dynamics simulation traces produced
by the Espresso MD simulation code.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#include "EspressoTrace.h"

#include <string>
#include <Misc/ThrowStdErr.h>
#include <Misc/MessageLogger.h>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>

// DEBUGGING
#include <iostream>

namespace MD {

/******************************
Methods of class EspressoTrace:
******************************/

EspressoTrace::EspressoTrace(const char* traceFileName)
	:traceFile(IO::openSeekableFile(traceFileName)),
	 bbox(Box::empty)
	{
	/* Scan the trace file to find the starting positions of all included time steps: */
	IO::ValueSource trace(traceFile);
	trace.setPunctuation("=()\n");
	trace.skipWs();
	try
		{
		while(!trace.eof())
			{
			/*****************************************************************
			Find the beginning of the next time step, identified by the
			"Entering Dynamics" keyword:
			*****************************************************************/
			
			/* Find the next dynamics section: */
			while(!trace.eof())
				{
				if(trace.isLiteral("Entering")&&trace.isLiteral("Dynamics:"))
					break;
				trace.skipLine();
				trace.skipWs();
				}
			if(trace.eof())
				break;
			
			bool sectionValid=true;
			unsigned int iterationNumber=0;
			double time=0.0;
			
			/* Parse the iteration number and time step: */
			sectionValid=sectionValid&&trace.isLiteral("iteration")&&trace.isLiteral('=');
			if(sectionValid)
				iterationNumber=trace.readInteger();
			sectionValid=sectionValid&&trace.isLiteral('\n')&&trace.isLiteral("time")&&trace.isLiteral('=');
			if(sectionValid)
				time=trace.readNumber();
			if(sectionValid)
				{
				trace.skipLine();
				trace.skipWs();
				}
			if(sectionValid)
				{
				/* Find the start of the atom position section: */
				while(!trace.eof())
					{
					if(trace.isLiteral("ATOMIC_POSITIONS"))
						break;
					trace.skipLine();
					trace.skipWs();
					}
				}
			sectionValid=sectionValid&&!trace.eof();
			sectionValid=sectionValid&&trace.isLiteral('(');
			if(sectionValid)
				trace.skipString();
			sectionValid=sectionValid&&trace.isLiteral(')')&&trace.isLiteral('\n');
			
			/* Check if we found a time step: */
			if(sectionValid)
				{
				/* Save the file position of the time step's atom positions: */
				TimeStep ts;
				ts.timePoint=time;
				ts.numAtoms=0;
				ts.bbox=Box::empty;
				ts.atomsFilePos=traceFile->getReadPos()-1; // Subtract one to account for IO::ValueSource's one-character read-ahead
				
				/* Check if this is the first time step: */
				if(timeSteps.empty())
					{
					/* Read atoms until the next blank line or eof: */
					while(!trace.eof()&&trace.peekc()!='\n')
						{
						/* Read the element name: */
						char elementName[3];
						elementName[0]=char(trace.getChar());
						elementName[1]='\0';
						if(!trace.isWs(trace.peekc()))
							elementName[1]=char(trace.getChar());
						elementName[2]='\0';
						trace.skipString();
						
						/* Create an atom type for this element: */
						Atom::Type atomType;
						atomType.element=Atom::parseElement(elementName);
						atomType.halfBondLength=Atom::getHalfBondLength(atomType.element);
						atomType.radius=Atom::getVanDerWaalsRadius(atomType.element);
						atomType.color=Atom::getColor(atomType.element);
						setAtomType(atomType.element,atomType);
						
						/* Read the atom position: */
						Point pos;
						for(int i=0;i<3;++i)
							pos[i]=Scalar(trace.readNumber())*Scalar(7.592405);
						ts.bbox.addPoint(pos);
						++ts.numAtoms;
						
						trace.skipLine();
						trace.skipWs();
						}
					}
				else
					{
					/* Read atoms until the next blank line or eof: */
					while(!trace.eof()&&trace.peekc()!='\n')
						{
						/* Skip the atom type: */
						trace.skipString();
						
						/* Read the atom position: */
						Point pos;
						for(int i=0;i<3;++i)
							pos[i]=Scalar(trace.readNumber())*Scalar(7.592405);
						ts.bbox.addPoint(pos);
						++ts.numAtoms;
						
						trace.skipLine();
						trace.skipWs();
						}
					}
				
				/* Save the time step: */
				timeSteps.push_back(ts);
				
				/* Add the time step's bounding box to the overall box: */
				bbox.addBox(ts.bbox);
				}
			}
		}
	catch(const std::runtime_error& err)
		{
		Misc::formattedLogWarning("EspressoTrace::EspressoTrace: Stopped reading trace file %s after %u time steps due to exception %s",traceFileName,(unsigned int)(timeSteps.size()),err.what());
		}
	
	// DEBUGGING
	std::cout<<"Trace file "<<traceFileName<<" contains "<<timeSteps.size()<<" time steps"<<std::endl;
	std::cout<<"Time steps contain "<<timeSteps[0].numAtoms<<" atoms"<<std::endl;
	}

size_t EspressoTrace::getNumTimeSteps(void) const
	{
	return timeSteps.size();
	}

double EspressoTrace::getTimePoint(size_t timeStepIndex) const
	{
	return timeSteps[timeStepIndex].timePoint;
	}

size_t EspressoTrace::getNumAtoms(size_t timeStepIndex) const
	{
	return timeSteps[timeStepIndex].numAtoms;
	}

Box EspressoTrace::calcBoundingBox(void) const
	{
	return bbox;
	}

void EspressoTrace::loadAtoms(size_t timeStepIndex,Atom atoms[])
	{
	const TimeStep& timeStep=timeSteps[timeStepIndex];
	
	/* Seek to the beginning of the time step's atom positions: */
	traceFile->setReadPosAbs(timeStep.atomsFilePos);
	
	/* Create a value source to read data: */
	IO::ValueSource trace(traceFile);
	trace.setPunctuation('\n',true);
	trace.skipWs();
	
	/* Read the time step's atoms: */
	Atom* atomPtr=atoms;
	for(size_t atomIndex=0;atomIndex<timeStep.numAtoms;++atomIndex,++atomPtr)
		{
		/* Read the element name: */
		char elementName[3];
		elementName[0]=char(trace.getChar());
		elementName[1]='\0';
		if(!trace.isWs(trace.peekc()))
			elementName[1]=char(trace.getChar());
		elementName[2]='\0';
		trace.skipString();
		atomPtr->setTypeId(Atom::TypeID(Atom::parseElement(elementName)));
		
		/* Read atom position: */
		Point pos;
		for(int i=0;i<3;++i)
			pos[i]=Scalar(trace.readNumber())*Scalar(7.592405);
		atomPtr->setPosition(pos);
		
		/* Skip end-of-line: */
		trace.skipLine();
		trace.skipWs();
		}
	}

}
