/***********************************************************************
CharmmTrace - Class for molecular dynamics simulation traces in CHARMM
format.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#include "CharmmTrace.h"

#include <ctype.h>
#include <stdio.h>
#include <algorithm>
#include <Misc/FileNameExtensions.h>
#include <Misc/PrintfTemplateTests.h>
#include <Misc/ThrowStdErr.h>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>

namespace MD {

/****************************
Methods of class CharmmTrace:
****************************/

CharmmTrace::CharmmTrace(const char* sTraceFileNameTemplate)
	:baseDirectory(IO::openFileDirectory(sTraceFileNameTemplate)),
	 traceFileNameTemplate(Misc::getFileName(sTraceFileNameTemplate)),
	 numAtoms(0),bbox(Box::empty)
	{
	/* Check if the file name template is valid: */
	unsigned int conversionStart,conversionLength;
	if(!Misc::isValidTemplate(traceFileNameTemplate,'u',1024,&conversionStart,&conversionLength))
		Misc::throwStdErr("CharmmTrace: %s is not a valid trace file name template",sTraceFileNameTemplate);
	
	/* Find all files in the base directory that match the file name template: */
	while(baseDirectory->readNextEntry())
		if(baseDirectory->getEntryType()==Misc::PATHTYPE_FILE)
			{
			/* Check if the template prefix matches the file name: */
			const char* fileName=baseDirectory->getEntryName();
			if(traceFileNameTemplate.compare(0,conversionStart,fileName,conversionStart)==0)
				{
				/* Extract the file's index: */
				unsigned int index=0;
				const char* cPtr=fileName+conversionStart;
				for(;isspace(*cPtr);++cPtr)
					;
				bool haveDigits=isdigit(*cPtr);
				for(;isdigit(*cPtr);++cPtr)
					index=index*10+(*cPtr-'0');
				if(haveDigits)
					{
					/* Check if the template suffix matches the file name: */
					if(traceFileNameTemplate.compare(conversionStart+conversionLength,std::string::npos,cPtr)==0)
						{
						/* Have a match; save the file number: */
						timeStepNumbers.push_back(index);
						}
					}
				}
			}
	
	/* Sort the list of time step file numbers: */
	std::sort(timeStepNumbers.begin(),timeStepNumbers.end());
	
	#if 0
	
	/* Create a standard type map: */
	Atom::createTypeMap(modifyAtomTypeMap());
	
	#endif
	
	/* Parse the first trace file to retrieve the number of atoms and calculate a bounding box: */
	char traceFileName[1024];
	snprintf(traceFileName,sizeof(traceFileName),traceFileNameTemplate.c_str(),timeStepNumbers.front());
	IO::ValueSource traceFile(baseDirectory->openFile(traceFileName));
	traceFile.setPunctuation('\n',true);
	traceFile.skipWs();
	
	/* Skip the file header: */
	while(traceFile.peekc()=='*')
		{
		traceFile.skipLine();
		traceFile.skipWs();
		}
	
	/* Read the number of atoms: */
	numAtoms=traceFile.readUnsignedInteger();
	traceFile.skipLine();
	traceFile.skipWs();
	
	/* Read all atoms and calculate the bounding box: */
	size_t i;
	for(i=0;i<numAtoms&&!traceFile.eof();++i)
		{
		/* Read the atom: */
		traceFile.readUnsignedInteger(); // Skip atom index
		traceFile.readUnsignedInteger(); // Skip residue index
		traceFile.skipString(); // Skip residue name
		
		/* Read the element name: */
		char elementName[4];
		elementName[0]=char(traceFile.getChar());
		for(int i=1;i<4;++i)
			elementName[i]='\0';
		traceFile.skipString();
		
		/* Create an atom type for this element: */
		Atom::Type atomType;
		atomType.element=Atom::parseElement(elementName);
		atomType.halfBondLength=Atom::getHalfBondLength(atomType.element);
		atomType.radius=Atom::getVanDerWaalsRadius(atomType.element);
		atomType.color=Atom::getColor(atomType.element);
		setAtomType(atomType.element,atomType);
		
		/* Read the atom position: */
		Point position;
		for(int j=0;j<3;++j)
			position[j]=Scalar(traceFile.readNumber());
		
		/* Skip the rest of the line: */
		traceFile.skipLine();
		traceFile.skipWs();
		
		/* Add the atom to the bounding box: */
		bbox.addPoint(position);
		}
	if(i<numAtoms)
		Misc::throwStdErr("CharmmTrace: Trace file %u is not a valid CHARMM trace file",timeStepNumbers.front());
	}

CharmmTrace::~CharmmTrace(void)
	{
	}

size_t CharmmTrace::getNumTimeSteps(void) const
	{
	return timeStepNumbers.size();
	}

double CharmmTrace::getTimePoint(size_t timeStepIndex) const
	{
	/* Assume that time step file numbers are time points in some arbitrary unit: */
	return double(timeStepNumbers[timeStepIndex]);
	}

size_t CharmmTrace::getNumAtoms(size_t timeStepIndex) const
	{
	return numAtoms;
	}

Box CharmmTrace::calcBoundingBox(void) const
	{
	return bbox;
	}

void CharmmTrace::loadAtoms(size_t timeStepIndex,Atom atoms[])
	{
	/* Open the trace file of the given index: */
	char traceFileName[1024];
	snprintf(traceFileName,sizeof(traceFileName),traceFileNameTemplate.c_str(),timeStepNumbers[timeStepIndex]);
	IO::ValueSource traceFile(baseDirectory->openFile(traceFileName));
	traceFile.setPunctuation('\n',true);
	traceFile.skipWs();
	
	/* Skip the file header: */
	while(traceFile.peekc()=='*')
		{
		traceFile.skipLine();
		traceFile.skipWs();
		}
	
	/* Read the number of atoms: */
	size_t fileNumAtoms=traceFile.readUnsignedInteger();
	if(fileNumAtoms!=numAtoms)
		Misc::throwStdErr("CharmmTrace: Trace file %u has mismatching number of atoms",timeStepNumbers[timeStepIndex]);
	traceFile.skipLine();
	traceFile.skipWs();
	
	/* Read all atoms: */
	size_t i;
	for(i=0;i<numAtoms&&!traceFile.eof();++i)
		{
		/* Read the atom: */
		traceFile.readUnsignedInteger(); // Skip atom index
		traceFile.readUnsignedInteger(); // Skip residue index
		traceFile.skipString(); // Skip residue name
		
		/* Read and parse the element name: */
		char elementName[4];
		elementName[0]=char(traceFile.getChar());
		for(int i=1;i<4;++i)
			elementName[i]='\0';
		traceFile.skipString();
		atoms[i].setTypeId(Atom::TypeID(Atom::parseElement(elementName)));
		
		/* Read the atom position: */
		Point pos;
		for(int j=0;j<3;++j)
			pos[j]=Scalar(traceFile.readNumber());
		atoms[i].setPosition(pos);
		
		/* Skip the rest of the line: */
		traceFile.skipLine();
		traceFile.skipWs();
		}
	if(i<numAtoms)
		Misc::throwStdErr("CharmmTrace: Trace file %u is not a valid CHARMM trace file",timeStepNumbers[timeStepIndex]);
	}

}
