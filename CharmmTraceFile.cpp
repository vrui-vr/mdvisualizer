/***********************************************************************
CharmmTraceFile - Class to encapsulate operations on MD trace files in
CHARMM format.
Copyright (c) 2005 Oliver Kreylos
***********************************************************************/

#include <stdio.h>
#include <string.h>

#include "CharmmTraceFile.h"

namespace MD {

/********************************
Methods of class CharmmTraceFile:
********************************/

void CharmmTraceFile::readAtomCache(void)
	{
	/* Allocate the atom cache: */
	atomCache=new Atom[numAtoms];
	
	/* Read the atoms from the file: */
	file.seekSet(atomsStart);
	char line[256];
	for(unsigned int i=0;i<numAtoms;++i)
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Parse the line: */
		char elementName[4];
		elementName[0]=line[16];
		for(int j=1;j<4;++j)
			elementName[j]='\0';
		Atom::Element type=Atom::parseType(elementName);
		Point pos;
		char coord[11];
		for(int j=0;j<3;++j)
			{
			strncpy(coord,line+j*10+20,10);
			coord[10]='\0';
			pos[j]=Scalar(atof(coord));
			}
		
		/* Save the atom: */
		atomCache[i]=Atom(type,pos);
		}
	}

CharmmTraceFile::CharmmTraceFile(const char* fileName)
	:file(fileName,"rt"),
	 atomCache(0)
	{
	/* Read the file header: */
	char line[256];
	do
		{
		file.gets(line,sizeof(line));
		}
	while(line[0]=='*');
	
	/* Extract the number of atoms: */
	numAtoms=atoi(line);
	
	/* Save the current file position: */
	atomsStart=file.tell();
	}

CharmmTraceFile::~CharmmTraceFile(void)
	{
	/* Delete any cached atoms: */
	delete[] atomCache;
	}

Box CharmmTraceFile::getBoundingBox(void)
	{
	/* Check if the atoms need to be cached: */
	if(atomCache==0)
		readAtomCache();
	
	/* Calculate the bounding box: */
	Box result=Box::empty;
	for(unsigned int i=0;i<numAtoms;++i)
		result.addPoint(atomCache[i].getPosition());
	return result;
	}

void CharmmTraceFile::readAtoms(Atom atoms[])
	{
	/* Check if the atoms need to be cached: */
	if(atomCache==0)
		readAtomCache();
	
	/* Copy the cached atoms: */
	for(unsigned int i=0;i<numAtoms;++i)
		atoms[i]=atomCache[i];
	}

}
