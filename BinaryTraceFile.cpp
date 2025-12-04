/***********************************************************************
BinaryTraceFile - Class to encapsulate operations on MD trace files in
compact binary format.
Copyright (c) 2005 Oliver Kreylos
***********************************************************************/

#include "BinaryTraceFile.h"

namespace MD {

/********************************
Methods of class BinaryTraceFile:
********************************/

BinaryTraceFile::BinaryTraceFile(const char* fileName)
	:file(fileName,"rb",Misc::File::LittleEndian)
	{
	/* Read the trace file header: */
	numAtoms=file.read<unsigned int>();
	Point min,max;
	file.read(min.getComponents(),3);
	file.read(max.getComponents(),3);
	boundingBox=Box(min,max);
	atomsStart=file.tell();
	}

BinaryTraceFile::~BinaryTraceFile(void)
	{
	}

Box BinaryTraceFile::getBoundingBox(void)
	{
	return boundingBox;
	}

void BinaryTraceFile::readAtoms(Atom atoms[])
	{
	/* Read all atoms from the file: */
	for(unsigned int i=0;i<numAtoms;++i)
		{
		Atom::Element type=Atom::Element(file.read<unsigned char>());
		Point pos;
		file.read(pos.getComponents(),3);
		atoms[i]=Atom(type,pos);
		}
	}

}
