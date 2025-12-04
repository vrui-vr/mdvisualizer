/***********************************************************************
BinaryTraceFile - Class to encapsulate operations on MD trace files in
compact binary format.
Copyright (c) 2005-2025 Oliver Kreylos

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
