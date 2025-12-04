/***********************************************************************
Trace - Base class for molecular dynamics simulation traces, i.e., time
series of atom positions.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#include "Trace.h"

namespace MD {

/**********************
Methods of class Trace:
**********************/

Trace::Trace(void)
	:atomTypeMap(151)
	{
	}

Trace::~Trace(void)
	{
	}

Scalar Trace::calcMaxHalfBondLength(void) const
	{
	/* Find the maximum half bond length among all defined atom types: */
	Scalar result(0);
	for(Atom::TypeMap::ConstIterator atmIt=atomTypeMap.begin();!atmIt.isFinished();++atmIt)
		if(result<atmIt->getDest().halfBondLength)
			result=atmIt->getDest().halfBondLength;
	
	return result;
	}

}
