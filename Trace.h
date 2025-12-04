/***********************************************************************
Trace - Base class for molecular dynamics simulation traces, i.e., time
series of atom positions.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#ifndef MD_TRACE_INCLUDED
#define MD_TRACE_INCLUDED

#include <stddef.h>

#include "MDGeometry.h"
#include "Atom.h"

namespace MD {

class Trace
	{
	/* Elements: */
	private:
	Atom::TypeMap atomTypeMap; // Map from atom type IDs to atom types
	
	/* Protected methods: */
	protected:
	Atom::TypeMap& modifyAtomTypeMap(void) // Returns the atom type map for modification
		{
		return atomTypeMap;
		}
	void setAtomType(Atom::TypeID atomTypeId,const Atom::Type& atomType) // Adds a mapping from an atom type ID to an atom type
		{
		/* Put the association into the hash table: */
		atomTypeMap.setEntry(Atom::TypeMap::Entry(atomTypeId,atomType));
		}
	
	/* Constructors and destructors: */
	public:
	Trace(void);
	virtual ~Trace(void);
	
	/* Methods: */
	const Atom::TypeMap& getAtomTypeMap(void) const // Returns the atom type map
		{
		return atomTypeMap;
		}
	Scalar calcMaxHalfBondLength(void) const; // Returns the maximum half bond length of all atom types
	const Atom::Type& getAtomType(const Atom& atom) const // Returns the type definition of the given atom
		{
		return atom.getType(atomTypeMap);
		}
	virtual size_t getNumTimeSteps(void) const =0; // Returns the number of time steps in the trace
	virtual double getTimePoint(size_t timeStepIndex) const =0; // Returns the time point of the given time step in some arbitrary unit
	virtual size_t getNumAtoms(size_t timeStepIndex) const =0; // Returns the number of atoms in the given time step
	virtual Box calcBoundingBox(void) const =0; // Calculates the bounding box of all atoms in all time steps
	virtual void loadAtoms(size_t timeStepIndex,Atom atoms[]) =0; // Loads the given time step's atoms into the provided array
	};

}

#endif
