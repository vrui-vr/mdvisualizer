/***********************************************************************
Atom - Class describing an atom inside molecular dynamics trace data.
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

#ifndef MD_ATOM_INCLUDED
#define MD_ATOM_INCLUDED

#include <Misc/SizedTypes.h>
#include <Misc/StandardHashFunction.h>
#include <Misc/HashTable.h>
#include <Math/Math.h>

#include "MDGeometry.h"
#include "Color.h"

namespace MD {

class Atom
	{
	/* Embedded classes: */
	public:
	typedef unsigned int TypeID; // Type to identify an atom type used in an MD simulation
	
	enum Element // Enumerated type for elements (standard names)
		{
		H=0,He,
		Li,Be,B,C,N,O,F,Ne,
		Na,Mg,Al,Si,P,S,Cl,Ar,
		K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,
		Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,
		Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hq,Tl,Pb,Bi,Po,At,Rn,
		Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og,
		ElementsEnd
		};
	
	struct Type // Structure representing a type of atom used in an MD simulation
		{
		/* Elements: */
		public:
		Element element; // The atom type's chemical element
		Scalar halfBondLength; // Maximum distance from atom position to midpoint of a potential bond with another atom
		Scalar radius; // Radius of spheres used to visualize atoms of this type
		Color color; // Color used to visualize atoms of this type
		};
	
	typedef Misc::HashTable<TypeID,Type> TypeMap; // Type for hash tables mapping atom type IDs to atom types
	
	/* Static elements: */
	private:
	static const char elementNames[ElementsEnd][3]; // Table of standard element names
	static const Scalar elementRadii[ElementsEnd]; // Table of atomic radii (in Angstrom)
	static const Scalar elementCovalentRadii[ElementsEnd]; // Table of covalent atomic radii (in Angstrom)
	static const Scalar elementVanDerWaalsRadii[ElementsEnd]; // Table of van-der-Waals atomic radii (in Angstrom)
	static const Color elementColors[ElementsEnd]; // Table of standard colors used to visualize atoms
	
	/* Elements: */
	protected:
	TypeID typeId; // Identifier of the atom's type
	Point pos; // The atom's position (in Angstrom)
	
	/* Constructors and destructors: */
	public:
	Atom(void) // Default constructor
		:typeId(ElementsEnd)
		{
		}
	Atom(TypeID sTypeId,const Point& sPos)
		:typeId(sTypeId),pos(sPos)
		{
		}
	
	/* Static methods: */
	static const char* getElementName(Element element) // Returns the name of an element
		{
		return elementNames[element];
		}
	static Scalar getRadius(Element element) // Returns radius of an element
		{
		return elementRadii[element];
		}
	static Scalar getCovalentRadius(Element element) // Returns radius of an element
		{
		return elementCovalentRadii[element];
		}
	static Scalar getVanDerWaalsRadius(Element element) // Returns radius of an element
		{
		return elementVanDerWaalsRadii[element];
		}
	static Scalar getHalfBondLength(Element element) // Returns half bond length of an element
		{
		return elementCovalentRadii[element]*Scalar(1.25); // Grant a 25% overshoot over covalent radius
		}
	static const Color& getColor(Element element) // Returns standard color of an element
		{
		return elementColors[element];
		}
	static Element parseElement(const char* elementName); // Converts a standard 2-character element name into an element type
	static void createTypeMap(TypeMap& typeMap); // Creates a default type map where atom types are identified by element number
	
	/* Methods: */
	void setTypeId(TypeID newTypeId) // Sets an atom's type identifier
		{
		typeId=newTypeId;
		}
	TypeID getTypeId(void) const // Returns an atom's type identifier
		{
		return typeId;
		}
	const Type& getType(const TypeMap& typeMap) const // Returns an atom's type definition from the given type map
		{
		return typeMap.getEntry(typeId).getDest();
		}
	Element getElement(const TypeMap& typeMap) const // Returns an atom's element from the given type map
		{
		return typeMap.getEntry(typeId).getDest().element;
		}
	Scalar getHalfBondLength(const TypeMap& typeMap) const // Returns an atom's half bond length from the given type map
		{
		return typeMap.getEntry(typeId).getDest().halfBondLength;
		}
	Scalar getRadius(const TypeMap& typeMap) const // Returns an atom's visualization sphere radius from the given type map
		{
		return typeMap.getEntry(typeId).getDest().halfBondLength;
		}
	const Color& getColor(const TypeMap& typeMap) const // Returns an atom's visualization color from the given type map
		{
		return typeMap.getEntry(typeId).getDest().color;
		}
	const Point& getPosition(void) const // Returns the position of an atom
		{
		return pos;
		}
	void setPosition(const Point& newPosition) // Sets the position of an atom
		{
		pos=newPosition;
		}
	friend Scalar dist(const Atom& atom1,const Atom& atom2) // Returns the Euclidean distance between two atoms
		{
		return Scalar(Geometry::dist(atom1.pos,atom2.pos));
		}
	friend bool canBond(const Atom& atom1,const Atom& atom2,const TypeMap& typeMap) // Returns true if two atoms can form a bond based on the given type map
		{
		return Geometry::sqrDist(atom1.pos,atom2.pos)<=Math::sqr(atom1.getHalfBondLength(typeMap)+atom2.getHalfBondLength(typeMap));
		}
	template <class SourceParam>
	void read(SourceParam& source) // Reads an atom from a binary source
		{
		typeId=TypeID(source.template read<Misc::UInt32>());
		for(int i=0;i<3;++i)
			pos[i]=Scalar(source.template read<Misc::Float32>());
		}
	template <class SinkParam>
	void write(SinkParam& sink) const // Writes an atom to a binary sink
		{
		sink.template write<Misc::UInt32>(typeId);
		for(int i=0;i<3;++i)
			sink.template write(Misc::Float32(pos[i]));
		}
	};

}

#endif
