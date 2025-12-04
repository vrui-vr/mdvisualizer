/***********************************************************************
Atom - Class describing an atom inside molecular dynamics trace data.
Copyright (c) 2005-2019 Oliver Kreylos
***********************************************************************/

#include "Atom.h"

#include <string.h>
#include <Misc/ThrowStdErr.h>

namespace MD {

/*****************************
Static elements of class Atom:
*****************************/

const char Atom::elementNames[Atom::ElementsEnd][3]=
	{"H","He",
	"Li","Be","B","C","N","O","F","Ne",
	"Na","Mg","Al","Si","P","S","Cl","Ar",
	"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
	"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
	"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
	"Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};

const Scalar Atom::elementRadii[Atom::ElementsEnd]=
	{
	/* Row 1: Elements H-He: */
	0.53,0.31,
	
	/* Row 2: Elements Li-Ne: */
	1.67,1.12,0.87,0.67,0.56,0.48,0.42,0.38,
	
	/* Row 3: Elements Na-Ar: */
	1.90,1.45,1.18,1.11,0.98,0.88,0.79,0.71,
	
	/* Row 4: Elements K-Kr: */
	2.43,1.94,1.84,1.76,1.71,1.66,1.61,1.56,1.52,1.49,1.45,1.42,1.36,1.25,1.14,1.03,0.94,0.88,
	
	/* Row 5: Elements Rb-Xe: */
	2.65,2.19,2.12,2.06,1.98,1.90,1.83,1.78,1.73,1.69,1.65,1.61,1.56,1.45,1.33,1.23,1.15,1.08,
	
	/* Row 6: Elements Cs-Rn: */
	2.98,2.53,0.00,0.00,2.47,2.06,2.05,2.38,2.31,2.33,2.25,2.28,2.26,2.26,2.22,2.22,2.17,2.08,2.00,1.93,1.88,1.85,1.80,1.77,1.74,1.71,1.56,1.54,1.43,1.35,1.27,1.20,
	
	/* Row 7: Elements Fr-Og: */
	0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
	};

const Scalar Atom::elementCovalentRadii[Atom::ElementsEnd]=
	{
	/* Row 1: Elements H-He: */
	0.37,0.32,
	
	/* Row 2: Elements Li-Ne: */
	1.34,0.90,0.82,0.77,0.75,0.73,0.71,0.69,
	
	/* Row 3: Elements Na-Ar: */
	1.54,1.30,1.18,1.11,1.06,1.02,0.99,0.97,
	
	/* Row 4: Elements K-Kr: */
	1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,1.10,
	
	/* Row 5: Elements Rb-Xe: */
	2.11,1.92,1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.30,
	
	/* Row 6: Elements Cs-Rn: */
	2.25,1.98,1.69,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.60,1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49,1.48,1.47,1.46,0.00,0.00,1.45,
	
	/* Row 7: Elements Fr-Og: */
	0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
	};

const Scalar Atom::elementVanDerWaalsRadii[Atom::ElementsEnd]=
	{
	/* Row 1: Elements H-He: */
	1.20,0.32,
	
	/* Row 2: Elements Li-Ne: */
	1.34,0.90,0.82,1.70,1.55,1.52,0.71,0.69,
	
	/* Row 3: Elements Na-Ar: */
	1.54,1.30,1.18,1.11,1.80,1.80,1.75,0.97,
	
	/* Row 4: Elements K-Kr: */
	1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.39,1.26,1.22,1.19,1.16,1.14,1.10,
	
	/* Row 5: Elements Rb-Xe: */
	2.11,1.92,1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.30,
	
	/* Row 6: Elements Cs-Rn: */
	2.25,1.98,1.69,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.60,1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49,1.48,1.47,1.46,0.00,0.00,1.45,
	
	/* Row 7: Elements Fr-Og: */
	0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
	};

const Color Atom::elementColors[Atom::ElementsEnd]=
	{
	/* Row 1: Elements H-He: */
	Color(1.0f,1.0f,1.0f), // Hydrogen
	Color::magenta,
	
	/* Row 2: Elements Li-Ne: */
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color(0.0f,0.0f,0.0f), // Carbon
	Color(0.0f,0.0f,1.0f), // Nitrogen
	Color(1.0f,0.0f,0.0f), // Oxygen
	Color::magenta,
	Color::magenta,
	
	/* Row 3: Elements Na-Ar: */
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color(0.5f,0.0f,1.0f), // Phosphorus
	Color(1.0f,1.0f,0.0f), // Sulfur
	Color(0.0f,0.5f,0.0f), // Chlorine
	Color::magenta,
	
	/* Row 4: Elements K-Kr: */
	Color(0.75f,0.25f,1.0f), // Potassium
	Color(0.0f,1.0f,0.0f), // Calcium
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color(1.0f,0.5f,0.0f), // Iron
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	
	/* Row 5: Elements Rb-Xe: */
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color(0.0f,0.25f,0.0f), // Iodine
	Color::magenta,
	
	/* Row 6: Elements Cs-Rn: */
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color(0.75f,0.25f,0.0f), // Gold
	Color::magenta,
	Color::magenta,
	Color(0.5f,0.5f,0.5f), // Lead
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	
	/* Row 7: Elements Fr-Og: */
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta,
	Color::magenta
	};

/*********************
Methods of class Atom:
*********************/

namespace {

/**************
Helper classes:
**************/

class ElementName // Wrapper class to provide operator!= and hash function for element names
	{
	/* Elements: */
	private:
	char name[2]; // The two standard element name characters
	
	/* Constructors and destructors: */
	public:
	ElementName(const char* sName)
		{
		name[0]=sName[0];
		name[1]=sName[1];
		};
	
	/* Methods: */
	friend bool operator!=(const ElementName& en1,const ElementName& en2)
		{
		return en1.name[0]!=en2.name[0]||en1.name[1]!=en2.name[1];
		};
	static size_t hash(const ElementName& value,size_t tableSize)
		{
		size_t hashValue=(size_t(value.name[0])<<8)|size_t(value.name[1]);
		return hashValue%tableSize;
		};
	};

typedef Misc::HashTable<ElementName,Atom::Element,ElementName> ElementNameMap;

/**************
Helper objects:
**************/

bool initialized=false;
ElementNameMap elementNameMap(151);

}

Atom::Element Atom::parseElement(const char* elementName)
	{
	if(!initialized)
		{
		/* Create the element name hash table: */
		for(int i=0;i<118;++i)
			elementNameMap.setEntry(ElementNameMap::Entry(ElementName(elementNames[i]),Element(i)));
		initialized=true;
		}
	
	/* Find the element type matching the given name: */
	ElementNameMap::Iterator ehIt=elementNameMap.end();
	if(elementName[0]!='\0'&&(elementName[1]=='\0'||elementName[2]=='\0'))
		ehIt=elementNameMap.findEntry(ElementName(elementName));
	if(ehIt.isFinished())
		Misc::throwStdErr("Atom::parseElement: Unknown element name \"%s\"",elementName);
	return ehIt->getDest();
	}

void Atom::createTypeMap(Atom::TypeMap& typeMap)
	{
	/* Create a map entry for each element, using the element number as ID: */
	for(TypeID typeId=H;typeId<ElementsEnd;++typeId)
		{
		Type atomType;
		atomType.element=Element(typeId);
		atomType.halfBondLength=getHalfBondLength(Element(typeId));
		atomType.radius=elementVanDerWaalsRadii[typeId];
		atomType.color=elementColors[typeId];
		typeMap.setEntry(TypeMap::Entry(typeId,atomType));
		}
	}

}
