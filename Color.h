/***********************************************************************
Color - Class defining a color as an RGB triple.
Copyright (c) 2019 Oliver Kreylos
***********************************************************************/

#ifndef MD_COLOR_INCLUDED
#define MD_COLOR_INCLUDED

namespace MD {

class Color
	{
	/* Elements: */
	private:
	float rgb[3]; // RGB color components in the range [0.0, 1.0]
	
	/* Some oft-used colors: */
	public:
	static const Color black;
	static const Color darkGrey;
	static const Color grey;
	static const Color lightGrey;
	static const Color white;
	static const Color red;
	static const Color yellow;
	static const Color green;
	static const Color cyan;
	static const Color blue;
	static const Color magenta;
	
	/* Constructors and destructors: */
	Color(void) // Dummy constructor
		{
		}
	Color(float sR,float sG,float sB) // Creates a color from RGB components
		{
		rgb[0]=sR;
		rgb[1]=sG;
		rgb[2]=sB;
		}
	Color(const float sRgb[3]) // Ditto, from array
		{
		for(int i=0;i<3;++i)
			rgb[i]=sRgb[i];
		}
	
	/* Methods: */
	const float* getRgb(void) const // Returns the color array
		{
		return rgb;
		}
	float operator[](int componentIndex) const // Returns a color component
		{
		return rgb[componentIndex];
		}
	float& operator[](int componentIndex) // Ditto
		{
		return rgb[componentIndex];
		}
	};

}

#endif
