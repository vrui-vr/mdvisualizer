/***********************************************************************
Color - Class defining a color as an RGB triple.
Copyright (c) 2019-2025 Oliver Kreylos

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

#include "Color.h"

namespace MD {

/******************************
Static elements of class Color:
******************************/

const Color Color::black(0.0f,0.0f,0.0f);
const Color Color::darkGrey(0.25f,0.25f,0.25f);
const Color Color::grey(0.5f,0.5f,0.5f);
const Color Color::lightGrey(0.75f,0.75f,0.75f);
const Color Color::white(1.0f,1.0f,1.0f);
const Color Color::red(1.0f,0.0f,0.0f);
const Color Color::yellow(1.0f,1.0f,0.0f);
const Color Color::green(0.0f,1.0f,0.0f);
const Color Color::cyan(0.0f,1.0f,1.0f);
const Color Color::blue(0.0f,0.0f,1.0f);
const Color Color::magenta(1.0f,0.0f,1.0f);

}
