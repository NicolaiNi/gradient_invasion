//============================================================================
// Emanuel A. Fronhofer, Nicolai Nitsche & Florian Altermatt
// Information use shapes range expansion dynamics into environmental gradients
// Global Ecology and Biogepgraphy
// 2016
//
// classes and constants
//============================================================================

/*
	Copyright (C) 2016  Emanuel A. Fronhofer

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//____________________________________________________________________________
//----------------------------------------------------------- define constants

const int WORLDDIM_X = 100;														// world dimensions, x
const int RS = 2;																// random seed
const int INIT_ROWS = 5;														// how many rows to initialize with individuals
const int CM_AREA = 5;															// definition of x-distance of core and marginal area for emigration rate analysis

//____________________________________________________________________________
//------------------------------------------------------------- define classes

// one individual ------------------------------------------------------------
class TIndiv {
public:
	TIndiv();
	//float dispRate;
	float alpha;
	float lambda;
};

TIndiv::TIndiv() { //constructor for TIndiv
	//dispRate = 0;
	alpha = 0;
	lambda = 0;
}

// one patch -----------------------------------------------------------------
class TPatch {
public:
	TPatch();
	vector<TIndiv> females;
	vector<TIndiv> newFemales;
};

TPatch::TPatch() {
	females.clear();
	newFemales.clear();
}
