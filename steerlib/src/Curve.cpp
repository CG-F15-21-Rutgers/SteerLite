//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;
int a; Point newPosition,oldPosition,temp;
float t;
Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
	
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	Point c, b;
	window = 1;
	//std::cout << "\tWork";
	Util::DrawLib::glColor(curveColor);
	sortControlPoints();
	c = (*controlPoints.begin()).position;
	std::vector<CurvePoint>::iterator it = controlPoints.begin();
	std::vector<CurvePoint>::iterator it2 = controlPoints.begin();
	it2++;
	float j = (*it).time;
	float l = (*it2).time;
	//std::cout << "\tMax Size : "<<controlPoints.size();
	for (int k = 1; k < controlPoints.size(); k++)
	{
		//std::cout << "\tWork2";
		//std::cout << "\nj :" << j << "  l:" << l;

		for (; j < l; j = j + window)  
		{
			Util::DrawLib::drawLine(c, useHermiteCurve(k, j));
			c = useHermiteCurve(k, j);

			//std::cout << "\tWork3";

			//std::cout << "j :"<<j<<"  l:"<<l;
		}
	it2++;
	it++;
	j = (*it).time;
	l = (*it2).time;
	}
	
	
	//std::cout << "\n O: " << oldPosition << "\tN: " << newPosition;
	
	//if ( (int)t % window==0)
	//Util::DrawLib::drawLine(a, newPosition);
	

	
	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	
	return;
#endif
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::sort(controlPoints.begin(), controlPoints.end());
	
	//testing sort function
	//std::cout << "Sorted List:";
	//for (std::vector<CurvePoint>::iterator it = controlPoints.begin(); it != controlPoints.end();it++)
		//std::cout << "x: "<<(*it).position.x<<"\n";
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;
	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{   
	int p=0;
	sortControlPoints();
	//int p = (*controlPoints.end()).time;
	for (std::vector<CurvePoint>::iterator it = controlPoints.begin(); it != controlPoints.end(); it++)
		p++;
	if (p < 2 )
	{
		return false;
	}
	
	
		return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	t = time;
	a = 0;
	for (std::vector<CurvePoint>::iterator it = controlPoints.begin();; it++)
	{
		
		if ((*it).time > time)
			break;
		a++;
		if (it == controlPoints.end())
		{
			a++;
			break;
		}
	}
	nextPoint = a;

	std::cout << "\nCurrent time:"<<time<<" Next Point:" << a;
	return true;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{   
	oldPosition = newPosition;

	float normalTime, intervalTime,u;
	//Hermite Curve Implementation
	std::vector<CurvePoint>::iterator it2 = controlPoints.begin();
	std::vector<CurvePoint>::iterator it = controlPoints.begin();
	int z = 1;
	for (; z < nextPoint || z == nextPoint; z++)
	{
		if (z < nextPoint)
			it++;
		it2++;
	}
	//std::cout << "\nTime :" << time;
	normalTime = float((*it2).time - (*it).time);
	intervalTime = time - (*it).time;
	u = intervalTime / normalTime;
	float ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
	float sx = (*it).tangent.x;
	float sx2 = (*it2).tangent.x;

	float sy = (*it).tangent.y;
	float sy2 = (*it2).tangent.y;

	float sz = (*it).tangent.z;
	float sz2 = (*it2).tangent.z;

	ax = (-1 * 2 * ((*it2).position.x - (*it).position.x) / pow(((*it2).time - (*it).time), 3)) + ((sx + sx2) / pow(((*it2).time - (*it).time), 2));
	bx = (3 * ((*it2).position.x - (*it).position.x) / pow(((*it2).time - (*it).time), 2)) - ((2 * sx + sx2) / pow(((*it2).time - (*it).time), 1));
	cx = sx;
	dx = (*it).position.x;
    ay = (-1 * 2 * ((*it2).position.y - (*it).position.y) / pow(((*it2).time - (*it).time), 3)) + ((sy + sy2) / pow(((*it2).time - (*it).time), 2));
	by = (3 * ((*it2).position.y - (*it).position.y) / pow(((*it2).time - (*it).time), 2)) - ((2 * sy + sy2) / pow(((*it2).time - (*it).time), 1));
	cy = sy;
	dy = (*it).position.y;
    az = (-1 * 2 * ((*it2).position.z - (*it).position.z) / pow(((*it2).time - (*it).time), 3)) + ((sz + sz2) / pow(((*it2).time - (*it).time), 2));
	bz = (3 * ((*it2).position.z - (*it).position.z) / pow(((*it2).time - (*it).time), 2)) - ((2 * sz + sz2) / pow(((*it2).time - (*it).time), 1));
	cz = sz;
	dz = (*it).position.z;

	newPosition.x = ax*pow((time - (*it).time), 3) + bx*pow((time - (*it).time), 2) + cx*(time - (*it).time) + dx;
	newPosition.y = ay*pow((time - (*it).time), 3) + by*pow((time - (*it).time), 2) + cy*(time - (*it).time) + dy;
	newPosition.z = az*pow((time - (*it).time), 3) + bz*pow((time - (*it).time), 2) + cz*(time - (*it).time) + dz;

	
//	if (newPosition.x == (*controlPoints.end()).position.x && newPosition.y == (*controlPoints.end()).position.y && newPosition.z == (*controlPoints.end()).position.z)
	//	std::cout << "Hello";
	
	if (a == controlPoints.size())
		return (*controlPoints.begin()).position;
	
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	std::vector<CurvePoint>::iterator it0 = controlPoints.begin();
	std::vector<CurvePoint>::iterator it = controlPoints.begin();
	std::vector<CurvePoint>::iterator it2 = controlPoints.begin();
	std::vector<CurvePoint>::iterator it3 = controlPoints.begin();
    int z = 1;
	float u,normalTime,intervalTime,h1,h2,h3,h4;
	for (; z < nextPoint; z++)
	{
		it++;
	}
	it0 = it; it0--;
	it2 = it; it2++;
	it3 = it2; it3++;
	//std::cout << "\ntime:" << time << "\tnextpoint :" << nextPoint << "\tz:" << z << " \tit :" << (*it).time << " \tit2 :" << (*it2).time << " \tit3 :" << (*it3).time << " \tit0 :" << (*it0).time;

	normalTime = (*it2).time - (*it).time;
	if (it2 == controlPoints.end())
	{
		
		(*it3).position.x*=0;
		(*it3).position.y*=0;
		(*it3).position.z*=0;
		
	}
	
	if (it == controlPoints.begin())
	{     
		(*it).tangent.x = (0.5*((*it2).position.x - 0));
		(*it).tangent.y = (0.5*((*it2).position.y - 0));
		(*it).tangent.z = (0.5*((*it2).position.z - 0));
	}

	else if (it == controlPoints.end())
	{
		(*it).tangent.x = (0.5*(0 - (*it0).position.x));
		(*it).tangent.y = (0.5*(0 - (*it0).position.y));
		(*it).tangent.z = (0.5*(0 - (*it0).position.z));
	}
	else
	{
		(*it).tangent.x = (0.5*((*it2).position.x - (*it0).position.x));
		(*it).tangent.y = (0.5*((*it2).position.y- (*it0).position.y));
		(*it).tangent.z = (0.5*((*it2).position.z - (*it0).position.z));
	}
	intervalTime = time - (*it).time;
	u = intervalTime / normalTime;

	h1 = (-0.5 * pow(u, 3)) + (pow(u, 2)) - (0.5 * u);
	h2 = (1.5* pow(u, 3)) + (-2.5 * pow(u, 2)) + 1;
	h3 = (-1.5*pow(u, 3)) + (2 * pow(u, 2)) + (0.5 * u);
	h4 = (pow(u, 3)* 0.5) - (0.5*pow(u, 2));

	newPosition.x = ((*it0).position.x * h1) + ((*it).position.x * h2) + ((*it2).position.x * h3) + ((*it3).position.x * h4);
	newPosition.y = ((*it0).position.y * h1) + ((*it).position.y * h2) + ((*it2).position.y * h3) + ((*it3).position.y * h4);
	newPosition.z = ((*it0).position.z * h1) + ((*it).position.z * h2) + ((*it2).position.z * h3) + ((*it3).position.z * h4);

	return newPosition;
}