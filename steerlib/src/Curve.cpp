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

	// Robustness: make sure there is at least two control point: start and end points
	if (checkRobust())
	{


		// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
		float TimeInterval, deltaTime, StartTime, EndTime;
		Point EndPoint;
		Point StartPoint;
		unsigned int CurrentIndex;
		for (unsigned int i = 1; i < controlPoints.size(); i++)
		{
			TimeInterval = controlPoints[i].time - controlPoints[i - 1].time;
			deltaTime = window;
			StartTime = controlPoints[i - 1].time;
			EndTime = controlPoints[i].time;
			CurrentIndex = i;
			StartPoint = controlPoints[i - 1].position;
			for (float CurrTime = StartTime; CurrTime <= EndTime; CurrTime = CurrTime + deltaTime)
			{

				if (type == hermiteCurve)
				{
					EndPoint = useHermiteCurve(CurrentIndex, CurrTime);
				}
				else if (type == catmullCurve)
				{
					EndPoint = useCatmullCurve(CurrentIndex, CurrTime);
				}

				DrawLib::drawLine(StartPoint, EndPoint, curveColor, curveThickness);
				StartPoint = EndPoint;
			}

		}
	}
	return;
#endif
}

bool compare (CurvePoint a, CurvePoint b) {
	return (a.time < b.time);
}

bool equalCheck(CurvePoint a, CurvePoint b)
{
	return a.time == b.time;
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	
	std::sort(controlPoints.begin(), controlPoints.end(), compare);
	controlPoints.erase(unique(controlPoints.begin(), controlPoints.end(), equalCheck), controlPoints.end());
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

// Check Roboustness make sure there is at least two control point: start and end points
bool Curve::checkRobust()
{
	if (type == hermiteCurve) {
		if (controlPoints.size() > 1)
		{
			return true;
		}
	}
	else if (type == catmullCurve){
		if (controlPoints.size() >= 3)
		{
			return true;
		}
	}
	else {
		return false;
	}
	
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	if (time < 0 ) //negative checks 
	{
		return false;
	}

	for (int i = 0; i < controlPoints.size(); i++)
	{
		CurvePoint points = controlPoints[i];
		if(points.time < 0)
		{
			return false;
		}	
		if (points.time > time )
		{
			nextPoint = i;
			return true;
		}
	}

	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	
	Point newPosition, a, b, c, d;
	float normalTime, intervalTime, coEfficient;
	Util::Vector vector;

	// Calculate time interval, and normal time required for later curve calculations
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;

	// Calculate position at t = time on Hermite curve
	
	coEfficient = 1 - (3 * pow(normalTime, 2) - 2 * pow(normalTime, 3));
	a = controlPoints[nextPoint - 1].position * coEfficient;

	coEfficient = (3 * pow(normalTime, 2) - 2 * pow(normalTime, 3));
	b = controlPoints[nextPoint].position * coEfficient;

	coEfficient = (pow(normalTime, 3) - pow(normalTime, 2) - pow(normalTime, 2) + normalTime) * intervalTime;
	vector = controlPoints[nextPoint - 1].tangent * coEfficient;
	c = Point(vector[0], vector[1], vector[2]);

	coEfficient = (pow(normalTime, 3) - pow(normalTime, 2))*intervalTime;
	vector = controlPoints[nextPoint].tangent * coEfficient;
	d = Point(vector[0], vector[1], vector[2]);

	newPosition = a + b + c + d;

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	float a, b, c,d, deltaT, rate;
	Vector s0, s1;
	Point p1, p2, p3;

	// Calculate time interval, and normal time required for later curve calculations
	deltaT = (Curve::controlPoints.at(nextPoint).time - Curve::controlPoints.at(nextPoint - 1).time);  /* time interval */
	rate = (time - Curve::controlPoints.at(nextPoint - 1).time) / deltaT; /* the given time minus controlPoints.at, over the time interval for a usable scale */
	
	// Calculate position at t = time on Catmull-Rom curve
	if (nextPoint == 1) {/* boundary */
		s0 = (controlPoints[nextPoint].position - controlPoints[nextPoint - 1].position) / deltaT;
		s1 = (controlPoints[nextPoint + 1].position - controlPoints[nextPoint - 1].position) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time);
	}
	else if (nextPoint == controlPoints.size() - 1) { /* not boundary */
		s0 = (controlPoints[nextPoint].position - controlPoints[nextPoint - 2].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint - 2].time);
		s1 = (controlPoints[nextPoint].position - controlPoints[nextPoint - 1].position) / deltaT;
	}
	else { /* not boundary */
		s0 = (controlPoints[nextPoint].position - controlPoints[nextPoint - 2].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint - 2].time);
		s1 = (controlPoints[nextPoint + 1].position - controlPoints[nextPoint - 1].position) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time);
	}

	/* formulae from lecture slides/the google, a/b/c/d coefficients in space */
	a = (2 * pow(rate, 3)) - (3 * pow(rate, 2)) + 1;
	b = (-2 * pow(rate, 3)) + (3 * pow(rate, 2));
	c = (pow(rate, 3)) - (2 * pow(rate, 2)) + rate;
	d = (pow(rate, 3)) - (pow(rate, 2));

	// Return result
	return (a * controlPoints[nextPoint - 1].position) + (b*controlPoints[nextPoint].position) + (c*s0*deltaT) + (d*s1*deltaT);

}