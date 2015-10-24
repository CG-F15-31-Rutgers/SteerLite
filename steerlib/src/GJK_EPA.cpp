/*!
*
* \author VaHiD AzIzI
*
*/


#include "obstacles/GJK_EPA.h"


SteerLib::GJK_EPA::GJK_EPA()
{
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB)
{
	
	bool collisionCheck = GJK(shapeA, shapeB); 
	return collisionCheck;
}

bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB) {
	std::vector<Util::Vector> simplex;
	Util::Vector centerA = shapeCenter(shapeA);
	Util::Vector centerB = shapeCenter(shapeB);
	Util::Vector searchDirection = centerB - centerA;
	simplex.push_back(support(shapeA, shapeB, searchDirection));// get the first Minkowski Difference point
	searchDirection = -searchDirection;// negate direction for the next point
	while (true) {
		simplex.push_back(support(shapeA, shapeB, searchDirection));// add a new point to the simplex because we haven't terminated yet
		if (dotProduct(simplex.back(), searchDirection) <= 0) {
			// if the point added last was not past the origin in the direction of d
			// then the Minkowski Sum cannot possibly contain the origin since
			// the last point added is on the edge of the Minkowski Difference
			return false;
		}
		else if (containsOrigin(simplex, searchDirection)) {//if true, collision detected
			
			return true;
		}
	}
}

Util::Vector SteerLib::GJK_EPA::support(const std::vector<Util::Vector>& shape1, const std::vector<Util::Vector>& shape2, const Util::Vector& direction) {
	Util::Vector point1 = furthestPointInDirection(shape1, direction);
	Util::Vector point2 = furthestPointInDirection(shape2, -direction);
	Util::Vector point3 = point1 - point2;
	return point3;
}

Util::Vector SteerLib::GJK_EPA::shapeCenter(const std::vector<Util::Vector>& shape) {
	Util::Vector center(0, 0, 0);
	for (int i = 0; i < shape.size(); i++) {
		Util::Vector point = shape[i];
		center[0] += point[0];
		center[2] += point[2];
	}
	// since we are in 2d we do not change y
	center[0] = center[0] / (float)shape.size();
	center[2] = center[2] / (float)shape.size();
	return center;
}

bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector>& simplex, Util::Vector& direction) {
	Util::Vector pointA = simplex.back();
	Util::Vector AO = -pointA;
	if (simplex.size() == 3) { //is a triangle
		Util::Vector pointB = simplex[1];
		Util::Vector pointC = simplex[0];
		Util::Vector AB = pointB - pointA;
		Util::Vector AC = pointC - pointA;
		Util::Vector abPerp = AB * dotProduct(AB, AC) - AC * (dotProduct(AB, AB));//computing normals
		Util::Vector acPerp = AC * dotProduct(AC, AB) - AB * (dotProduct(AC, AC));//computing normals
		if (dotProduct(abPerp, AO) > 0) {
			simplex.erase(simplex.begin());//remove point c
			direction = abPerp; //set new direction to abPerp
		}
		else if (dotProduct(acPerp, AO) > 0) {
			simplex.erase(simplex.begin() + 1);//remove point b
			direction = acPerp; //set new direction to acPerp
		}
		else {
			return true;
		}
	}
	else { //simplex is a line segment
		Util::Vector pointB = simplex.at(0);
		Util::Vector AB = pointB - pointA;
		Util::Vector abPerp = AO * dotProduct(AB, AB) - AB * dotProduct(AB, AO);// get the perp to AB in the direction of the origin
		direction = abPerp; // set the direction to abPerp
		if (dotProduct(abPerp, AO) == 0) {
			float ABO = dotProduct(AB, AO);
			if (ABO >= 0 && ABO < dotProduct(AB, AB)) {
				return true;
			}
		}
	}
	return false;
}

Util::Vector SteerLib::GJK_EPA::furthestPointInDirection(const std::vector<Util::Vector>& shape, const Util::Vector& direction) {
	Util::Vector furthestPoint(0, 0, 0);
	float lastPoint = dotProduct(shape[0], direction);
	int lastIndex = 0;
	for (int i = 1; i < shape.size(); i++) {
		float dotProd = dotProduct(shape[i], direction);
		if (dotProd > lastPoint) {
			lastPoint = dotProd;
			lastIndex = i;
		}
	}
	furthestPoint[0] = shape[lastIndex][0];
	furthestPoint[1] = shape[lastIndex][1];
	furthestPoint[2] = shape[lastIndex][2];
	return furthestPoint;
}

float SteerLib::GJK_EPA::dotProduct(const Util::Vector& vector1, const Util::Vector& vector2) {
	float product = 0;
	for (int i = 0; i < 3; i++) {
		product += vector1[i] * vector2[i];
	}
	return product;
}