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
	std::vector<Util::Vector> simplex;
	bool collisionCheck = GJK(shapeA, shapeB, simplex);
	if (collisionCheck == true) {//found collision
		EPA(return_penetration_depth, return_penetration_vector, shapeA, shapeB, simplex);
		return true;
	}
	else {// no collision
		return false;
	}
}

//Begin GJK algorithm
bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB, std::vector<Util::Vector>& simplex) {
	
	Util::Vector centerA = shapeCenter(shapeA);
	Util::Vector centerB = shapeCenter(shapeB);
	Util::Vector searchDirection = centerB - centerA;
	simplex.push_back(support(shapeA, shapeB, searchDirection));// get the first Minkowski Difference point
	searchDirection = -searchDirection;// negate direction for the next point
	while (true) {
		simplex.push_back(support(shapeA, shapeB, searchDirection));// add a new point to the simplex because we haven't terminated yet
		if (dot(simplex.back(), searchDirection) <= 0) {
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
Util::Vector crossProduct(Util::Vector v1, Util::Vector v2) {
	Util::Vector prod;
	prod.x = v1.y*v2.z - v1.z*v2.y;
	prod.y = v1.z*v2.x - v1.x*v2.z;
	prod.z = v1.x*v2.y - v1.y*v2.x;
	return prod;
}

Util::Vector tripleProduct(Util::Vector v1, Util::Vector v2, Util::Vector v3) {
	return crossProduct(crossProduct(v1, v2), v3);
}
bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector>& simplex, Util::Vector& direction) {
	
	Util::Vector pointA = simplex.back();
	Util::Vector AO = -pointA;
	if (simplex.size() == 3) { //is a triangle
		Util::Vector pointB = simplex.at(0);
		Util::Vector pointC = simplex.at(1);

		Util::Vector AB = pointB - pointA;
		Util::Vector AC = pointC - pointA;
		Util::Vector abPerp = tripleProduct(AC, AB, AB);
		if (dot(abPerp, pointC) >= 0) { // new direction might be going "wrong way", check against other side
			abPerp = -abPerp;
		}
		
		if (dot(abPerp, AO) > 0) {
			simplex.erase(std::find(simplex.begin(), simplex.end(), pointC));//remove point c
			direction = abPerp; //set new direction to abPerp
		}
		else {
			Util::Vector AC = pointC - pointA;
			Util::Vector acPerp = tripleProduct(AB, AC, AC);
			if (dot(acPerp, pointB) >= 0) {
				acPerp = -acPerp;
			}
			if (dot(acPerp, AO) <= 0) {
				return true;
			}
			simplex.erase(std::find(simplex.begin(), simplex.end(), pointB));//remove point b
			direction = acPerp; //set new direction to acPerp
		}
		
	}
	else if (simplex.size() == 2) { //simplex is a line segment
		Util::Vector pointB = simplex.at(0);
		Util::Vector AB = pointB - pointA;
		Util::Vector abPerp = tripleProduct(AB, AO, AB);
		if (dot(abPerp, AO) < 0) {
			abPerp = -abPerp;
		}
		// set direction to the perpendicular vector, this is the direction where we want to look for a third point
		direction = abPerp;
	}
	return false;
}
//End GJK algorithm

Util::Vector SteerLib::GJK_EPA::shapeCenter(const std::vector<Util::Vector>& shape) {
	Util::Vector center(0, 0, 0);
	for (int i = 0; i < shape.size(); i++) {
		Util::Vector point = shape[i];
		center[0] += point[0];
		center[1] += point[1];
		center[2] += point[2];
	}
	
	center[0] = center[0] / (float)shape.size();
	center[1] = center[1] / (float)shape.size();
	center[2] = center[2] / (float)shape.size();
	return center;
}

Util::Vector SteerLib::GJK_EPA::support(const std::vector<Util::Vector>& shape1, const std::vector<Util::Vector>& shape2, const Util::Vector& direction) {
	Util::Vector point1 = furthestPointInDirection(shape1, direction);
	Util::Vector point2 = furthestPointInDirection(shape2, -direction);
	Util::Vector point3 = point1 - point2;
	return point3;
}

Util::Vector SteerLib::GJK_EPA::furthestPointInDirection(const std::vector<Util::Vector>& shape, const Util::Vector& direction) {
	Util::Vector furthestPoint(0, 0, 0);
	float lastPoint = dot(shape[0], direction);
	int lastIndex = 0;
	for (int i = 1; i < shape.size(); i++) {
		float dotProd = dot(shape[i], direction);
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

//start EPA algorithm
 struct Edge {
	Util::Vector normal;
	float distance;
	size_t index;
};

 Edge FindClosestEdge(const std::vector<Util::Vector>& simplex) {
	 Edge closest;

	 // set distance to maximum of a float
	 closest.distance = FLT_MAX;

	 for (size_t i = 0; i < simplex.size(); ++i) {
		 size_t j = i + 1;
		 if (j == simplex.size()) {
			 j = 0;
		 }
		 // i is the current point in the simplex and j is the next

		 const Util::Vector& a = simplex.at(i);
		 const Util::Vector& b = simplex.at(j);

		 // find our edge
		 Util::Vector e = b - a;

		 // find the normal of the edge
		 Util::Vector n = tripleProduct(e, a, e);
		 if (dot(n, -a) >= 0) { // new direction might be going "wrong way", check against other side
			 n = -n;
		 }

		 n *= 1 / sqrt(dot(n, n)); // normalize

		 float d = dot(a, n); // distance between orgin and normal (a is the vector between origin and a since: a - origin = a - [0,0] = a)

							  // if d is less than the previous closest
		 if (d < closest.distance) {
			 closest.distance = d;
			 closest.normal = n;
			 closest.index = j; 
		 }
	 }

	 return closest;
 }

void SteerLib::GJK_EPA::EPA(float& penetration_depth, Util::Vector& penetration_vector, const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB, std::vector<Util::Vector>& simplex) {
	
	while (true) {
		Edge e = FindClosestEdge(simplex);

		// find the furthest minkowski difference point in the direction of the normal
		const Util::Vector& p = support(shapeA, shapeB, e.normal);

		// find the distance between the point and the edge
		float d = dot(p, e.normal);
		
		// find if we've hit the border of the minkowski difference
		if (fabs(d - e.distance) < 0.01f) {
			penetration_depth = d;
			penetration_vector = e.normal;
			return;
		}
		else {
			// add the point inbetween the points where it was found (due to the need of correct winding)
			simplex.insert(simplex.begin() + e.index, p);
		}
	}

	return;
}



