#include "obstacles/GJK_EPA.h"
#include<math.h>
#define TOLERANCE 0.1
struct Edge
{
	float distance;
	Util::Vector normal;
	int index;
};
SteerLib::GJK_EPA::GJK_EPA()
{
}

Util::Vector farthestpoint(const std::vector<Util::Vector>& shape, Util::Vector d)
{
	float temp2 = d*shape[0];
	int k = 0, temp1;
	for (int i = 1; i < shape.size(); i++)
	{
		temp1 = d*shape[i];
		if (temp1 > temp2)
		{
			temp2 = temp1;
			k = i;
		}
	}
	return (shape[k]);
}

Util::Vector support(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, Util::Vector d)
{
	Util::Vector first = farthestpoint(_shapeA, d);
	Util::Vector nd = -1 * d;
	Util::Vector second = farthestpoint(_shapeB, nd);
	Util::Vector support = first - second;

	return support;
}
bool containsOrigin(Util::Vector &d, std::vector<Util::Vector>& simplex)
{
	Util::Vector a, ao, b, c, ab, ac;
	a = simplex.back();
	ao = -1 * a;
	if (simplex.size() == 3)
	{
		b = simplex[1];
		c = simplex[0];
		ab = b - a;
		ac = c - a;
		d = Util::Vector(ab.z, ab.y, -1 * ab.x);
		if (d*c>0)
		{
			d = -1 * d;
		}
		if (d*ao > 0)
		{
			simplex.erase(simplex.begin());
			return false;
		}
		d = Util::Vector(ac.z, ac.y, -1 * ac.x);
		if (d*b > 0)
		{

			d = -1 * d;
		}
		if (d*ao > 0)
		{
			simplex.erase(simplex.begin());
			return false;
		}
		return true;
	}
	else
	{
		b = simplex[0];
		ab = b - a;
		d = Util::Vector(ab.z, ab.y, -1 * ab.x);

		if (d*ao < 0)
		{
			d = -1 * d;

		}

	}
	return false;
}

bool GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{



	Util::Vector d(1, 0, -1);

	simplex.push_back(support(_shapeA, _shapeB, d));
	Util::Vector nd = -1 * d;
	while (1)
	{
		simplex.push_back(support(_shapeA, _shapeB, nd));
		if (simplex.back() * nd <= 0)
		{
			return false;
		}
		else
		{
			if (containsOrigin(nd, simplex))
			{
				return true;
			}

		}
	}

}
void normalize(Util::Vector &n)
{
	float length = sqrt(((n.x*n.x) + (n.z*n.z)));

	if (length != 0) {
		n.x = n.x / length;
		n.z = n.z / length;
	}

}


Edge findClosestEdge(std::vector<Util::Vector>& simplex)
{
	Util::Vector a, b, e, oa, n, v1, v2;
	float d, d1, d2;
	Edge closest;
	closest.distance = 999999;
	for (int i = 0, j; i < simplex.size(); i++)
	{
		if (i + 1 == simplex.size())
			j = 0;
		else
			j = i + 1;
		a = simplex[i];
		b = simplex[j];
		e = b - a;
		oa = a;

		Util::Vector n = oa*(e*e) - e*(e*oa); //n = vectorTripleProduct(e, oa, e);
		normalize(n);
		d = abs(dot(n, a));
		if (d < closest.distance)
		{
			closest.distance = d;
			closest.normal.x = n.x;
			closest.normal.z = n.z;
			closest.normal.y = 0;
			closest.index = j;
		}
	}
	return closest;
}

float EPA(Util::Vector &penetrationvector, const std::vector<Util::Vector>& shapeA, const std::vector<Util::Vector>& shapeB, std::vector<Util::Vector>& simplex)
{
	while (1)
	{

		Edge e = findClosestEdge(simplex);
		Util::Vector p = support(shapeA, shapeB, e.normal);
		float d = dot(p, e.normal);
		if (d - e.distance <= TOLERANCE)
		{
			penetrationvector = e.normal;
			return d;
		}
		else
		{
			simplex.insert(simplex.begin() + e.index, p);
		}
	}
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> simplex;
	bool is_colliding;
	is_colliding = GJK(_shapeA, _shapeB, simplex);
	if (is_colliding)
	{
		return_penetration_depth = EPA(return_penetration_vector, _shapeA, _shapeB, simplex);
		/*	if (return_penetration_depth == 0)
		return false;
		*/	return true;
	}
	else
	{
		return false;
	}
	// There is no collision 
}