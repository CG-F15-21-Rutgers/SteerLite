//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"
#include <limits.h>
#include <cmath>

#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MANHATTAN 1                // 1 for Manhattan distance, 0 for Eucledian
namespace SteerLib
{
	AStarPlanner::AStarPlanner() {}

	AStarPlanner::~AStarPlanner() {}

	bool AStarPlanner::canBeTraversed(int id)
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x, z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = MAX(x - OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z - OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());


		for (int i = x_range_min; i <= x_range_max; i += GRID_STEP)
		{
			for (int j = z_range_min; j <= z_range_max; j += GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);
				traversal_cost += gSpatialDatabase->getTraversalCost(index);

			}
		}

		if (traversal_cost > COLLISION_COST)
			return false;
		return true;
	}



	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}
	double heuristic_cost_estimate(AStarPlannerNode Astart, AStarPlannerNode Agoal)
	{
		double cost;
		if (MANHATTAN)
		{
			cost = abs(Astart.point.x - Agoal.point.x) + abs(Astart.point.z - Agoal.point.z);
		}
		else
		{
			cost = sqrt(pow((Astart.point.x - Agoal.point.x), 2) + pow((Astart.point.z - Agoal.point.z), 2));
		}
		return cost;
	}

	AStarPlannerNode lowestFScore(std::vector<AStarPlannerNode> &open, int &index)
	{
		std::vector<AStarPlannerNode>::iterator itr = open.begin();
		index = INT_MAX;
		AStarPlannerNode temp;
		for (int i = 0; itr != open.end(); itr++, i++)
		{
			if (temp > *itr) {
				index = i;
				temp = *itr;
			}
			if (temp.f == (*itr).f)
			{
				if (temp.g < (*itr).g)
				{
					temp = *itr;
					index = i;
				}
			}

		}

		return temp;
	}
	std::vector<AStarPlannerNode> AStarPlanner::getNeighbours(AStarPlannerNode current, AStarPlannerNode goal)
	{    
		int currentIndex = gSpatialDatabase->getCellIndexFromLocation(current.point);
		std::vector<AStarPlannerNode> neighbours;
		unsigned int currentPointx, currentPointz;
		gSpatialDatabase->getGridCoordinatesFromIndex(currentIndex, currentPointx, currentPointz);
		Util::Point currentPoint(currentPointx, 0, currentPointz);
		Util::Point temp(0, 0, 0);
		for (int i = 0; i < 8; i++)
		{
			temp = currentPoint;
			if (i == 0)
				temp.x++;
			if (i == 1)
				temp.x--;
			if (i == 2)
				temp.z++;
			if (i == 3)
				temp.z--;
			if (i == 4)
			{
				temp.x++; temp.z++;
			}
			if (i == 5)
			{
				temp.x++; temp.z--;
			}
			if (i == 6)
			{
				temp.x--; temp.z--;
			}
			if (i == 7)
			{
				temp.x--; temp.z++;
			}
			int index = gSpatialDatabase->getCellIndexFromGridCoords(temp.x, temp.z);
			temp = getPointFromGridIndex(index);
			if (canBeTraversed(index))
			{
				//std::cout << "\t\t Going in";
				
				AStarPlannerNode x;
				x.point = temp;
				//			if(i<4)
				x.g = current.g + 1;
				//		else
				//	x.g = current.g + sqrt(2);
				x.f = x.g + (8 * heuristic_cost_estimate(x, goal));
				x.parentx = current.point.x;
				x.parentz = current.point.z;
				neighbours.push_back(x);
			}

		}

		return neighbours;
	}
	bool inClosedList(AStarPlannerNode check, std::vector<AStarPlannerNode> closed)
	{
		for (std::vector<AStarPlannerNode>::iterator itr = closed.begin(); itr != closed.end(); itr++)
		{
			if (check == (*itr))
				return true;
		}
		return false;
	}
	int inOpenList(AStarPlannerNode check, std::vector<AStarPlannerNode> closed)
	{
		int i = 0;
		for (std::vector<AStarPlannerNode>::iterator itr = closed.begin(); itr != closed.end(); itr++, i++)
		{
			if (check == (*itr))
				return i;
		}
		i = -1;
		return i;
	}
	AStarPlannerNode findclosed(std::vector<AStarPlannerNode> closed, AStarPlannerNode find)
	{
		for (std::vector<AStarPlannerNode>::iterator itr = closed.begin(); itr != closed.end(); itr++)
		{
			if ((*itr).point.x == find.parentx && (*itr).point.z == find.parentz)
				return (*itr);
		}
	}


	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		gSpatialDatabase = _gSpatialDatabase;

		int startIndex = _gSpatialDatabase->getCellIndexFromLocation(start);
		start = getPointFromGridIndex(startIndex);
		AStarPlannerNode Astart(start, 0, 0, NULL);
		int goalIndex = _gSpatialDatabase->getCellIndexFromLocation(goal);
		goal = getPointFromGridIndex(goalIndex);
		AStarPlannerNode Agoal(goal, DBL_MAX, DBL_MAX, NULL);
		std::vector<AStarPlannerNode> open;
		std::vector<AStarPlannerNode> closed;
		Astart.g = 0;
		Astart.f = Astart.g + heuristic_cost_estimate(Astart, Agoal);
		open.emplace_back(Astart);
		int index;
		std::vector<AStarPlannerNode> neighbours;
		AStarPlannerNode current;
		while (!open.empty())
		{
			neighbours.clear();
			current = lowestFScore(open, index);
			if (current.point == Agoal.point)
			{
				//std::cout << "\nPath length:" << current.g;

				while (!(current.parentx == 0 && current.parentz == 0))
				{
					agent_path.push_back(current.point);
					current = findclosed(closed, current);
				}
				std::reverse(agent_path.begin(), agent_path.end());
				//std::cout << "\n\tNumber of nodes in Closed List: " << closed.size();
				return TRUE;
			}
			open.erase(open.begin() + index);

			closed.push_back(current);

			neighbours = getNeighbours(current, Agoal);
			double tentative_gscore = 0;
			int i; int z = 0;
			for (std::vector<AStarPlannerNode>::iterator itr = neighbours.begin(); itr != neighbours.end(); itr++, z++)
			{
				if (inClosedList(*itr, closed))
					continue;
				tentative_gscore = current.g + 1;
				i = inOpenList(*itr, open);
				if (i != -1)
				{
					if (tentative_gscore < open.at(i).g)
					{
						open.at(i).g = tentative_gscore;
						open.at(i).parentx = current.point.x;
						open.at(i).parentz = current.point.z;
						open.at(i).f = open.at(i).g + heuristic_cost_estimate(open.at(i), Agoal);
					}
				}
				else
				{
					open.push_back(neighbours.at(z));
				}
			}

		}
		std::cout << "\n No path";
		return false;
	}
}