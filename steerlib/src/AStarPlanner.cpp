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


#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 0
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


/*Use the below parameters to adjust output for each part of the assignment. 
*/
#define USE_MANHATTAN_DISTANCE false
#define TIE_BREAKER 1
#define WEIGHT 0
#define DIAGNOLS false
#define PRINT_RESULTS true


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

	int AStarPlanner::getIndexFromPoint(Util::Point p){
		return gSpatialDatabase->getCellIndexFromLocation(p);
	}

	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		gSpatialDatabase = _gSpatialDatabase;

		start = getPointFromGridIndex(getIndexFromPoint(start));
		goal = getPointFromGridIndex(getIndexFromPoint(goal));
		agent_path.clear();

		std::map<Util::Point, SteerLib::AStarPlannerNode, comparePoints> node_map;//will hold our mapping of nodes to points
		std::map<SteerLib::AStarPlannerNode, SteerLib::AStarPlannerNode, compareNode> came_from_map;//will keep track of our path

		std::vector<Util::Point> closed_list;
		std::vector<Util::Point> open_list;
		
		int min_f_index = 0;

		SteerLib::AStarPlannerNode start_node(start, 0, heuristics(start, goal), nullptr);
		node_map.emplace(start, start_node);
		open_list.push_back(start_node.point);//add start node to open list

		while (!open_list.empty()){
			//find the node with the lowest f value
			min_f_index = findMinF(node_map, open_list);
			SteerLib::AStarPlannerNode curr_node = node_map.at(open_list[min_f_index]);

			if (curr_node.point == goal){//we found our goal
				reconstructPath(curr_node, agent_path, came_from_map, start, closed_list);
				return true;
			}

			closed_list.push_back(open_list[min_f_index]);// add to closed list
			open_list.erase(open_list.begin() + min_f_index);//remove from open list
			getSuccessors(curr_node.point, goal, node_map, closed_list, open_list, came_from_map);//get list of nodes successors
		}
		return false;//goal not found
	}

	//Function to mind the node in the open list with the smallest F value
	int AStarPlanner::findMinF(std::map<Util::Point, SteerLib::AStarPlannerNode, comparePoints> node_map, std::vector<Util::Point> open_list) {
		int min_f_index = 0;
		int min_f = node_map.at(open_list[0]).f;
		double g_score = node_map.at(open_list[0]).g;

		for (int i = 0; i < open_list.size(); i++){
			if (node_map.at(open_list[i]).f < min_f){
				min_f = node_map.at(open_list[i]).f;
				g_score = node_map.at(open_list[i]).g;
				min_f_index = i;
			}
			else if (node_map.at(open_list[i]).f == min_f){//both nodes have the same f value, break tie based on value of G
				if (TIE_BREAKER == 1){//break ties in favor of larger G value
					if (node_map.at(open_list[i]).f > g_score){
						min_f = node_map.at(open_list[i]).f;
						g_score = node_map.at(open_list[i]).g;
						min_f_index = i;
					}
				}else if (TIE_BREAKER == -1){//break ties in favor of smaller G value
					if (node_map.at(open_list[i]).f < g_score){
						min_f = node_map.at(open_list[i]).f;
						g_score = node_map.at(open_list[i]).g;
						min_f_index = i;
					}
				}
			}
		}
		return min_f_index;
	}

	//builds a path from start to goal once a path has been found
	void AStarPlanner::reconstructPath(SteerLib::AStarPlannerNode curr_node, std::vector<Util::Point>& agent_path, std::map<SteerLib::AStarPlannerNode, SteerLib::AStarPlannerNode, compareNode> came_from_map, Util::Point start, std::vector<Util::Point> closed_list) {
		double total_distance = curr_node.g;
		agent_path.push_back(curr_node.point);

		while (curr_node.point != start){
			curr_node = came_from_map.at(curr_node);
			agent_path.push_back(curr_node.point);

		}
		agent_path.push_back(start);
		std::reverse(agent_path.begin(), agent_path.end());

		if (PRINT_RESULTS){
			std::cout << "\n Length of Solution Path " << total_distance << '\n';
			std::cout << "Number of Expanded Nodes " << closed_list.size() << '\n';
		}
	}

	void AStarPlanner::getSuccessors(Util::Point curr_node, Util::Point goal, std::map<Util::Point, SteerLib::AStarPlannerNode, comparePoints>& node_map, std::vector<Util::Point>& closed_list, std::vector<Util::Point>& open_list, std::map<SteerLib::AStarPlannerNode, SteerLib::AStarPlannerNode, compareNode>& came_from_map) {
		int x;
		int z;
		double initial_g_cost;

		SteerLib::AStarPlannerNode current_node_to_expand = node_map.at(curr_node);
		initial_g_cost = std::numeric_limits<double>::infinity();
		std::vector<Util::Point> successor_points;

		successor_points = getSuccessorPoints(curr_node);

		for (int i = 0; i < successor_points.size(); i++)
		{
			addNodes(successor_points[i], initial_g_cost, current_node_to_expand, goal, node_map, closed_list, open_list, came_from_map);

		}
	}

	//Finds all the points adjacent to the current node
	std::vector<Util::Point> AStarPlanner::getSuccessorPoints(Util::Point current_node_to_expand){
		int node_index;
		unsigned int x, z;
		node_index = getIndexFromPoint(current_node_to_expand);

		std::vector<Util::Point> successorPoints;
		Util::Point successor_point;

		gSpatialDatabase->getGridCoordinatesFromIndex(node_index, x, z);

		int XRangeMin, XRangeMax, ZRangeMin, ZRangeMax;

		XRangeMin = MAX(x - 1, 0);
		XRangeMax = MIN(x + 1, gSpatialDatabase->getNumCellsX());
		ZRangeMin = MAX(z - 1, 0);
		ZRangeMax = MIN(z + 1, gSpatialDatabase->getNumCellsZ());

		for (int i = XRangeMin; i <= XRangeMax; i += GRID_STEP){
			for (int j = ZRangeMin; j <= ZRangeMax; j += GRID_STEP){
				int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);

				if (index != node_index) {
					successor_point = getPointFromGridIndex(index);
				
					if (USE_MANHATTAN_DISTANCE){
						if (successor_point.x == current_node_to_expand.x || successor_point.z == current_node_to_expand.z){
							successorPoints.push_back(successor_point);
						}
					}else{
						successorPoints.push_back(successor_point);
					}
				}
			}
		}

		return successorPoints;
	}

	//Attempt to add each successor to the open list
	void AStarPlanner::addNodes(Util::Point curr_point, double cost, SteerLib::AStarPlannerNode previous_node, Util::Point goal, std::map<Util::Point, SteerLib::AStarPlannerNode, comparePoints>& node_map, std::vector<Util::Point>& closed_list, std::vector<Util::Point>& open_list, std::map<SteerLib::AStarPlannerNode, SteerLib::AStarPlannerNode, compareNode>& came_from_map){
		int node_index; 
		double tentative_g_score;
		node_index = gSpatialDatabase->getCellIndexFromLocation(curr_point);

		if (!canBeTraversed(node_index)){//check to see if adjacent points are obstacles
			return;
		}

		if (node_map.count(curr_point) == 0){//current point is not in the Node Map, so we will add to the Map
			SteerLib::AStarPlannerNode new_node(curr_point, cost, cost, &previous_node);
			node_map.emplace(curr_point, new_node);
		}

		if (std::find(closed_list.begin(), closed_list.end(), curr_point) != closed_list.end()){//if successor is already in the open list, ignore
			return;
		}

		tentative_g_score = previous_node.g + distanceBetween(previous_node.point, curr_point);// length of this path

		if (DIAGNOLS){//add weight to the diagnols for part3 of the assignment
			if (curr_point.x != previous_node.point.x && curr_point.z != previous_node.point.z) //if diagonal
			{
				tentative_g_score = previous_node.g + 2 * distanceBetween(previous_node.point, curr_point);
			}
		}

		if (std::find(open_list.begin(), open_list.end(), curr_point) == open_list.end()){//check to see if successor is in the open list
			open_list.push_back(curr_point);//discovers a new point, add to open list
		}else if (tentative_g_score >= node_map.at(curr_point).g){
			return; //not a better path
		}

		SteerLib::AStarPlannerNode new_node(curr_point, tentative_g_score, tentative_g_score + heuristics(curr_point, goal), &previous_node);
		node_map.erase(curr_point);
		node_map.emplace(curr_point, new_node);//update our map with the new node

		if (came_from_map.count(new_node) != 0){
			came_from_map.erase(new_node);
		}

		came_from_map.emplace(new_node, previous_node);//update our came from map
	}

	//Calculates either the Manhattan distance or Euclidean distance
	double AStarPlanner::heuristics(Util::Point a, Util::Point b)
	{
		if (USE_MANHATTAN_DISTANCE) {
			Util::Vector difference = a - b;
			return WEIGHT*(abs(difference.x) + abs(difference.z));
		}
		else {
			return WEIGHT*(double)distanceBetween(a, b);//default use Euclidean distance
		}
	}
}