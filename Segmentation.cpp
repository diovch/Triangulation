#include "Segmentation.h"
#include "Triangulation.h"
#include <algorithm>
#include <deque>
//#include <set>
#include <iterator>

std::vector<std::set<Triangulation::Triangle>> Segmentation(Triangulation& t, std::map<int, std::set<int>>& Vertexneighbours)
{
	if (t.triangles.size() == 0)
		return {};

	std::sort(t.triangles.begin(), t.triangles.end());

	std::vector<std::set<Triangulation::Triangle>> ConnectedSpaces;
	
	std::map<int, std::set<int>> TriangleEdgeNeighbours;

	// Fill TriangleEdgeNeighbours from t.triangles and Vertexneighbours
	Fill(TriangleEdgeNeighbours, t.triangles, Vertexneighbours);
	
	ConnectedSpaces = DequeAlgo(t.triangles, TriangleEdgeNeighbours);
	
	return ConnectedSpaces;
}

void Fill(std::map<int, std::set<int>>& TriangleEdgeNeighbours,
	std::vector<Triangulation::Triangle>& triangles, std::map<int, std::set<int>>& Vertexneighbours)
{
	for (int TrianInd = 0; TrianInd < triangles.size(); ++TrianInd)
	{
		if (TriangleEdgeNeighbours[TrianInd].size() == 3) // Triangle doesn't have more than 3 edge neighbours
			continue;

		for (int i = 0; i < 3; ++i)
		{
			Triangulation::Triangle& triangle = triangles[TrianInd];

			int e1 = triangle.indices[i];// edge verticies
			int e2 = triangle.indices[(i + 1) % 3];

			int other = triangle.indices[(i + 2) % 3];// other vertex in triangle

			std::set<int>& e1Neighbours = Vertexneighbours[e1];
			std::set<int>& e2Neighbours = Vertexneighbours[e2];

			std::vector<int> intersection;

			// common neighbours for e1 and e2
			std::set_intersection(e1Neighbours.begin(), e1Neighbours.end(),
				e2Neighbours.begin(), e2Neighbours.end(), std::back_inserter(intersection));

			for (const auto& EdgeNeighbour : intersection)
			{
				if (EdgeNeighbour != other) // is a vertex from another triangle found?
				{
					auto it = std::find(triangles.begin(), triangles.end(),
						Triangulation::Triangle(e1, e2, EdgeNeighbour)); // t.triangles has sorted by indicies
					int ind = it - triangles.begin();
					TriangleEdgeNeighbours[TrianInd].insert(ind);
					TriangleEdgeNeighbours[ind].insert(TrianInd);
				}
			}
		}
	}
}

std::vector<std::set<Triangulation::Triangle>> DequeAlgo(std::vector<Triangulation::Triangle>& triangles,
	std::map<int, std::set<int>>& TriangleEdgeNeighbours)
{
	std::vector<std::set<Triangulation::Triangle>> ConnectedSpaces;
	int NumberOfConnctedSpaces = 0;
	//	list of considered triangles
	std::vector<bool> added(triangles.size(), false);
	for (int TrianInd = 0; TrianInd < triangles.size(); ++TrianInd)
	{
		if (added[TrianInd] == true)
			continue;
		added[TrianInd] = true; // triagle with number TrianInd is considered

		std::deque<Triangulation::Triangle> deq;
		Triangulation::Triangle& triangle = triangles[TrianInd];

		deq.push_back(triangle);
		ConnectedSpaces.push_back(std::set<Triangulation::Triangle>{});
		ConnectedSpaces[NumberOfConnctedSpaces].insert(triangle);

		while (!deq.empty())
		{
			Triangulation::Triangle& current = deq.front();
			deq.pop_front();

			auto it = std::find(triangles.begin(), triangles.end(), current); // t.triangles has sorted by indicies
			int ind = it - triangles.begin();
			std::vector<int> neighbours(TriangleEdgeNeighbours[ind].begin(), TriangleEdgeNeighbours[ind].end());

			for (int i = 0; i < neighbours.size(); ++i)
			{
				if (added[neighbours[i]] == true)
					continue;
				added[neighbours[i]] = true;
				deq.push_back(triangles[neighbours[i]]);
				ConnectedSpaces[NumberOfConnctedSpaces].insert(triangles[neighbours[i]]);
			}
		}

		NumberOfConnctedSpaces++;
	}
	return ConnectedSpaces;
}