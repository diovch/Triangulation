#include "Segmentation.h"
#include "Triangulation.h"
#include <algorithm>
#include <deque>
#include <iterator>

std::vector<std::set<Triangulation::Triangle>> Segmentation(Triangulation& t, std::map<int, std::set<int>>& Vertexneighbours)
{
	if (t.triangles.size() == 0)
		return {};

	std::sort(t.triangles.begin(), t.triangles.end());

	std::vector<std::set<Triangulation::Triangle>> ConnectedSpaces;
	std::map<int, std::set<int>> TriangleEdgeNeighbours;

	
	Fill(TriangleEdgeNeighbours, t.triangles, Vertexneighbours);
	Research(TriangleEdgeNeighbours, t);

	ConnectedSpaces = DequeAlgo(t.triangles, TriangleEdgeNeighbours);
	
	return ConnectedSpaces;
}

void Fill(std::map<int, std::set<int>>& TriangleEdgeNeighbours,
	std::vector<Triangulation::Triangle>& triangles, std::map<int, std::set<int>>& Vertexneighbours)
{
	for (int TrianInd = 0; TrianInd < triangles.size(); ++TrianInd)
	{
		//if (TriangleEdgeNeighbours[TrianInd].size() == 3) // Triangle doesn't have more than 3 edge neighbours
		//	continue;

		for (int VertexInd = 0; VertexInd < 3; ++VertexInd)
		{
			Triangulation::Triangle& triangle = triangles[TrianInd];

			int Vertex1 = triangle.indices[VertexInd];
			int Vertex2 = triangle.indices[(VertexInd + 1) % 3];

			int OtherVertex = triangle.indices[(VertexInd + 2) % 3];

			std::set<int>& Vertex1Neighbours = Vertexneighbours[Vertex1];
			std::set<int>& Vertex2Neighbours = Vertexneighbours[Vertex2];

			std::vector<int> CommonNeighbours;

			
			std::set_intersection(Vertex1Neighbours.begin(), Vertex1Neighbours.end(),
				Vertex2Neighbours.begin(), Vertex2Neighbours.end(), std::back_inserter(CommonNeighbours));

			for (const auto& EdgeNeighbour : CommonNeighbours)
			{
				if (EdgeNeighbour != OtherVertex)
				{
					auto it = std::find(triangles.begin(), triangles.end(),
						Triangulation::Triangle(Vertex1, Vertex2, EdgeNeighbour)); // t.triangles has sorted by indicies
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
	int NumberOfConnectedSpaces = 0;
	//	list of considered triangles
	std::vector<bool> IsConsidered(triangles.size(), false);
	for (int TrianInd = 0; TrianInd < triangles.size(); ++TrianInd)
	{
		if (IsConsidered[TrianInd] == true)
			continue;
		IsConsidered[TrianInd] = true; 

		std::deque<Triangulation::Triangle> TriangleQueue;
		Triangulation::Triangle& triangle = triangles[TrianInd];

		TriangleQueue.push_back(triangle);
		ConnectedSpaces.push_back(std::set<Triangulation::Triangle>{});// it requires initialization for further set insertion
		ConnectedSpaces[NumberOfConnectedSpaces].insert(triangle);

		while (!TriangleQueue.empty())
		{
			Triangulation::Triangle& CurrentTriangle = TriangleQueue.front();
			TriangleQueue.pop_front();

			auto it = std::find(triangles.begin(), triangles.end(), CurrentTriangle); // t.triangles has sorted by indicies
			int ind = it - triangles.begin();
			std::vector<int> TriangleNeighbours(TriangleEdgeNeighbours[ind].begin(), TriangleEdgeNeighbours[ind].end()); //coping

			for (int i = 0; i < TriangleNeighbours.size(); ++i)
			{
				if (IsConsidered[TriangleNeighbours[i]] == true)
					continue;
				IsConsidered[TriangleNeighbours[i]] = true;
				TriangleQueue.push_back(triangles[TriangleNeighbours[i]]);
				ConnectedSpaces[NumberOfConnectedSpaces].insert(triangles[TriangleNeighbours[i]]);
			}
		}
		NumberOfConnectedSpaces++;
	}
	return ConnectedSpaces;
}

void Research(std::map<int, std::set<int>>& neighbours, Triangulation& t)
{
	std::map<int, std::set<int>> slice;
	for (auto& x : neighbours)
	{
		if (x.second.size() > 3)
		{
			slice[x.first] = x.second;
		}
	}

	//for(auto& x)
	return;
}