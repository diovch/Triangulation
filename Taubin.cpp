#include "Taubin.h"
#include "R3Graph.h"

#include <algorithm>
#include <iterator>

void Taubin(Triangulation& t, std::map<int, std::set<int>>& neighbours, double lambda, double mu, int iterations)
{
	for (int it = 0; it < iterations; ++it)
	{
		Update(t, neighbours, lambda);
		Update(t, neighbours, mu);
	}

	//NormalsUpdate(t);
	
	return;
}

void Update(Triangulation& t, std::map<int, std::set<int>>& neighbours, double factor)
{
	std::vector<R3Graph::R3Vector> meshBias(t.vertices.size(), { 0.,0.,0. });

	for (int i = 0; i < t.vertices.size(); ++i)
	{
		double weight = 1. / neighbours[i].size();
		R3Graph::R3Vector average_vector = { 0., 0., 0. };

		for (const auto ind : neighbours[i])
			average_vector += (t.vertices[ind].point - t.vertices[i].point) * weight;
		meshBias[i] = average_vector;
	}

	for (int i = 0; i < t.vertices.size(); ++i)
		t.vertices[i].point += meshBias[i] * factor;
}

void TaubinSkeleton(Triangulation& t, std::map<int, std::set<int>>& neighbours, double lambda, double mu, int iterations)
{
	std::map<int, std::set<Triangulation::Edge>> skeletonEdges;

	for (int i = 0; i < t.vertices.size(); ++i)
	{
		for (int const n : neighbours[i])
		{
			std::vector<int> intersection;
			std::set_intersection(neighbours[i].begin(), neighbours[i].end(),
				neighbours[n].begin(), neighbours[n].end(),
				std::back_inserter(intersection));

			for (auto const skeletonNeighbour : intersection)
			{
				skeletonEdges[i].insert(Triangulation::Edge( n, skeletonNeighbour ));
			}
		}
	}
	neighbours.clear();

	for (int it = 0; it < iterations; ++it)
	{
		UpdateSkeleton(t, skeletonEdges, lambda);
		UpdateSkeleton(t, skeletonEdges, mu);
	}
}

void UpdateSkeleton(Triangulation& t, std::map<int, std::set<Triangulation::Edge>>& skeletonEdges, double factor)
{
	std::vector<R3Graph::R3Point> meshBias(t.vertices.size(), { 0.,0.,0. });
	for (int i = 0; i < t.vertices.size(); ++i)
	{
		R3Graph::R3Point sumPoint;
		double sumLenght = 0.;
		for (auto const& edge : skeletonEdges[i])
		{
			int v1 = edge.vertIdx[0], v2 = edge.vertIdx[1];
			double lenght = (t.vertices[v1].point - t.vertices[v2].point).length();
			R3Graph::R3Point middle = (t.vertices[v1].point + t.vertices[v2].point) * 0.5;
			sumPoint += middle * lenght;
			sumLenght += lenght;
		}
		double oppLenght = 1 / sumLenght;
		R3Graph::R3Point baseMassCenter = sumPoint * oppLenght;
		meshBias[i] = baseMassCenter - t.vertices[i].point;
	}

	for (int i = 0; i < t.vertices.size(); ++i)
	{
		t.vertices[i].point += meshBias[i] * factor;
	}
}

void NormalsUpdate(Triangulation& triangulation)
{
	for (int i = 0; i < triangulation.triangles.size(); ++i)
	{
		Triangulation::Triangle& t = triangulation.triangles.at(i);

		R3Graph::R3Point& p0 = triangulation.vertices.at(t.indices[0]).point;
		R3Graph::R3Point& p1 = triangulation.vertices.at(t.indices[1]).point;
		R3Graph::R3Point& p2 = triangulation.vertices.at(t.indices[2]).point;

		R3Graph::R3Vector v1 = p1 - p0; v1.normalize();
		R3Graph::R3Vector v2 = p2 - p0; v2.normalize();
		R3Graph::R3Vector n = v1.vectorProduct(v2);
		n.normalize();
		if (t.Normal.scalarProduct(n) < 0)
			n *= (-1.);
		t.Normal = n;
	}
}