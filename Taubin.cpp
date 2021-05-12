#include "Taubin.h"
#include "R3Graph.h"
// TODO: Implement another Taubin smoothing
void Taubin(Triangulation& t, std::map<int, std::set<int>>& neighbours, double lambda, double mu, int iterations)
{
	for (int it = 0; it < iterations; ++it)
	{
		Update(t, neighbours, lambda);
		Update(t, neighbours, mu);
	}

	NormalsUpdate(t);
	
	return;
}

void Update(Triangulation& t, std::map<int, std::set<int>>& neighbours, double factor)
{
	for (int i = 0; i < t.vertices.size(); ++i)
	{
		double weight = 1. / neighbours[i].size();
		R3Graph::R3Vector average_vector = { 0., 0., 0. };

		for (const auto ind : neighbours[i])
			average_vector += (t.vertices[ind].point - t.vertices[i].point) * weight;


		t.vertices[i].point += average_vector * factor;
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