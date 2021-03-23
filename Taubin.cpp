#include "Taubin.h"

void Taubin(Triangulation& t, std::map<int, std::set<int>>& neighbours, double lambda, double mu, int iterations)
{
	for (int it = 0; it < iterations; ++it)
	{
		Update(t, neighbours, lambda);
		Update(t, neighbours, mu);
	}
	
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