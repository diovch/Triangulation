#pragma once

#include "Triangulation.h"

void Taubin(Triangulation&, std::map<int, std::set<int>>&, double, double, int);

void Update(Triangulation&, std::map<int, std::set<int>>&, double);

void NormalsUpdate(Triangulation&);

void TaubinSkeleton(Triangulation&, std::map<int, std::set<int>>&, double, double, int);
void UpdateSkeleton(Triangulation& t, std::map<int, std::set<Triangulation::Edge>>& skeletonEdges, double factor);