#pragma once

#include "Triangulation.h"

void Taubin(Triangulation&, std::map<int, std::set<int>>&, double, double, int);

void Update(Triangulation&, std::map<int, std::set<int>>&, double);