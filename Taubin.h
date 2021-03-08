#pragma once

#include "Triangulation.h"

void Taubin(Triangulation&, double, double, int);

void Update(Triangulation&, std::map<int, std::set<int>>&, double);