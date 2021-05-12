#pragma once
#include "Triangulation.h"
#include <string>

void WriteBinarySTL(const Triangulation&, std::string&);


void WriteStlFile(const Triangulation& triangulation, std::ofstream& out);

void WriteStlFile(const Triangulation& triangulation, std::ofstream& out, std::vector<std::set<Triangulation::Triangle>>& s);