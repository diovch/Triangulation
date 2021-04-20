#pragma once
#include "Triangulation.h"
#include <string>

void WriteSTLBinary(const Triangulation&, std::string&);

void WriteStlASCII(const Triangulation& triangulation, std::string& filename);