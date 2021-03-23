#pragma once 
#include "Triangulation.h"

std::vector<std::set<Triangulation::Triangle>> Segmentation(Triangulation&, std::map<int, std::set<int>>&);

void Fill(std::map<int, std::set<int>>&, std::vector<Triangulation::Triangle>&, std::map<int, std::set<int>>&);

std::vector<std::set<Triangulation::Triangle>> DequeAlgo(std::vector<Triangulation::Triangle>&, 
														std::map<int, std::set<int>>&);