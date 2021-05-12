#include "STL.h"
#include "R3Graph.h"
#include <fstream>
//TODO Tune WriteBinarySTL()
void WriteBinarySTL(const Triangulation& triangulation, std::string& filename)
{
	std::ofstream out;
	filename = "C:\\Users\\owchi\\source\\repos\\TEST\\bin\\" + filename;
	out.open(filename);
	if (out.is_open())
	{
		out << "header" << std::endl;
		out << triangulation.triangles.size() << std::endl;
		out << "foreach triangle" << std::endl;
		for (int i = 0; i < triangulation.triangles.size(); ++i)
		{
			Triangulation::Triangle t = triangulation.triangles.at(i);

			R3Graph::R3Point p0 = triangulation.vertices.at(t.indices[0]).point;
			R3Graph::R3Point p1 = triangulation.vertices.at(t.indices[1]).point;
			R3Graph::R3Point p2 = triangulation.vertices.at(t.indices[2]).point;

			out << t.Normal.x << " " << t.Normal.y << " " << t.Normal.z << std::endl;

			out << p0.x << " " << p0.y << " " << p0.z << std::endl;
			out << p1.x << " " << p1.y << " " << p1.z << std::endl;
			out << p2.x << " " << p2.y << " " << p2.z << std::endl;

			out << 0 << std::endl;
		}
		out << "end";
	}
	out.close();
}

void WriteStlFile(const Triangulation& triangulation, std::ofstream& out)
{
	out << "solid iso" << std::endl;

	for (int i = 0; i < triangulation.triangles.size(); ++i) {
		out << "    " << "facet normal ";

		Triangulation::Triangle t = triangulation.triangles.at(i);

		R3Graph::R3Point p0 = triangulation.vertices.at(t.indices[0]).point;
		R3Graph::R3Point p1 = triangulation.vertices.at(t.indices[1]).point;
		R3Graph::R3Point p2 = triangulation.vertices.at(t.indices[2]).point;

		out << t.Normal.x << " " << t.Normal.y << " " << t.Normal.z << std::endl;
		out << "    " << "outer loop" << std::endl;

		out << "    " << "vertex " << p0.x << " " << p0.y << " " << p0.z << std::endl;
		out << "    " << "vertex " << p1.x << " " << p1.y << " " << p1.z << std::endl;
		out << "    " << "vertex " << p2.x << " " << p2.y << " " << p2.z << std::endl;

		out << "    " << "endloop" << std::endl;
		out << "    " << "endfacet" << std::endl;
	}

	out << "endsolid iso";
}

