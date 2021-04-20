#include "STL.h"
#include "R3Graph.h"
#include <fstream>
//TODO Tune WriteBinarySTL()
void WriteSTLBinary(const Triangulation& triangulation, std::string& filename)
{
	std::ofstream out;
	filename = "C:\\Users\\owchi\\source\\repos\\TEST\\bin\\" + filename;
	out.open(filename, std::ofstream::binary);
	if (out.is_open())
	{

		char header[80];
		std::memset(header, 0, 80);
		out.write(header, 80);

		int numTriang = int(triangulation.triangles.size());
		out.write((const char*)(&numTriang), 4);


		for (int i = 0; i < triangulation.triangles.size(); ++i)
		{
			Triangulation::Triangle t = triangulation.triangles.at(i);
			
			R3Graph::R3Point p[4];
			p[0] = triangulation.vertices.at(t.indices[0]).point;
			p[1] = triangulation.vertices.at(t.indices[1]).point;
			p[2] = triangulation.vertices.at(t.indices[2]).point;
			p[3] = t.Normal;

			for (int i = 0; i < 4; ++i)
			{
				float vec[3];
				vec[0] = (float) p[i].x;
				vec[1] = (float) p[i].y;
				vec[2] = (float) p[i].z;
				out.write((const char*) vec, 12);
			}
		
			unsigned short dummy = 48;
			out.write((const char*)&dummy, 2);
		}
	}
	out.close();
}

void WriteStlASCII(const Triangulation& triangulation, std::string& filename)
{
	std::ofstream out;
	filename = "C:\\Users\\owchi\\source\\repos\\TEST\\bin\\" + filename;
	out.open(filename);
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

