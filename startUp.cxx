#include "itkImage.h"
#include "itkImageFileReader.h"
//#include "itkSTLMeshIO.h"
#include "voxelset.h"
#include "Triangulation.h"
#include "R3Graph.h"
#include <iostream>
#include <fstream>

double VoxelDensity(const Voxel&, short*);
Voxel SearchSeed(short*, int, VoxelBox&);
void WriteStlFile(const Triangulation&, std::ofstream&);

int main(int argc, char* argv[])
{
	if (0) // 6 necessary arguments
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " -i head_cta.nii -m head_cta_segm.nii -t threshold -o out.nii" << std::endl;
		return EXIT_FAILURE;
	}

	const char* inputFileName = "";
	char* outputFileName = "";
	int threshold = 0;

	for (int i = 1; i < argc; ++i) {
		auto parametr = std::string(argv[i]);
		if (parametr == "-i") {
			inputFileName = argv[++i];
			continue;
		}
		else if (parametr == "-t") {
			threshold = std::stod(argv[++i]);
			continue;
		}
		else if (parametr == "-o") {
			outputFileName = argv[++i];
			continue;
		}
		else 
		{
			//std::cout << "Incorrect label " << argv[i] << std::endl;
		}
	}
	constexpr unsigned int Dimension = 3;
	using PixelType = short;
	using ImageType = itk::Image<PixelType, Dimension>;
	using ReaderType = itk::ImageFileReader<ImageType>;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName(inputFileName);
	reader->Update();
	auto image = reader->GetOutput();
	auto pointer = image->GetBufferPointer();

	auto region = image->GetBufferedRegion();
	auto size = region.GetSize();
	int xMax = size[0], yMax = size[1], MaxSlices = size[2];
	// 512, 512, 389
	auto scale = image->GetSpacing();
	auto x_sc = scale[0], y_sc = scale[1], z_sc = scale[2];

	VoxelBox voxelBoxOfImage(Voxel (0, (0,0)), xMax, yMax, MaxSlices); 
	Voxel seed = SearchSeed(pointer, threshold, voxelBoxOfImage);
	VoxelSet voxelSet;

	detectVoxelSetFromCta(
		(*VoxelDensity),
		threshold,
		voxelBoxOfImage,
		seed,
		pointer,
		voxelSet);
	
	Triangulation triangulation;
	
	computeTriangulationOfVoxelSet_MY(
		pointer, 
		threshold,
		triangulation,
		voxelSet,
		{ 0,0,0 },
		x_sc, y_sc, z_sc
	);

	std::ofstream out;
	std::string FileName = "C:\\Users\\owchi\\source\\repos\\TEST\\bin\\IsoSurface.stl";
	out.open(FileName);
	if(out.is_open())
		WriteStlFile(triangulation, out);
	out.close();

	return EXIT_SUCCESS;
}

double VoxelDensity(const Voxel& v, short* p) {// поменять тип функции на short
	auto VoxelPointer = p + (v.point.x + v.point.y * 512 + v.slice * 512 * 512);
	return *VoxelPointer;
}

Voxel SearchSeed(short* pointer, int threshold, VoxelBox& voxelBox) {
	Voxel seed;
	int num_seed = 0;
	int xMax = voxelBox.width, yMax = voxelBox.depth, MaxSlices = voxelBox.height;
//#pragma omp parallel for 
	for (int k = MaxSlices * 3 / 8; k < MaxSlices * 4 / 8; ++k) {
		for (int i = xMax * 4 / 10; i < xMax * 5 / 10; ++i) {
			for (int j = yMax * 4 / 10; j < yMax * 5 / 10; ++j) {
				auto voxel = pointer + (i + j * xMax + k * xMax * yMax);
				if (*voxel > threshold) {
				
					seed = { k, {i, j} };
					return seed;
				}
			}
		}
	}
	return seed;
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