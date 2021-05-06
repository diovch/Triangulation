#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRecursiveGaussianImageFilter.h"
//#include "itkSTLMeshIO.h"
#include "voxelset.h"
#include "Triangulation.h"
#include "R3Graph.h"
#include "Taubin.h"
#include "Segmentation.h"
#include "STL.h"
#include <iostream>
#include <fstream>
#include <algorithm>

short VoxelDensity(const Voxel&, short*);
Voxel SearchSeed(short*, int, VoxelBox&);
void WriteStlASCII(const Triangulation&, std::ofstream&);

int main(int argc, char* argv[])
{
	if (0) // 6 necessary arguments
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " -i head_cta.nii -s sigma -t threshold" << std::endl;
		return EXIT_FAILURE;
	}

	const char* inputFileName = "";
	char* outputFileName = "";
	short threshold = 0;
	double sigma = 0.;

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
		else if (parametr == "-s")
		{
			sigma = std::atof(argv[++i]);
			continue;
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

	if (1)
	{
		reader->SetFileName(inputFileName);
		reader->Update();
		using FilterType = itk::RecursiveGaussianImageFilter<ImageType, ImageType>;
		FilterType::Pointer smoothFilterX = FilterType::New();
		FilterType::Pointer smoothFilterY = FilterType::New();
		FilterType::Pointer smoothFilterZ = FilterType::New();
		smoothFilterX->SetDirection(0);
		smoothFilterY->SetDirection(1);
		smoothFilterZ->SetDirection(2);
		smoothFilterX->SetOrder(itk::GaussianOrderEnum::ZeroOrder);
		smoothFilterY->SetOrder(itk::GaussianOrderEnum::ZeroOrder);
		smoothFilterZ->SetOrder(itk::GaussianOrderEnum::ZeroOrder);
		smoothFilterX->SetNormalizeAcrossScale(true);
		smoothFilterY->SetNormalizeAcrossScale(true);
		smoothFilterZ->SetNormalizeAcrossScale(true);
		smoothFilterX->SetInput(reader->GetOutput());
		smoothFilterY->SetInput(smoothFilterX->GetOutput());
		smoothFilterZ->SetInput(smoothFilterY->GetOutput());
		smoothFilterX->SetSigma(sigma); smoothFilterY->SetSigma(sigma); smoothFilterZ->SetSigma(sigma);
		smoothFilterZ->Update();

		auto image = smoothFilterZ->GetOutput();
	}
	
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
	
	if(0)
		FillVoids(voxelSet);

	ContourSegmentation(voxelSet, seed);
	
	Triangulation triangulation;
	std::map<int, std::set<int>> VertexNeighbours;
	
	if (0)
		computeTriangulationOfVoxelSet_MY(
			pointer,
			threshold,
			triangulation,
			voxelSet,
			{ 0,0,0 },
			x_sc, y_sc, z_sc
		);
	else
		Triangulate_Custom(
			VertexNeighbours,
			triangulation,
			voxelSet,
			{ 0,0,0 },
			x_sc, y_sc, z_sc
		);
	
	

	if (0)
		Taubin(triangulation, VertexNeighbours, 0.33, -0.331, 25);
	//else
	//	triangulation.taubinSmoothing(1, 0.33, 0.331, false);
	
	//std::vector<std::set<Triangulation::Triangle>> s = Segmentation(triangulation, VerxteNeighbours);

	std::string filename = "IsoSurface.stl";
	WriteStlASCII(triangulation, filename);
	

	return EXIT_SUCCESS;
}

short VoxelDensity(const Voxel& v, short* p) {// �������� ��� ������� �� short
	auto VoxelPointer = p + (v.point.x + v.point.y * 512 + v.slice * 512 * 512);
	return *VoxelPointer;
}

Voxel SearchSeed(short* pointer, int threshold, VoxelBox& voxelBox) {
	Voxel seed;
	int num_seed = 0;
	int xMax = voxelBox.width, yMax = voxelBox.depth, MaxSlices = voxelBox.height;
	
	int aortaX = 278, aortaY = 239, aortaZ = 436;
	auto voxel = pointer + (aortaX + aortaY * xMax + aortaZ * xMax * yMax);
	if (*voxel > threshold) 
	{
		seed = { aortaZ, {aortaX, aortaY} };
		return seed;
	}

	for (int k = MaxSlices * 3 / 8; k < MaxSlices * 4 / 8; ++k) {
		for (int i = xMax * 4 / 10; i < xMax * 5 / 10; ++i) {
			for (int j = yMax * 4 / 10; j < yMax * 5 / 10; ++j) {
	//for (int k = 0; k < MaxSlices; ++k) {
	//	for (int i = 0; i < xMax; ++i) {
	//		for (int j = 0; j < yMax; ++j) {
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

