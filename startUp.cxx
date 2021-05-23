#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "voxelset.h"
#include "Triangulation.h"
#include "R3Graph.h"
#include "Taubin.h"
#include "Segmentation.h"
#include "STL.h"


Voxel SearchSeed(short*, unsigned char* , unsigned char , int, VoxelBox&);
void WriteStlASCII(const Triangulation&, std::ofstream&);

int main(int argc, char* argv[])
{
	if (0) // 6 necessary arguments
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " -i inputFileName.nii -m maskFileName.nii -l maskLabel -s sigma -t threshold -o outputFileName.stl" << std::endl;
		return EXIT_FAILURE;
	}

	const char* inputFileName = "";
	const char* maskFileName = "";
	unsigned char maskLabel{};
	char* outputFileName = "";
	short threshold = 0;
	double sigma = 0.;

	for (int i = 1; i < argc; ++i) 
	{
		auto parameter = std::string(argv[i]);
		if (parameter == "-i") 
		{
			inputFileName = argv[++i];
			continue;
		}
		else if (parameter == "-t") 
		{
			threshold = std::stod(argv[++i]);
			continue;
		}
		else if (parameter == "-m") 
		{
			maskFileName = argv[++i];
			continue;
		}
		else if (parameter == "-o") 
		{
			outputFileName = argv[++i];
			continue;
		}
		else if (parameter == "-s")
		{
			sigma = std::atof(argv[++i]);
			continue;
		}
		else if (parameter == "-l")
		{
			maskLabel = std::atoi(argv[++i]);
			continue;
		}
	}
	constexpr unsigned int Dimension = 3;
	using PixelType = short;

	using ImageType = itk::Image<PixelType, Dimension>;
	using MaskType = itk::Image<unsigned char, Dimension>;
	using ReaderType = itk::ImageFileReader<ImageType>;
	using ReaderMaskType = itk::ImageFileReader<MaskType>;
	ReaderType::Pointer reader = ReaderType::New();
	ReaderMaskType::Pointer mask_reader = ReaderMaskType::New();

	reader->SetFileName(inputFileName);
	reader->Update();
	mask_reader->SetFileName(maskFileName);
	mask_reader->Update();

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

	auto pointer = image->GetBufferPointer();
	auto mask_pointer = (mask_reader->GetOutput())->GetBufferPointer();

	auto region = image->GetBufferedRegion();
	auto size = region.GetSize();
	int xMax = size[0], yMax = size[1], MaxSlices = size[2];
	// 512, 512, 389
	auto scale = image->GetSpacing();
	auto x_sc = scale[0], y_sc = scale[1], z_sc = scale[2];

	VoxelBox voxelBoxOfImage(Voxel (0, (0,0)), xMax, yMax, MaxSlices); 
	Voxel seed = SearchSeed(pointer, mask_pointer, maskLabel, threshold, voxelBoxOfImage);
	VoxelSet voxelSet;
	
	detectVoxelSetFromCta(
		threshold,
		voxelBoxOfImage,
		seed,
		pointer,
		mask_pointer,
		maskLabel,
		voxelSet);
	
	if(0)
		FillVoids(voxelSet);
	
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
	

	if (1)
		Taubin(triangulation, VertexNeighbours, 0.33, -0.331, 15);

	std::string filename = "";

	if (outputFileName == "")
		filename = "IsoSurface.stl";
	else
		filename = outputFileName;

	WriteStlASCII(triangulation, filename);

	return EXIT_SUCCESS;
}




Voxel SearchSeed(short* pointer, unsigned char* mask_pointer, unsigned char maskLabel, int threshold, VoxelBox& voxelBox) 
{
	Voxel seed;
	int num_seed = 0;
	int xMax = voxelBox.width, yMax = voxelBox.depth, MaxSlices = voxelBox.height;
	
	//for (int k = MaxSlices * 1 / 10; k < MaxSlices * 9 / 10; ++k)
	//{
	//	for (int i = xMax * 1 / 10; i < xMax * 9 / 10; ++i)
	//	{
	//		for (int j = yMax * 1 / 10; j < yMax * 9 / 10; ++j)
	//		{
	for (int k = 0; k < MaxSlices; ++k) {
		for (int i = 0; i < xMax; ++i) {
			for (int j = 0; j < yMax; ++j) {
				auto voxel = pointer + (i + j * xMax + k * xMax * yMax);
				auto isROI = mask_pointer + (i + j * xMax + k * xMax * yMax);
				if (*voxel > threshold && *isROI == maskLabel)
				{
					seed = { k, {i, j} };
					return seed;
				}
			}
		}
	}
	return seed;
}

