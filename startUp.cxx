#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "voxelset.h"
#include "Triangulation.h"
#include "R3Graph.h"
#include "Taubin.h"
#include "Segmentation.h"
#include "STL.h"

int main(int argc, char* argv[])
{
	if (0) // 6 necessary arguments
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " -i head_cta.nii -m head_cta_segm.nii -s sigma -t threshold" << std::endl;
		return EXIT_FAILURE;
	}

	const char* inputFileName = "";
	char* outputFileName = "";
	short threshold = 0;
	const char* maskFileName = "";
	double sigma = 0.;

	for (int i = 1; i < argc; ++i) {
		auto parameter = std::string(argv[i]);
		if (parameter == "-i") {
			inputFileName = argv[++i];
			continue;
		}
		else if (parameter == "-t") {
			threshold = std::stod(argv[++i]);
			continue;
		}
		else if (parameter == "-m") {
			maskFileName = argv[++i];
			continue;
		}
		else if (parameter == "-o") {
			outputFileName = argv[++i];
			continue;
		}
		else if (parameter == "-s")
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

	
	using MaskType = itk::Image<unsigned char, Dimension>;
	using ReaderMaskType = itk::ImageFileReader<MaskType>;
	ReaderMaskType::Pointer mask_reader = ReaderMaskType::New();

	reader->SetFileName(inputFileName);
	reader->Update();
	auto image = reader->GetOutput();

	mask_reader->SetFileName(maskFileName);
	mask_reader->Update();
	auto mask_image = mask_reader->GetOutput();


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
	auto mask_pointer = mask_image->GetBufferPointer();

	auto region = image->GetBufferedRegion();
	auto size = region.GetSize();
	int xMax = size[0], yMax = size[1], MaxSlices = size[2];
	// 512, 512, 389
	auto scale = image->GetSpacing();
	auto x_sc = scale[0], y_sc = scale[1], z_sc = scale[2];

	VoxelBox voxelBoxOfImage(Voxel (0, (0,0)), xMax, yMax, MaxSlices); 
	Voxel seed = SearchSeed(pointer, mask_pointer, threshold, voxelBoxOfImage);
	VoxelSet voxelSet;
	
	detectVoxelSetFromCta(
		threshold,
		voxelBoxOfImage,
		seed,
		pointer,
		mask_pointer,
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
		Taubin(triangulation, VertexNeighbours, 0.33, -0.331, 35);

	std::string filename = "IsoSurface.stl";
	WriteStlASCII(triangulation, filename);

	return EXIT_SUCCESS;
}
