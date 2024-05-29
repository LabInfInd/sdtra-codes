#include <iostream>
#include <vector>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkOpenCVImageBridge.h>
#include <itkOrientImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkViewImage.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>

#include <itkConfidenceConnectedImageFilter.h>
#include <itkChangeInformationImageFilter.h>

#include <itkMaskImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <vtkMarchingCubes.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <itkRegionOfInterestImageFilter.h>
#include <vtkSTLWriter.h>
#include <vtkMassProperties.h>

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

using namespace std;
using namespace cv;
using namespace itk;

#if defined(__APPLE__) || defined(__MACH__)
#define PLATFORM_NAME "macos"
#else
#define PLATFORM_NAME "windows"
#endif

const int Dimension2D = 2;
const int Dimension3D = 3;
using PixelType = short;
using PixelTypeUC = unsigned char;
using ImageType = itk::Image<PixelType, Dimension2D>;
using ImageTypeUC = itk::Image<PixelTypeUC, Dimension2D>;
using SeriesType = itk::Image<PixelType, Dimension3D>;
using SeriesTypeUC = itk::Image<PixelTypeUC, Dimension3D>;

// Seme per RG
SeriesTypeUC::IndexType seed;

// Serie DICOM in UCHAR
SeriesTypeUC::Pointer serie;

// Maschera RG
SeriesTypeUC::Pointer mask;

// VOI TUMORE
SeriesTypeUC::Pointer voi;

void visualizzaDICOMSerie();
void visualizzaMaschera();
void regionGrowing();
void estraiVOI();
void filtraVOI();
void renderObjects();
cv::Mat estraiImmagine(int, SeriesTypeUC::Pointer);

void CallBackFuncAxial(int event, int x, int y, int flags, void* userdata)
{
    //USerdata contiene valore attuale e massimo
    if (strcmp(PLATFORM_NAME, "windows") == 0) {
        // WINDOWS SECTION
        // USerdata contiene valore attuale e massimo
        if (event == cv::EVENT_MOUSEWHEEL)
        {
            int maxVal = ((int*)userdata)[1];
            if (cv::getMouseWheelDelta(flags) > 0)
            {
                if (((int*)userdata)[0] < maxVal - 1)
                    ((int*)userdata)[0] += 1;
            }
            else
            {
                if (((int*)userdata)[0] > 0)
                    ((int*)userdata)[0] -= 1;
            }
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, serie);
            imshow("immagine assiale", img);
        }
        if (event == cv::EVENT_LBUTTONDOWN && flags == cv::EVENT_FLAG_SHIFTKEY + cv::EVENT_LBUTTONDOWN)
        {
            int val = ((int*)userdata)[0];
            seed[0] = x;
            seed[1] = y;
            seed[2] = val;
            return;
        }
    }
    else {
        if (event == cv::EVENT_LBUTTONDOWN && flags == cv::EVENT_FLAG_SHIFTKEY + cv::EVENT_LBUTTONDOWN)
        {
            int val = ((int*)userdata)[0];
            seed[0] = x;
            seed[1] = y;
            seed[2] = val;
            return;
        }
        if (event == cv::EVENT_LBUTTONDOWN)
        {
            int maxVal = ((int*)userdata)[1];
            if (((int*)userdata)[0] < maxVal - 1)
                ((int*)userdata)[0] += 1;
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, serie);
            imshow("immagine assiale", img);
        }
        if (event == cv::EVENT_RBUTTONDOWN)
        {
            int maxVal = ((int*)userdata)[1];
            if (((int*)userdata)[0] > 0)
                ((int*)userdata)[0] -= 1;
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, serie);
            imshow("immagine assiale", img);
        }
    }
}

void CallBackFuncSegmentation(int event, int x, int y, int flags, void* userdata)
{
    //USerdata contiene valore attuale e massimo
    if (strcmp(PLATFORM_NAME, "windows") == 0) {
        // WINDOWS SECTION
        // USerdata contiene valore attuale e massimo
        if (event == cv::EVENT_MOUSEWHEEL)
        {
            int maxVal = ((int*)userdata)[1];
            if (cv::getMouseWheelDelta(flags) > 0)
            {
                if (((int*)userdata)[0] < maxVal - 1)
                    ((int*)userdata)[0] += 1;
            }
            else
            {
                if (((int*)userdata)[0] > 0)
                    ((int*)userdata)[0] -= 1;
            }
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, mask);
            imshow("segmentazione", img);
        }
    }
    else {
        if (event == cv::EVENT_LBUTTONDOWN)
        {
            int maxVal = ((int*)userdata)[1];
            if (((int*)userdata)[0] < maxVal - 1)
                ((int*)userdata)[0] += 1;
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, mask);
            imshow("segmentazione", img);
        }
        if (event == cv::EVENT_RBUTTONDOWN)
        {
            int maxVal = ((int*)userdata)[1];
            if (((int*)userdata)[0] > 0)
                ((int*)userdata)[0] -= 1;
            int val = ((int*)userdata)[0];
            std::cout << "Slice: " << (val + 1) << " of " << maxVal << std::endl;
            Mat img = estraiImmagine(val, mask);
            imshow("segmentazione", img);
        }
    }
}

int main()
{
    // Lettura serie DICOM
    using ReaderType = itk::ImageFileReader<ImageType>;
    using SeriesReaderType = itk::ImageSeriesReader<SeriesType>;

    // Specifico il tipo di file da leggere (utile nel caso in cui non sia possibile dedurlo automaticamente)
    using ImageIOType = itk::GDCMImageIO;
    ImageIOType::Pointer dcmImageIO = ImageIOType::New();

    std::string seriesPath = "C:\\Users\\Brunetti\\Downloads\\anonimoComplete\\anonimoComplete";

    // Generiamo i nomi dei file
    using NameGeneratorType = itk::GDCMSeriesFileNames;
    NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
    nameGenerator->SetInputDirectory(seriesPath);
    nameGenerator->Update();

    std::vector<std::string> UIDs = nameGenerator->GetSeriesUIDs();

    SeriesReaderType::Pointer dicomSeriesReader = SeriesReaderType::New();
    dicomSeriesReader->SetFileNames(nameGenerator->GetFileNames(UIDs[0]));
    dicomSeriesReader->SetImageIO(dcmImageIO);
    dicomSeriesReader->Update();

    SeriesType::Pointer dcmSerie = dicomSeriesReader->GetOutput();

    // Orientare il volume
    using OrientImageFilter = itk::OrientImageFilter<SeriesType, SeriesType>;
    OrientImageFilter::Pointer filtroO = OrientImageFilter::New();
    filtroO->UseImageDirectionOn();
    filtroO->SetDesiredCoordinateOrientation(itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RAS);
    filtroO->SetInput(dicomSeriesReader->GetOutput());
    filtroO->Update();

    int level = -600;
    int larghezza = 1500;

    int min = level - (int)(larghezza / 2);
    int max = level + (int)(larghezza / 2);


    using IntensityFilterType = itk::IntensityWindowingImageFilter<SeriesType, SeriesTypeUC>;
    IntensityFilterType::Pointer filtro = IntensityFilterType::New();
    filtro->SetInput(filtroO->GetOutput());
    filtro->SetWindowMinimum(min);
    filtro->SetWindowMaximum(max);
    filtro->SetOutputMinimum(0);
    filtro->SetOutputMaximum(255);
    filtro->Update();

    serie = SeriesTypeUC::New();
    serie = filtro->GetOutput();

    visualizzaDICOMSerie();
    regionGrowing();
    visualizzaMaschera();
    estraiVOI();
    filtraVOI();
    renderObjects();

}

void visualizzaDICOMSerie() {
    int sliceNumber = 100;
    SeriesTypeUC::RegionType inputRegion = serie->GetLargestPossibleRegion();
    SeriesTypeUC::SizeType regionSize = inputRegion.GetSize();
    int numSlicesSagittal = regionSize[0];
    int numSlicesCoronal = regionSize[1];
    int numSlicesAxial = regionSize[2];

    int userD[2] = { sliceNumber , numSlicesAxial };
    cv::Mat img = estraiImmagine(sliceNumber, serie);
    imshow("immagine assiale", img);
    cv::setMouseCallback("immagine assiale", CallBackFuncAxial, userD);
    cv::waitKey();
}

cv::Mat estraiImmagine(int ind, SeriesTypeUC::Pointer ser) {

    SeriesTypeUC::RegionType inputRegion = ser->GetLargestPossibleRegion();
    SeriesTypeUC::SizeType regionSize = inputRegion.GetSize();
    int numSlicesSagittal = regionSize[0];
    int numSlicesCoronal = regionSize[1];
    int numSlicesAxial = regionSize[2];

    regionSize[2] = 0;
    SeriesTypeUC::IndexType start = inputRegion.GetIndex();
    start[2] = ind;

    SeriesTypeUC::RegionType desiredRegion;
    desiredRegion.SetSize(regionSize);
    desiredRegion.SetIndex(start);

    using FilterType = itk::ExtractImageFilter<SeriesTypeUC, ImageTypeUC >;
    FilterType::Pointer extractFilter = FilterType::New();
    extractFilter->SetExtractionRegion(desiredRegion);
    extractFilter->SetInput(ser);
    extractFilter->SetDirectionCollapseToIdentity();
    extractFilter->Update();

    return itk::OpenCVImageBridge::ITKImageToCVMat<ImageTypeUC>(extractFilter->GetOutput());
}

void regionGrowing()
{
    // input : SERIE
    // output: MASK

    // itk::ConfidenceConnectedImageFilter
    using FiltroRG = itk::ConfidenceConnectedImageFilter<SeriesTypeUC, SeriesTypeUC>;
    FiltroRG::Pointer filtro = FiltroRG::New();
    filtro->SetInitialNeighborhoodRadius(3);
    filtro->SetMultiplier(2.5);
    filtro->SetNumberOfIterations(1);
    filtro->SetReplaceValue(255);
    filtro->SetInput(serie);
    filtro->SetSeed(seed);
    filtro->Update();

    mask = filtro->GetOutput();
}

void visualizzaMaschera() {
    int sliceNumber = 100;
    SeriesTypeUC::RegionType inputRegion = mask->GetLargestPossibleRegion();
    SeriesTypeUC::SizeType regionSize = inputRegion.GetSize();
    int numSlicesSagittal = regionSize[0];
    int numSlicesCoronal = regionSize[1];
    int numSlicesAxial = regionSize[2];

    int userD[2] = { sliceNumber , numSlicesAxial };
    cv::Mat img = estraiImmagine(sliceNumber, mask);
    imshow("segmentazione", img);
    cv::setMouseCallback("segmentazione", CallBackFuncSegmentation, userD);
    cv::waitKey();
}

void renderObjects()
{
    // Costruiamo modello 3D dei polmoni
    using vtkConvert = itk::ImageToVTKImageFilter<SeriesTypeUC>;
    vtkConvert::Pointer convert = vtkConvert::New();
    convert->SetInput(mask);
    convert->Update();

    // Costruiamo modello 3D della lesione
    vtkConvert::Pointer convertLesion = vtkConvert::New();
    convertLesion->SetInput(voi);
    convertLesion->Update();

    vtkNew<vtkMarchingCubes> cube;
    cube->SetInputData(convert->GetOutput());
    cube->SetValue(0, 1);
    cube->Update();
    // Oggetto 3d cube->GetOutput();
    vtkNew<vtkSTLWriter> stlWriter;
    string filename = "polmoni.stl";
    stlWriter->SetFileName(filename.c_str());
    stlWriter->SetInputData(cube->GetOutput());
    stlWriter->Update();
    stlWriter->Write();

    vtkNew<vtkMarchingCubes> cubeLesion;
    cubeLesion->SetInputData(convertLesion->GetOutput());
    cubeLesion->SetValue(0, 2);
    cubeLesion->Update();
    // Oggetto 3d cubeLesion->GetOutput();
    filename = "lesione.stl";
    stlWriter->SetFileName(filename.c_str());
    stlWriter->SetInputData(cubeLesion->GetOutput());
    stlWriter->Update();
    stlWriter->Write();

    vtkNew<vtkMassProperties> props;
    props->SetInputConnection(cube->GetOutputPort());
    props->Update();
    std::cout << "Volume: " << props->GetVolume() << std::endl;
    std::cout << "Superficie: " << props->GetSurfaceArea() << std::endl;


    // Genero le proprietà necessarie al motore grafico del PC per mostrare graficamente l'oggetto 3d
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(cube->GetOutput());
    mapper->ScalarVisibilityOff();
    mapper->Update();

    vtkNew<vtkPolyDataMapper> mapperLesion;
    mapperLesion->SetInputData(cubeLesion->GetOutput());
    mapperLesion->ScalarVisibilityOff();
    mapperLesion->Update();

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0, 0, 1);
    actor->GetProperty()->SetOpacity(0.5);
    actor->ApplyProperties();

    vtkNew<vtkActor> actorLesion;
    actorLesion->SetMapper(mapperLesion);
    actorLesion->GetProperty()->SetColor(1, 0, 0);
    actorLesion->GetProperty()->SetOpacity(1);
    actorLesion->ApplyProperties();

    vtkNew<vtkRenderer> render;
    vtkNew<vtkRenderWindow> renWin;
    vtkNew<vtkRenderWindowInteractor> iren;

    renWin->AddRenderer(render);
    iren->SetRenderWindow(renWin);
    render->AddActor(actor);
    render->AddActor(actorLesion);
    render->SetBackground(0, 0, 0);

    renWin->Render();
    iren->Start();
}

void filtraVOI()
{
    using BinaryFilterType = itk::BinaryThresholdImageFilter<SeriesTypeUC, SeriesTypeUC>;
    BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();
    binaryFilter->SetInput(voi);
    binaryFilter->SetLowerThreshold(200);
    binaryFilter->SetUpperThreshold(255);
    binaryFilter->SetInsideValue(255);
    binaryFilter->SetOutsideValue(0);
    binaryFilter->Update();

    voi = binaryFilter->GetOutput();
}

void estraiVOI()
{
    /* ATTENZIONE: Utilizzare ExtractImageFilter al posto di RegionOfInterestImageFilter
    per mantenere nella VOI i riferimenti spaziali della serie di partenza (necessari per renderizzare la
    il volume nel posto giusto).

    RegionOfInterestImageFilter utile quando il processing da effettuare deve limitarsi ad una regione limitata
    nelle sue dimensioni.
    */
    using VOIExtractor = itk::ExtractImageFilter<SeriesTypeUC, SeriesTypeUC>;
    VOIExtractor::Pointer voiFilter = VOIExtractor::New();

    SeriesTypeUC::IndexType indexStart;
    indexStart[0] = 324;
    indexStart[1] = 307;
    indexStart[2] = 270;

    SeriesTypeUC::SizeType regionSize;
    regionSize[0] = 40;
    regionSize[1] = 32;
    regionSize[2] = 35;

    SeriesTypeUC::RegionType region;
    region.SetIndex(indexStart);
    region.SetSize(regionSize);

    voiFilter->SetInput(serie);
    voiFilter->SetExtractionRegion(region);
    voiFilter->Update();
    voi = voiFilter->GetOutput();
}
