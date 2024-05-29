#include <iostream>
#include <list>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;


// Funzione di regionGrowing: immagine, seme, soglia, tipo_connettività -> Restrituisce la maschera di segmentazione ottenuta
Mat regionGrowing4(Mat, vector<Point>, double);
void regionGrowing8(Mat, vector<Point>, double, Mat&);
double regionMeanValue(Mat, Mat);
void salvaImmagine(Mat, string);
void CallBackFunc(int event, int x, int y, int flags, void* userdata);


int main()
{
	// Aperuta Immagine
	string filename = "../images/liver.jpg";
	Mat image = imread(filename, CV_8UC1);

	// Variabili necessarie per region growing
	vector<Point> seeds;
	//seeds.push_back(Point(100,100));

	// TODO: Populate list of seed points from mous selection

		// Imposta funzione di callback per il click del mouse
	namedWindow("Seleziona seed");
	imshow("Seleziona seed", image);
	setMouseCallback("Seleziona seed", CallBackFunc, &seeds);
	waitKey(0);

	// Set threshold
	double threshold = (double)8 / 255;

	// Chimata alla funzione che fa region growing
	Mat masks;
	Mat blurred;
	blur(image, blurred, Size(5, 5));
	masks = regionGrowing4(blurred, seeds, threshold);
	regionGrowing8(blurred, seeds, threshold, masks);
	waitKey(0);

	cout << "Segmentazione terminata" << endl;

	// Modifica con operatori morfologici
	/*
	cv::MorphShapes {
	  cv::MORPH_RECT = 0,
	  cv::MORPH_CROSS = 1,
	  cv::MORPH_ELLIPSE = 2
	}
	*/

	Mat element = getStructuringElement(MORPH_CROSS, Size(7,7));
	Mat processed;
	erode(masks, processed, element);
	imshow("eroded", processed);
	waitKey(0);

	dilate(masks, processed, element);
	imshow("dilated", processed);
	waitKey(0);

	morphologyEx(masks, processed, MORPH_CLOSE, element, Point(-1, -1), 5);
	imshow("Closed", processed);
	waitKey(0);


	//salvaImmagine(mask, "./images/roi_m.bmp");
	return 0;
}

Mat regionGrowing4(Mat image, vector<Point> seeds, double thresh)
{
	list<Point> listaPunti;
	Mat mask = Mat::zeros(image.size(), image.type());
	for (int i = 0; i < seeds.size(); i++)
	{
		mask.at<uchar>(seeds.at(i)) = 255;
		listaPunti.push_back(seeds.at(i));
	}

	int counter = 0;

	Point neighbours[4];

	double region_mean = regionMeanValue(image, mask);

	Mat inputImage;
	cvtColor(image, inputImage, COLOR_GRAY2BGR);

	Mat maskImage;
	cvtColor(mask, maskImage, COLOR_GRAY2BGR);

	Mat blendedImage;
	addWeighted(inputImage, 0.5, maskImage, 0.5, 0.0, blendedImage);
	namedWindow("Segmentazione");
	imshow("Segmentazione", blendedImage);
	waitKey(1);

	// Eseguo region growing fino a quando non ho più punti da valutare
	while (!listaPunti.empty())
	{
		// Prendo il prossimo punto da valutare
		Point next = listaPunti.front();

		// Vedo se per questo punto ci sono punti nel vicinato da aggiungere alla regione
		// Caso 4-neighbors
		neighbours[0].x = next.x - 1;
		neighbours[0].y = next.y;

		neighbours[1].x = next.x + 1;
		neighbours[1].y = next.y;

		neighbours[2].x = next.x;
		neighbours[2].y = next.y - 1;

		neighbours[3].x = next.x;
		neighbours[3].y = next.y + 1;

		// Per ogni vicino, valuto se può appartenere o meno alla regione
		// In caso positivo, lo aggiungo alla maschera e ai prossimi punti da valutare

		for (int i = 0; i < 4; i++)
		{
			Rect range = Rect(0, 0, image.cols, image.rows);
			// 1. Devo prima valutare se il punto è ammissibile (le coordinate rientrano nelle dimensioni dell'immagine?)
			if (neighbours[i].inside(range))
			{
				// Devo valutare se il punto non appartiene gia alla maschera
				if (mask.at<uchar>(neighbours[i]) == 0)
				{
					// Calcolo la distanza
					double dist = abs((double)image.at<uchar>(neighbours[i]) - region_mean) / 255;
					if (dist < thresh)
					{
						// Aggiungo il punto alla maschera e alla lista dei punti da esaminare
						mask.at<uchar>(neighbours[i]) = 255;
						listaPunti.push_back(neighbours[i]);
					}
				}
			}
		}

		// Aggiorno la media della regione
		region_mean = regionMeanValue(image, mask);

		// Elimino il punto dalla lista
		listaPunti.pop_front();

		counter++;

		if (counter % 10000 == 0)
		{
			cvtColor(mask, maskImage, COLOR_GRAY2BGR);
			addWeighted(inputImage, 0.5, maskImage, 0.5, 0.0, blendedImage);
			imshow("Segmentazione", blendedImage);
			waitKey(1);
		}



	}
	return mask;
}

void regionGrowing8(Mat image, vector<Point> seeds, double thresh, Mat& masks)
{
	int counter = 0;
	list<Point> listaPunti;
	masks = Mat::zeros(image.size(), image.type());
	for (int i = 0; i < seeds.size(); i++)
	{
		masks.at<uchar>(seeds.at(i)) = 255;
		listaPunti.push_back(seeds.at(i));
	}

	vector<Point> neighbours;

	double region_mean;

	Mat inputImage;
	cvtColor(image, inputImage, COLOR_GRAY2BGR);
	
	Mat maskImage;
	cvtColor(masks, maskImage, COLOR_GRAY2BGR);

	Mat blendedImage;
	addWeighted(inputImage, 0.5, maskImage, 0.5, 0.0, blendedImage);
	namedWindow("Segmentazione");
	imshow("Segmentazione", blendedImage);
	waitKey(1);

	region_mean = regionMeanValue(image, masks);
	while (!listaPunti.empty())
	{
		// Aggiorno la media della regione
		region_mean = regionMeanValue(image, masks);
		counter = counter + 1;
		// Prendo il prossimo punto da valutare
		Point next = listaPunti.front();


		// Vedo se per questo punto ci sono punti nel vicinato da aggiungere alla regione
		// Caso 8-neighbors
		for (int x = next.x - 1; x < next.x + 2; x++)
		{
			for (int y = next.y - 1; y < next.y + 2; y++)
			{
				if (!(x == next.x && y == next.y))
				{
					neighbours.push_back(Point(x, y));
				}
			}
		}


		// Per ogni vicino, valuto se può appartenere o meno alla regione
		// In caso positivo, lo aggiungo alla maschera e ai prossimi punti da valutare

		for (int i = 0; i < neighbours.size(); i++)
		{
			Rect range = Rect(0, 0, image.cols, image.rows);
			// 1. Devo prima valutare se il punto è ammissibile (le coordinate rientrano nelle dimensioni dell'immagine?)
			if (neighbours[i].inside(range))
			{
				// Devo valutare se il punto non appartiene gia alla maschera
				if (masks.at<uchar>(neighbours[i]) == 0)
				{
					// Calcolo la distanza
					double dist = abs((double)image.at<uchar>(neighbours[i]) - region_mean) / 255;
					if (dist < thresh)
					{
						// Aggiungo il punto alla maschera e alla lista dei punti da esaminare
						masks.at<uchar>(neighbours[i]) = 255;
						listaPunti.push_back(neighbours[i]);
					}
				}
			}
		}

		// Elimino il punto dalla lista
		listaPunti.pop_front();
		neighbours.clear();

		if (counter % 100 == 0)
		{
			cvtColor(masks, maskImage, COLOR_GRAY2BGR);
			addWeighted(inputImage, 0.5, maskImage, 0.5, 0.0, blendedImage);
			imshow("Segmentazione", blendedImage);
			waitKey(1);
		}

	}

	Mat ch1, ch2, ch3; // declare three matrices 
	// "channels" is a vector of 3 Mat arrays:
	vector<Mat> channels(3);
	// split img:
	split(maskImage, channels);
	// get the channels (follow BGR order in OpenCV)
	channels[0] = Mat::zeros(maskImage.size(), CV_8UC1);
	channels[1] = Mat::zeros(maskImage.size(), CV_8UC1);
	// modify channel// then merge


	Mat image2;
	merge(channels, maskImage);
	addWeighted(inputImage, 0.5, maskImage, 0.5, 0.0, blendedImage);
	imshow("Segmentazione Finita", blendedImage);
	waitKey(1);
}

double regionMeanValue(Mat image, Mat mask)
{
	double mean = 0.0;
	int counter = 0;

	for (int i = 0; i < image.rows; i++)
	{
		for (int j = 0; j < image.cols; j++)
		{
			if (mask.at<uchar>(i, j) > 0)
			{
				mean = mean + image.at<uchar>(i, j);
				counter = counter + 1;
			}
		}
	}
	if (counter > 0)
	{
		mean = mean / counter;
	}
	else
		mean = 0;

	
	return mean;
}

void salvaImmagine(Mat image, string filename)
{
	imwrite(filename, image);

}

void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{
	vector<Point>* seeds = (vector<Point>*)userdata;
	if (event == EVENT_LBUTTONDOWN)
	{
		cout << "Selezionato il punto (" << x << ", " << y << ")" << endl;
		// cout << "Premere un tasto qualsiasi per proseguire." << endl;
		seeds->push_back(Point(x, y));
	}
}



