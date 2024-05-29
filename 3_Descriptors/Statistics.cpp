#include <iostream>
#include <fstream>
#include <cmath>
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"

using namespace std;
using namespace cv;

void drawHistogram(const Mat&, int);
int* calcolaIstogramma(Mat, int);

double imageMean(Mat);
double imageSigma(Mat);
double imageSkewness(Mat);
double imageKurtosis(Mat);

int main()
{
	string imagePath = "../images/Immagini Bitmap TEST/grayscale.bmp";

	Mat image;

	// Controllo che l'immagine sia stata effettivamente aperta
	image = imread(imagePath, IMREAD_GRAYSCALE);

	if (image.empty()) {
		cout << "Errore caricamento immagine." << endl;
		return -1;
	}
	cout << "Immagine letta correttamente." << endl;

	Mat noise(image.size(), CV_8UC1);
	randn(noise, 20, 10);

	add(image, noise, image);

	imshow("Immagine", image);
	waitKey(0);

	int bin = 16;
	int* hist2 = calcolaIstogramma(image, bin);


	int channel = 0;
	int histSize = 16;
	float grayranges[] = { 0, 256 };
	const float* range[] = { grayranges };

	Mat hist;
	calcHist(&image, 1, &channel, Mat(), hist, 1, &histSize, range, true, false);

	cout << "\n Opencv \n";
	for (int i = 0; i < histSize; i++)
	{
		cout << cvRound(hist.at<float>(i)) << " ";
	}

	// Normalizzazione dell'istogramma 0-1
	normalize(hist, hist, 0, 1, NORM_MINMAX, -1, Mat());
	drawHistogram(hist, histSize);


	cout << "media: " << imageMean(hist) << endl;
	cout << "std: " << imageSigma(hist) << endl;
	cout << "skew: " << imageSkewness(hist) << endl;
	cout << "kurt: " << imageKurtosis(hist) << endl;
	return 0;
}

int* calcolaIstogramma(Mat image, int bin) {
	int* v = new int[bin];
	for (int i = 0; i < bin; i++)
		v[i] = 0;
	int bin_dim = 256 / bin;

	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			int bidx = image.at<uchar>(i, j) / bin_dim;
			v[bidx]++;
		}
	}
	cout << "\n Codice custom \n";
	for (int i = 0; i < bin; i++)
		cout << v[i] << " ";

	return v;
}




void drawHistogram(const Mat& hist, int histSize) {
	// Dimensioni immagine istogramma
	int hist_w = 512;
	int hist_h = 512;
	// Larghezza colonna singolo bin
	int bin_w = cvRound((double)hist_w / histSize);
	// Immagine contenente l'istogramma
	Mat histImage(hist_h, hist_w, CV_8UC1, Scalar(0));

	for (int i = 0; i < histSize; i++)
	{
		// Estraggo i singoli valori dell'istograma e il scalo per poterli rappresentare
		int value = cvRound(hist.at<float>(i) * hist_h);

		rectangle(histImage,
			Point(i * bin_w, hist_h),
			Point(i * bin_w + bin_w, hist_h - value),
			Scalar::all(255),
			FILLED
		);
	}
	cout << endl;

	namedWindow("Istogramma", WINDOW_AUTOSIZE);
	imshow("Istogramma", histImage);
	waitKey(0);
}




double imageMean(Mat hist)
{
	double numEle = 0;
	double sum = 0.0;

	for (int i = 0; i < hist.rows; i++)
	{
		sum = sum + i * hist.at<float>(i);
		numEle = numEle + hist.at<float>(i);
	}

	return sum / numEle;
}

double imageSigma(Mat hist)
{
	double mean = imageMean(hist);
	double numEle = 0;
	double stdDev = 0.0;
	for (int i = 0; i < hist.rows; i++)
	{
		stdDev = stdDev + hist.at<float>(i) * pow((double)i - mean, 2);
		numEle = numEle + hist.at<float>(i);
	}

	stdDev = sqrt(stdDev / numEle);

	return stdDev;
}

double imageSkewness(Mat hist)
{
	double mean = imageMean(hist);
	double std = imageSigma(hist);
	double numEle = 0;

	double moment3 = 0.0;
	double skew = 0.0;

	for (int i = 0; i < hist.rows; i++)
	{
		moment3 = moment3 + pow((double)i - mean, 3) * hist.at<float>(i);
		numEle = numEle + hist.at<float>(i);
	}

	moment3 = moment3 / numEle;
	skew = moment3 / pow(std, 3.0);

	return skew;
}

double imageKurtosis(Mat hist)
{
	double mean = imageMean(hist);
	double std = imageSigma(hist);
	double numEle = 0;

	double moment4 = 0.0;
	double kurt = 0.0;

	for (int i = 0; i < hist.rows; i++)
	{
		moment4 = moment4 + pow((double)i - mean, 4) * hist.at<float>(i);
		numEle = numEle + hist.at<float>(i);
	}

	moment4 = moment4 / numEle;
	kurt = moment4 / pow(std, 4.0) - 3;

	return kurt;
}