#include <iostream>
#include <fstream>
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cmath>

#include "ImmagineBMP.h"

#define PI 3.14159265

using namespace std;
using namespace cv;

void openCVMoments(Mat);
void customMoments(ImmagineBMP);

double* centroidCoordinates(ImmagineBMP);
void axisInfo(ImmagineBMP, Mat);

int main()
{
	string imagePath = "../images/Immagini Momenti/lesion.bmp";
	string imageRotPath = "../images/Immagini Momenti/lesion_rotated.bmp";
	string imageTransPath = "../images/Immagini Momenti/lesion_translated.bmp";

	string imageScaledPath = "../images/Immagini Momenti/lesion_scaled.bmp";
	string imageRotTransPath = "../images/Immagini Momenti/lesion_rotated_translated.bmp";
	string imageRotTransScaledPath = "../images/Immagini Momenti/lesion_rotated_translated_scaled.bmp";
	string imageFlippedPath = "../images/Immagini Momenti/lesion_flipped.bmp";

	double* centroid;


	ImmagineBMP image(imagePath);
	//image.trasformaImmagine();
	centroid = centroidCoordinates(image);

	Mat imgOpen = imread(imagePath, CV_8UC1);
	Mat imgVisualization;

	// Converto l'immagine in scala di grigi a colori per mostrare sovrapposizione ellisse e assi sull'immagine
	cvtColor(imgOpen, imgVisualization, COLOR_GRAY2RGB);
	circle(imgVisualization, Point((int)centroid[0], (int)(centroid[1])), 10, Scalar(255, 0, 0));
	imshow("Image", imgVisualization);
	waitKey(0);

	axisInfo(image, imgVisualization);

	// Calcolo momenti immagine originale
	cout << imagePath << " Moments" << endl;
	customMoments(image);
	openCVMoments(imgOpen);


	ImmagineBMP image_r(imageRotPath);
	//image_r.trasformaImmagine();
	Mat imgOpen_r = imread(imageRotPath, CV_8UC1);

	imshow("Image Rotated", imgOpen_r);
	waitKey(0);

	// Calcolo momenti immagine ruotata
	cout << imageRotPath << " Moments" << endl;
	customMoments(image_r);
	openCVMoments(imgOpen_r);

	ImmagineBMP image_trans(imageTransPath);
	//image_trans.trasformaImmagine();
	Mat imgOpen_trans = imread(imageTransPath, CV_8UC1);

	imshow("Image Translated", imgOpen_trans);
	waitKey(0);

	cout << imageTransPath << " Moments" << endl;
	// Calcolo momenti immagine traslata
	customMoments(image_trans);
	openCVMoments(imgOpen_trans);


	ImmagineBMP image_scaled(imageScaledPath);
	//image_scaled.trasformaImmagine();
	Mat imgOpen_scaled = imread(imageScaledPath, CV_8UC1);

	imshow("Image Scaled", imgOpen_scaled);
	waitKey(0);

	cout << imageScaledPath << " Moments" << endl;
	// Calcolo momenti immagine scalata
	customMoments(image_scaled);
	openCVMoments(imgOpen_scaled);

	ImmagineBMP image_rot_trans(imageRotTransPath);
	//image_scaled.trasformaImmagine();
	Mat imgOpen_rot_trans = imread(imageRotTransPath, CV_8UC1);

	imshow("Image Rotated and Translated", imgOpen_rot_trans);
	waitKey(0);

	cout << imageRotTransPath << " Moments" << endl;
	// Calcolo momenti immagine traslata e ruotata
	customMoments(image_rot_trans);
	openCVMoments(imgOpen_rot_trans);

	ImmagineBMP image_rot_trans_scaled(imageRotTransScaledPath);
	//image_scaled.trasformaImmagine();
	Mat imgOpen_rot_trans_scaled = imread(imageRotTransScaledPath, CV_8UC1);

	imshow("Image Rotated, Translated, Scaled", imgOpen_rot_trans_scaled);
	waitKey(0);

	cout << imageRotTransScaledPath << " Moments" << endl;
	// Calcolo momenti immagine traslata, ruotata e scalata
	customMoments(image_rot_trans_scaled);
	openCVMoments(imgOpen_rot_trans_scaled);

	ImmagineBMP image_flipped(imageFlippedPath);
	//image_scaled.trasformaImmagine();
	Mat imgOpen_flipped = imread(imageFlippedPath, CV_8UC1);

	imshow("Image Flipped", imgOpen_flipped);
	waitKey(0);

	cout << imageFlippedPath << " Moments" << endl;
	// Calcolo momenti immagine flippata
	customMoments(image_flipped);
	openCVMoments(imgOpen_flipped);

#ifdef _WIN64
	system("pause");
#endif

	return 0;
}

void openCVMoments(Mat image)
{
	double humm[7];
	Moments moment = moments(image, false);
	HuMoments(moment, humm);


	//cout << "OpenCV Raw Moments:\t" << moment.m00 << "\t" << moment.m01 << "\t" << moment.m10 << "\t"
	//	<< moment.m11 << "\t" << moment.m02 << "\t" << moment.m20 <<
	//	"\t" << moment.m12 << "\t" << moment.m21 << "\t" << moment.m03 << "\t" << moment.m30 << endl;

	//cout << "OpenCV Central Moments:\t" << moment.mu11 << "\t" << moment.mu02 << "\t" << moment.mu20 <<
	//	"\t" << moment.mu12 << "\t" << moment.mu21 << "\t" << moment.mu03 << "\t" << moment.mu30 << endl;

	//cout << "OpenCV Norm. Central Moments:\t" << moment.nu11 << "\t" << moment.nu02 << "\t" << moment.nu20 <<
	//	"\t" << moment.nu12 << "\t" << moment.nu21 << "\t" << moment.nu03 << "\t" << moment.nu30 << endl;

	// Log scale hu moments 
	for (int i = 0; i < 7; i++) {
		humm[i] = -1 * copysign(1.0, humm[i]) * log10(abs(humm[i]));
	}

	cout << "OpenCV Hu Moments:\t" << humm[0] << "\t" << humm[1] << "\t" << humm[2] <<
		"\t" << humm[3] << "\t" << humm[4] << "\t" << humm[5] << "\t" << humm[6] << endl;

}

double* centroidCoordinates(ImmagineBMP image)
{
	// Passo 1: Calcolare i momenti semplici (raw moments)
	double m00 = 0.0, m01 = 0.0, m10 = 0.0;

	for (int i = 0; i < image.getRows(); i++)
		for (int j = 0; j < image.getColums(); j++)
		{
			m00 += (double)image.getPixVal(i, j);
			m01 += (double)i * image.getPixVal(i, j);
			m10 += (double)j * image.getPixVal(i, j);
		}

	// Centrodi
	double x_sign = 0.0, y_sign = 0.0;
	x_sign = m10 / m00;
	y_sign = m01 / m00;

	double* centroid = new double[2];
	centroid[0] = x_sign;
	centroid[1] = y_sign;
	return centroid;
}

void axisInfo(ImmagineBMP image, Mat drawingImage)
{
	double m00 = 0.0, m01 = 0.0, m10 = 0.0;

	for (int i = 0; i < image.getRows(); i++)
		for (int j = 0; j < image.getColums(); j++)
		{
			m00 += (double)image.getPixVal(i, j);
			m01 += (double)i * image.getPixVal(i, j);
			m10 += (double)j * image.getPixVal(i, j);
		}

	//cout << "Custom Raw Moments:\t" << m00 << "\t" << m01 << "\t" << m10 << "\t"
	//	<< m11 << "\t" << m02 << "\t" << m20 << 
	//	"\t" << m12 << "\t" << m21 << "\t" << m03 << "\t" << m30 << endl;

	// Passo 2: Calcolo i Momenti Centrali
	double mu00 = 0.0, mu01 = 0.0, mu10 = 0.0,
		mu11 = 0.0, mu02 = 0.0, mu20 = 0.0;

	// Centrodi
	double x_sign = 0.0, y_sign = 0.0;
	x_sign = m10 / m00;
	y_sign = m01 / m00;

	// By definition
	mu00 = m00;
	mu01 = 0;
	mu10 = 0;


	for (int i = 0; i < image.getRows(); i++)
		for (int j = 0; j < image.getColums(); j++)
		{
			double i_val = (i - y_sign);
			double j_val = (j - x_sign);

			mu11 += (double)i_val * j_val * image.getPixVal(i, j);
			mu02 += (double)i_val * i_val * image.getPixVal(i, j);
			mu20 += (double)j_val * j_val * image.getPixVal(i, j);
		}

	double I1 = (0.5 * (mu20 + mu02)) + (0.5 * sqrt(4 * pow(mu11, 2) + pow((mu20 - mu02), 2)));
	double I2 = (0.5 * (mu20 + mu02)) - (0.5 * sqrt(4 * pow(mu11, 2) + pow((mu20 - mu02), 2)));
	double theta = 0.5 * atan(2 * mu11 / (mu20 - mu02));

	double temp;
	if (mu20 < mu02)
	{
		temp = I1;
		I1 = I2;
		I2 = temp;
	}

	double a1 = 2 * sqrt(I1 / mu00);
	double a2 = 2 * sqrt(I2 / mu00);

	double angle = 0.0;
	angle = theta * 180 / PI;


	ellipse(drawingImage, Point((int)x_sign, (int)y_sign), Size((int)a1, (int)a2), angle, 0, 360, Scalar(0, 0, 255));
	imshow("Ellipse", drawingImage);
	waitKey(0);
	Point start, end, start1, end1;

	// NB: Cos e Sin vogliono l'angolo in radianti
	start = Point((int)(x_sign + a1 * cos(theta)), (int)(y_sign + a1 * sin(theta)));
	end = Point((int)(x_sign - a1 * cos(theta)), (int)(y_sign - a1 * sin(theta)));

	start1 = Point((int)(x_sign + a2 * sin(theta)), (int)(y_sign - a2 * cos(theta)));
	end1 = Point((int)(x_sign - a2 * sin(theta)), (int)(y_sign + a2 * cos(theta)));



	line(drawingImage, start, end, Scalar(0, 255, 0));

	line(drawingImage, start1, end1, Scalar(255, 0, 0));

	imshow("Ellipse Axis", drawingImage);
	waitKey(0);

}

void customMoments(ImmagineBMP image)
{
	// Passo 1: Calcolare i momenti semplici (raw moments)
	double m00 = 0.0, m01 = 0.0, m10 = 0.0,
		m11 = 0.0, m02 = 0.0, m20 = 0.0,
		m30 = 0.0, m03 = 0.0, m21 = 0.0, m12 = 0.0;

	for (int i = 0; i < image.getRows(); i++)
		for (int j = 0; j < image.getColums(); j++)
		{
			m00 += (double)image.getPixVal(i, j);
			m01 += (double)i * image.getPixVal(i, j);
			m10 += (double)j * image.getPixVal(i, j);
			m11 += (double)i * j * image.getPixVal(i, j);
			m02 += (double)i * i * image.getPixVal(i, j);
			m20 += (double)j * j * image.getPixVal(i, j);
			m12 += (double)i * i * j * image.getPixVal(i, j);
			m21 += (double)i * j * j * image.getPixVal(i, j);
			m03 += (double)i * i * i * image.getPixVal(i, j);
			m30 += (double)j * j * j * image.getPixVal(i, j);
		}

	//cout << "Custom Raw Moments:\t" << m00 << "\t" << m01 << "\t" << m10 << "\t"
	//	<< m11 << "\t" << m02 << "\t" << m20 << 
	//	"\t" << m12 << "\t" << m21 << "\t" << m03 << "\t" << m30 << endl;

	// Passo 2: Calcolo i Momenti Centrali
	double mu00 = 0.0, mu01 = 0.0, mu10 = 0.0,
		mu11 = 0.0, mu02 = 0.0, mu20 = 0.0,
		mu30 = 0.0, mu03 = 0.0, mu21 = 0.0, mu12 = 0.0;

	// Centrodi
	double x_sign = 0.0, y_sign = 0.0;
	x_sign = m10 / m00;
	y_sign = m01 / m00;

	// By definition
	mu00 = m00;
	mu01 = 0;
	mu10 = 0;


	for (int i = 0; i < image.getRows(); i++)
		for (int j = 0; j < image.getColums(); j++)
		{
			double i_val = (i - y_sign);
			double j_val = (j - x_sign);

			mu11 += (double)i_val * j_val * image.getPixVal(i, j);
			mu02 += (double)i_val * i_val * image.getPixVal(i, j);
			mu20 += (double)j_val * j_val * image.getPixVal(i, j);
			mu12 += (double)i_val * i_val * j_val * image.getPixVal(i, j);
			mu21 += (double)i_val * j_val * j_val * image.getPixVal(i, j);
			mu03 += (double)i_val * i_val * i_val * image.getPixVal(i, j);
			mu30 += (double)j_val * j_val * j_val * image.getPixVal(i, j);
		}

	//cout << "Custom Central Moments:\t" << mu11 << "\t" << mu02 << "\t" << mu20 <<
	//	"\t" << mu12 << "\t" << mu21 << "\t" << mu03 << "\t" << mu30 << endl;

	// Passo 3: Calcolo dei Momenti Centrali Normalizzati
	double ni11 = 0.0, ni02 = 0.0, ni20 = 0.0,
		ni30 = 0.0, ni03 = 0.0, ni21 = 0.0, ni12 = 0.0;

	ni11 = mu11 / pow(m00, 2);
	ni02 = mu02 / pow(m00, 2);
	ni20 = mu20 / pow(m00, 2);
	ni30 = mu30 / pow(m00, 2.5);
	ni03 = mu03 / pow(m00, 2.5);
	ni21 = mu21 / pow(m00, 2.5);
	ni12 = mu12 / pow(m00, 2.5);

	//cout << "Custom Norm. Central Moments:\t" << ni11 << "\t" << ni02 << "\t" << ni20 <<
	//	"\t" << ni12 << "\t" << ni21 << "\t" << ni03 << "\t" << ni30 << endl;

	// Passo 4: Hu Moments
	double I1 = 0.0, I2 = 0.0, I3 = 0.0, I4 = 0.0, I5 = 0.0, I6 = 0.0, I7 = 0.0;
	I1 = ni20 + ni02;
	I2 = pow((ni20 - ni02), 2) + 4 * ni11 * ni11;
	I3 = pow((ni30 - 3 * ni12), 2) + pow((3 * ni21 - ni03), 2);
	I4 = pow((ni30 + ni12), 2) + pow((ni21 + ni03), 2);
	I5 = (ni30 - 3 * ni12) * (ni30 + ni12) * (pow((ni30 + ni12), 2) - 3 * pow((ni21 + ni03), 2)) +
		(3 * ni21 - ni03) * (ni21 + ni03) * (3 * pow((ni30 + ni12), 2) - pow((ni21 + ni03), 2));
	I6 = (ni20 - ni02) * (pow((ni30 + ni12), 2) - pow((ni21 + ni03), 2)) + 4 * ni11 * (ni30 + ni12) * (ni21 + ni03);
	I7 = (3 * ni21 - ni03) * (ni30 + ni12) * (pow((ni30 + ni12), 2) - 3 * pow((ni21 + ni03), 2)) -
		(ni30 - 3 * ni12) * (ni21 + ni03) * (3 * pow((ni30 + ni12), 2) - pow((ni21 + ni03), 2));

	// Note that hu[0] is not comparable in magnitude as hu[6].
	// We can use use a log transform given below to bring them in the same range.

	I1 = -1 * copysign(1.0, I1) * log10(abs(I1));
	I2 = -1 * copysign(1.0, I2) * log10(abs(I2));
	I3 = -1 * copysign(1.0, I3) * log10(abs(I3));
	I4 = -1 * copysign(1.0, I4) * log10(abs(I4));
	I5 = -1 * copysign(1.0, I5) * log10(abs(I5));
	I6 = -1 * copysign(1.0, I6) * log10(abs(I6));
	I7 = -1 * copysign(1.0, I7) * log10(abs(I7));


	cout << "Custom Hu Moments:\t" << I1 << "\t" << I2 << "\t" << I3 <<
		"\t" << I4 << "\t" << I5 << "\t" << I6 << "\t" << I7 << endl;

}
