#include <iostream>
#include <fstream>
#include <list>
#include <iomanip>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

void creaMatriceCooccorrenze(Mat, int**, int, vector<int>, bool);
void updateRanges(int*, int*, int*, int*, vector<int>);
void features(int**, int, double*, bool);
int sommaMatrice(int**, int);
double mean(double*, int);
double standardDeviation(double*, int, double);
double calculateMCC(double**, int);
double calculateHXY2(double*, double*, int);
double calculateHXY1(double**, double*, double*, int);
double calculateEntropy(double*, int);
void calculate_p_plus(double**, int, double*);
void calculate_p_minus(double**, int, double*);
void saveMatrixToFile(int**, int, string);

int main()
{
	string filename = "../images/Immagini Bitmap TEST/scacchiera.bmp";
	Mat inputImage = imread(filename, IMREAD_GRAYSCALE);

	if (inputImage.data == NULL)
	{
		// Immagine non letta correttamente;
		return 1;
	}

	// Se l'immagine è stata letta correttamente, vado avanti
	//int direction = 0; // 45, 90, 135
	list<vector<int>> offsets;
	int dist = 1;
	vector<int> offset = { 0, dist };
	offsets.push_back(offset);
	/*offset = { dist, dist };
	offsets.push_back(offset);
	offset = { -dist, 0 };
	offsets.push_back(offset);
	offset = { dist, -dist };
	offsets.push_back(offset);*/
	bool symmetric = true;

	int grayLevels = 16;

	for (int i = 0; i < inputImage.rows; i++)
		for (int j = 0; j < inputImage.cols; j++)
		{
			int temp = floor((double)inputImage.at<uchar>(i, j) / 255 * (grayLevels));
			if (temp > grayLevels - 1)
			{
				temp = grayLevels - 1;
			}
			inputImage.at<uchar>(i, j) = (unsigned char)temp;
		}

	// Alloco matrice delle co-occorrenze
	int** matrix;
	matrix = new int* [grayLevels];
	for (int i = 0; i < grayLevels; i++)
		matrix[i] = new int[grayLevels];

	vector<double*> res;

	while (!offsets.empty())
	{
		vector<int> offset = offsets.front();
		creaMatriceCooccorrenze(inputImage, matrix, grayLevels, offset, symmetric);
		saveMatrixToFile(matrix, grayLevels, "mat1.txt");

		double* feat_vect = new double[14];
		for (int i = 0; i < 14; i++)
			feat_vect[i] = 0.0;

		features(matrix, grayLevels, feat_vect, symmetric);

		res.push_back(feat_vect);
		offsets.pop_front();
	}

	for (int i = 0; i < res.size(); i++)
	{
		for (int j = 0; j < 14; j++)
		{
			cout << fixed;
			cout << setprecision(5);
			cout << (res[i])[j] << "\t";
		}
		cout << endl;
	}

	return 0;
}

void creaMatriceCooccorrenze(Mat image, int** matrix, int grey, vector<int> offset, bool symmetric)
{
	for (int i = 0; i < grey; i++)
		for (int j = 0; j < grey; j++)
			matrix[i][j] = 0;

	int x_min = 0;
	int y_min = 0;
	int x_max = image.cols;
	int y_max = image.rows;


	updateRanges(&x_min, &x_max, &y_min, &y_max, offset);

	for (int i = 0; i < grey; i++)
	{
		for (int j = 0; j < grey; j++)
		{
			for (int x = x_min; x < x_max; x++)
			{
				for (int y = y_min; y < y_max; y++)
				{
					if ((int)(image.at<uchar>(Point(x, y))) == i &&
						(int)(image.at<uchar>(Point(x + offset[1], y - offset[0]))) == j)
						matrix[i][j] += 1;
				}
			}
		}
	}

	if (symmetric)
	{
		x_min = 0;
		y_min = 0;
		x_max = image.cols;
		y_max = image.rows;

		vector<int> temp = { -offset[0], -offset[1] };
		updateRanges(&x_min, &x_max, &y_min, &y_max, temp);

		for (int i = 0; i < grey; i++)
		{
			for (int j = 0; j < grey; j++)
			{
				for (int x = x_min; x < x_max; x++)
				{
					for (int y = y_min; y < y_max; y++)
					{
						if ((int)(image.at<uchar>(Point(x, y))) == i &&
							(int)(image.at<uchar>(Point(x + temp[1], y - temp[0]))) == j)
							matrix[i][j] += 1;
					}
				}
			}
		}
	}
}

void updateRanges(int* x_min, int* x_max, int* y_min, int* y_max, vector<int> offset)
{
	if (offset[1] < 0)
		*x_min = *x_min - offset[1];
	else
		*x_max = *x_max - offset[1];

	if (offset[0] > 0)
		*y_max = *y_max - offset[0];
	else
		*y_min = *y_min - offset[0];
}

void features(int** matrix, int grey, double* feat, bool symmetric)
{
	// Calcolo Matrice cooccorrenze normalizzata
	double** normMatrix = new double* [grey];
	for (int i = 0; i < grey; i++)
		normMatrix[i] = new double[grey];

	double somma = (double)sommaMatrice(matrix, grey);

	for (int i = 0; i < grey; i++)
		for (int j = 0; j < grey; j++)
			normMatrix[i][j] = (double)matrix[i][j] / somma;

	// Marginal row prob (p_x(i))
	// Marginal col prob (p_y(j))
	double* mrp = new double[grey];
	double* mcp = new double[grey];

	for (int i = 0; i < grey; i++)
	{
		mrp[i] = 0.0;
		mcp[i] = 0.0;
	}

	for (int i = 0; i < grey; i++)
	{
		for (int j = 0; j < grey; j++)
		{
			mrp[i] += normMatrix[i][j];
			mcp[j] += normMatrix[i][j];
		}
	}

	double mu_x = mean(mrp, grey);
	double mu_y = mean(mcp, grey);
	double std_x = standardDeviation(mrp, grey, mu_x);
	double std_y = standardDeviation(mcp, grey, mu_y);

	// Calcolo p_x+y(k)
	double* p_plus = new double[grey * 2 - 1];
	for (int i = 0; i < grey * 2 - 1; i++)
		p_plus[i] = 0.0;
	calculate_p_plus(normMatrix, grey, p_plus);
	// Calcolo p_x-y(k)
	double* p_minus = new double[grey];
	for (int i = 0; i < grey; i++)
		p_minus[i] = 0.0;
	calculate_p_minus(normMatrix, grey, p_minus);

	double hx = 0.0;
	hx = calculateEntropy(mrp, grey);
	double hy = 0.0;
	hy = calculateEntropy(mcp, grey);

	double hxy1 = 0.0;
	hxy1 = calculateHXY1(normMatrix, mrp, mcp, grey);
	double hxy2 = 0.0;
	hxy2 = calculateHXY2(mrp, mcp, grey);

	double** Q;
	Q = new double* [grey];
	for (int i = 0; i < grey; i++)
		Q[i] = new double[grey];

	for (int i = 0; i < grey; i++)
	{
		for (int j = 0; j < grey; j++)
		{
			Q[i][j] = 0;
		}
	}

	for (int i = 0; i < grey; i++)
	{
		for (int j = 0; j < grey; j++)
		{
			for (int k = 0; k < grey; k++)
			{
				Q[i][j] += (normMatrix[i][k] * normMatrix[j][k]) / ((mrp[i] * mcp[k]) + DBL_EPSILON);
			}
		}
	}


	// 1. Angular second moment
	double angular_sm = 0.0;
	// 2. Contrast
	double contrast = 0.0;
	// 3. Correlation
	double correlation = 0.0;
	// 4. Variance
	double variance = 0.0;
	// 5. Inverse Diference Moment
	double idm = 0.0;
	// 6. Sum average
	double sum_avg = 0.0;
	// 7. Sum variance
	double sum_var = 0.0;
	// 8. Sum entropy
	double sum_entr = 0.0;
	// 9. Entropy
	double entr = 0.0;
	// 10. Diff variance
	double diff_var = 0.0;
	// 11. Diff entropy
	double diff_entr = 0.0;
	// 12. 13. Information Measures of Correlation
	double imc1 = 0.0;
	double imc2 = 0.0;
	// 14. Maximal Correlation Coefficient
	double mcc = 0.0;

	for (int i = 0; i < grey; i++)
	{
		for (int j = 0; j < grey; j++)
		{
			angular_sm += pow(normMatrix[i][j], 2.0);
			contrast += pow((i - j), 2.0) * normMatrix[i][j];
			correlation += (double)i * (double)j * normMatrix[i][j];
			variance += (pow((double)i - mu_x, 2.0) * normMatrix[i][j]);
			idm += normMatrix[i][j] / (1 + pow((i + j), 2));
			entr += normMatrix[i][j] * log(normMatrix[i][j] + DBL_EPSILON);
		}
	}

	for (int i = 0; i < grey * 2 - 1; i++)
	{
		sum_avg += i * p_plus[i];
		sum_entr += p_plus[i] * log(p_plus[i] + DBL_EPSILON);

	}

	sum_entr = -sum_entr;
	entr = -entr;

	for (int i = 0; i < grey * 2 - 1; i++)
	{
		sum_var += pow(((double)i - sum_avg), 2.0) * p_plus[i];
	}

	double  diff_avg = 0.0;
	for (int i = 0; i < grey; i++)
	{
		diff_avg += i * p_minus[i];
		diff_entr += p_minus[i] * log(p_minus[i] + DBL_EPSILON);
	}

	diff_entr = -diff_entr;

	for (int i = 0; i < grey; i++)
	{
		diff_var += pow(((double)i - diff_avg), 2.0) * p_minus[i];
	}

	imc1 = (entr - hxy1) / max(hx, hy);
	imc2 = sqrt(1 - exp(-2 * (hxy2 - entr)));

	mcc = calculateMCC(Q, grey);


	if (std_x == 0)
		std_x = 1;
	if (std_y == 0)
		std_y = 1;

	correlation = (correlation - (mu_x * mu_y)) / (std_x * std_y);

	feat[0] = angular_sm;
	feat[1] = contrast;
	feat[2] = correlation;
	feat[3] = variance;
	feat[4] = idm;
	feat[5] = sum_avg;
	feat[6] = sum_var;
	feat[7] = sum_entr;
	feat[8] = entr;
	feat[9] = diff_var;
	feat[10] = diff_entr;
	feat[11] = imc1;
	feat[12] = imc2;
	feat[13] = mcc;
}

int sommaMatrice(int** mat, int grey)
{
	int somma = 0;
	for (int i = 0; i < grey; i++)
		for (int j = 0; j < grey; j++)
			somma += mat[i][j];
	return somma;
}

double mean(double* probDistr, int grey)
{
	double mean = 0.0;
	for (int i = 0; i < grey; i++)
		mean += (double)i * probDistr[i];

	return mean;
}

void calculate_p_plus(double** matrix, int grey, double* plus)
{
	for (int i = 0; i < grey; i++)
		for (int j = 0; j < grey; j++)
			plus[i + j] += matrix[i][j];
}

void calculate_p_minus(double** matrix, int grey, double* minus)
{
	for (int i = 0; i < grey; i++)
		for (int j = 0; j < grey; j++)
			minus[abs(i - j)] += matrix[i][j];
}

double standardDeviation(double* probDistr, int grey, double mean)
{
	double sigma = 0;
	for (int i = 0; i < grey; i++)
		sigma += probDistr[i] * pow(((double)i - mean), 2.0);

	return sqrt(sigma);
}

double calculateEntropy(double* arr, int dim)
{
	double val = 0.0;
	for (int i = 0; i < dim; i++)
		val += arr[i] * log(arr[i] + DBL_EPSILON);
	return -val;
}

double calculateHXY1(double** mat, double* px, double* py, int dim)
{
	double val = 0.0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			val += mat[i][j] * log((px[i] * py[j]) + DBL_EPSILON);
		}
	}
	return -val;
}

double calculateHXY2(double* px, double* py, int dim)
{
	double val = 0.0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			val += px[i] * py[j] * log((px[i] * py[j]) + DBL_EPSILON);
		}
	}
	return -val;
}

double calculateMCC(double** mat, int dim)
{
	double* mat2 = new double[dim * dim];
	int k = 0;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
		{
			mat2[k] = mat[i][j];
			k = k + 1;
		}
	Mat matrix(dim, dim, CV_64FC1, mat2);
	Mat eigs;
	eigenNonSymmetric(matrix, eigs, noArray());
	cv::sort(eigs, eigs, SORT_EVERY_COLUMN + SORT_DESCENDING);

	// Prendo il secondo autovalore
	double val = eigs.at<double>(1, 0);

	return sqrt(val);
}

void saveMatrixToFile(int** mat, int dim, string filename)
{
	fstream file;
	file.open(filename, ios::out);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			file << mat[i][j] << "\t";
		}
		file << endl;
	}
	file.close();
}

