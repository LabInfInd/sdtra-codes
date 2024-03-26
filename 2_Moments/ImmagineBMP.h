#pragma once
#include <iostream>
#include <fstream>

using namespace std;

class ImmagineBMP
{
private:
	int rows, columns;
	unsigned char** mat;

public:
	ImmagineBMP(string);
	~ImmagineBMP();
	ImmagineBMP(const ImmagineBMP&);

	int getRows() { return rows; }
	int getColums() { return columns; }
	unsigned char getPixVal(int i, int j) { return mat[i][j]; }
	void trasformaImmagine();
};

