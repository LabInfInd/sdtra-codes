#include "ImmagineBMP.h"

ImmagineBMP::ImmagineBMP(string path)
{
	ifstream f(path, ios::in | ios::binary);

	// Mi muovo di 18 byte dall'inizio del file per leggere direttamente le informazioni 
	// di larghezza e altezza dell'immagine
	f.seekg(18, ios::beg);
	f.read((char*)&columns, sizeof(columns));
	f.read((char*)&rows, sizeof(rows));

	f.seekg(1078, ios::beg);


	mat = new unsigned char* [rows];
	for (int i = 0; i < rows; i++)
		mat[i] = new unsigned char[columns];

	int r = 4 - (columns % 4);
	unsigned char temp = 0;

	for (int i = rows - 1; i >= 0; i--) {
		for (int j = 0; j < columns; j++)
			f.read((char*)&mat[i][j], sizeof(mat[i][j]));

		if (columns % 4 != 0)
			for (int k = 0; k < r; k++)
				f.read((char*)&temp, sizeof(temp));
	}

	f.close();
}

ImmagineBMP::~ImmagineBMP() {

	for (int i = 0; i < rows; i++)
		delete mat[i];
	delete mat;
}

ImmagineBMP::ImmagineBMP(const ImmagineBMP& im)
{
	rows = im.rows;
	columns = im.columns;

	mat = new unsigned char* [rows];
	for (int i = 0; i < rows; i++)
		mat[i] = new unsigned char[columns];

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			mat[i][j] = im.mat[i][j];

}

void ImmagineBMP::trasformaImmagine()
{
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
		{
			if (mat[i][j] > 0)
			{
				mat[i][j] = 1;
			}
		}

}
