#pragma once

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

using namespace std;

#ifndef IMMAGINE_H_
#define IMMAGINE_H_

struct palette {

    unsigned char rgbBlue;
    unsigned char rgbGreen;
    unsigned char rgbRed;
    unsigned char rgbReserved;

};

class immagine {

    unsigned short int bfType;
    int bfSize;
    unsigned short int bfReserved1;
    unsigned short int bfReserved2;
    int bfOffBits;

    int biSize;
    int biWidth;
    int biHeight;

    unsigned short int   biPlanes;
    unsigned short int   biBitCount;
    int  biCompression;
    int  biSizeImage;
    int  biXPelsPerMeter;
    int  biYPelsPerMeter;
    int  biClrUsed;
    int  biClrImportant;

    struct palette pal[256];

    char nome[100];
    unsigned char** mat;

public:

    immagine() {}
    immagine(char*);
    immagine(const immagine&);
    ~immagine();

    int getWidth() { return biWidth; }
    int getHeight() { return biHeight; }
    unsigned char** getMatrix() { return mat; }

    void setDim(int w, int h) { biWidth = w; biHeight = h; }
    void setMatrix(unsigned char** mat) { this->mat = mat; }

    void salvaImmagine(char*);
    void salvaImmagine(char*, int);

    unsigned char** conv(int, int**);
    unsigned char** conv(int, float**, bool);

    void creaScacchiera(int, int);
};

#endif /* IMMAGINE_H_ */

void compressioneRLE(ofstream&, immagine);
void decompressioneRLE(ifstream&, immagine&);

