#include <iostream>
#include "Immagine.h"

using namespace std;

int main()
{
    /*
     Esercizio 1:
     - aprire immagine bitmap in scala di grigi a scelta
     - stampare a video le dimensioni dell'immagine
     */

    //TODO: Aggiornare con percorso relativo o assoluto dell'immagine da aprire
    const char* imagePath = "../images/Immagini Bitmap TEST/scacchiera.bmp";
    immagine im1((char*)imagePath);
    int w = im1.getWidth();
    int h = im1.getHeight();

    cout << "Larghezza: " << w << endl;
    cout << "Altezza: " << h << endl;


    immagine im_s(im1);
    im_s.creaScacchiera(8, 16);

    /*
     Esercizio 2:
     - definire kernel 3x3 (K)
     - eseguire la convoluzione dell'immagine aperta con il kernel K
     - salvare l'immagine convoluta
     */

    int dim = 3;
    int** k;
    k = new int* [dim];
    for (int i = 0; i < dim; i++)
        k[i] = new int[dim];

    // Kernel sobel X
    k[0][0] = 1; k[0][1] = 0; k[0][2] = -1;
    k[1][0] = 2; k[1][1] = 0; k[1][2] = -2;
    k[2][0] = 1; k[2][1] = 0; k[2][2] = -1;

    // Salvataggio della convoluzione in res
    unsigned char** res = im1.conv(dim, k);
    immagine im2(im1); // Costruttore di copia
    int newh = h;
    int neww = w;
    im2.setDim(neww, newh);
    im2.setMatrix(res);

    im2.salvaImmagine((char*)"conv3_x.bmp");

    // Kernel sobel Y
    k[0][0] = 1; k[0][1] = 2; k[0][2] = 1;
    k[1][0] = 0; k[1][1] = 0; k[1][2] = 0;
    k[2][0] = -1; k[2][1] = -2; k[2][2] = -1;
    
    res = im1.conv(dim, k);

    immagine im3(im1); // Costruttore di copia
    im3.setDim(neww, newh);
    im3.setMatrix(res);

    im3.salvaImmagine((char*)"conv3_y.bmp");

    /*
     Esercizio 2:
     - definire kernel 5x5 (K)
     - eseguire la convoluzione dell'immagine aperta con il kernel K
     - salvare l'immagine convoluta
     */

    dim = 5;
    k = new int* [dim];
    for (int i = 0; i < dim; i++)
        k[i] = new int[dim];

    k[0][0] = -2; k[0][1] = -1; k[0][2] = 0; k[0][3] = 1; k[0][4] = 2;
    k[1][0] = -2; k[1][1] = -1; k[1][2] = 0; k[1][3] = 1; k[1][4] = 2;
    k[2][0] = -2; k[2][1] = -1; k[2][2] = 0; k[2][3] = 1; k[2][4] = 2;
    k[3][0] = -2; k[3][1] = -1; k[3][2] = 0; k[3][3] = 1; k[3][4] = 2;
    k[4][0] = -2; k[4][1] = -1; k[4][2] = 0; k[4][3] = 1; k[4][4] = 2;

    res = im1.conv(dim, k);
    immagine im4(im1);
    newh = h;
    neww = w;
    im4.setDim(neww, newh);
    im4.setMatrix(res);

    im4.salvaImmagine((char*)"convoluzione_dimk.bmp");

    /*
     Esercizio 3:
     - definire kernel media (K1)
     - eseguire la convoluzione dell'immagine aperta con il kernel K1
     - salvare l'immagine convoluta
     */

    dim = 3;
    float** k1;
    k1 = new float* [dim];
    for (int i = 0; i < dim; i++)
        k1[i] = new float[dim];

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            k1[i][j] = 1.0 / (dim * dim);

    // Chiamo la funzione con 0-padding
    res = im1.conv(dim, k1, true);
    immagine im5(im1);
    newh = h;
    neww = w;
    im5.setDim(neww, newh);
    im5.setMatrix(res);

    im5.salvaImmagine((char*)"convoluzione_blurk.bmp");

    // Chiamo la funzione con 0-padding
    res = im1.conv(dim, k1, true);
    immagine im6(im1);
    newh = h-dim+1;
    neww = w-dim+1;
    im6.setDim(neww, newh);
    im6.setMatrix(res);

    im6.salvaImmagine((char*)"convoluzione_blurk_no_pad.bmp");


    return 0;
}







