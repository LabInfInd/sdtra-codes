/*
    #include "immagine.h" in questo caso dal momento che nella directory il file .h si chiama Immagine.h, 
    è necessario includere come di seguito
*/
#include "Immagine.h" // sul nome del file non viene fatto il controllo case sensitive

immagine::immagine(char* filename) {

    ifstream f(filename, ios::in | ios::binary);

    f.read((char*)&bfType, sizeof(bfType));
    f.read((char*)&bfSize, sizeof(bfSize));
    f.read((char*)&bfReserved1, sizeof(bfReserved1));
    f.read((char*)&bfReserved2, sizeof(bfReserved2));
    f.read((char*)&bfOffBits, sizeof(bfOffBits));
    f.read((char*)&biSize, sizeof(biSize));
    f.read((char*)&biWidth, sizeof(biWidth));
    f.read((char*)&biHeight, sizeof(biHeight));
    f.read((char*)&biPlanes, sizeof(biPlanes));
    f.read((char*)&biBitCount, sizeof(biBitCount));
    f.read((char*)&biCompression, sizeof(biCompression));
    f.read((char*)&biSizeImage, sizeof(biSizeImage));
    f.read((char*)&biXPelsPerMeter, sizeof(biXPelsPerMeter));
    f.read((char*)&biYPelsPerMeter, sizeof(biYPelsPerMeter));
    f.read((char*)&biClrUsed, sizeof(biClrUsed));
    f.read((char*)&biClrImportant, sizeof(biClrImportant));

    f.read((char*)&pal, sizeof(pal));

    mat = new unsigned char* [biHeight];
    for (int i = 0; i < biHeight; i++)
        mat[i] = new unsigned char[biWidth];

    int r = 4 - (biWidth % 4);
    unsigned char temp = 0;

    for (int i = biHeight - 1; i >= 0; i--) {
        for (int j = 0; j < biWidth; j++)
            f.read((char*)&mat[i][j], sizeof(mat[i][j]));

        if (biWidth % 4 != 0)
            for (int k = 0; k < r; k++)
                f.read((char*)&temp, sizeof(temp));
    }

    f.close();

}

immagine::immagine(const immagine& im) {

    bfType = im.bfType;
    bfSize = im.bfSize;
    bfReserved1 = im.bfReserved1;
    bfReserved2 = im.bfReserved2;
    bfOffBits = im.bfOffBits;
    biSize = im.biSize;
    biWidth = im.biWidth;
    biHeight = im.biHeight;
    biPlanes = im.biPlanes;
    biBitCount = im.biBitCount;
    biCompression = im.biCompression;
    biSizeImage = im.biSizeImage;
    biXPelsPerMeter = im.biXPelsPerMeter;
    biYPelsPerMeter = im.biYPelsPerMeter;
    biClrUsed = im.biClrUsed;
    biClrImportant = im.biClrImportant;

    for (int i = 0; i < 256; i++)
        pal[i] = im.pal[i];

    this->mat = new unsigned char* [biHeight];
    for (int i = 0; i < biHeight; i++)
        this->mat[i] = new unsigned char[biWidth];

    for (int i = 0; i < biHeight; i++)
        for (int j = 0; j < biWidth; j++)
            this->mat[i][j] = im.mat[i][j];

}

immagine::~immagine() {

    for (int i = 0; i < biHeight; i++)
        delete mat[i];
    delete mat;

}

void immagine::salvaImmagine(char* name) {

    ofstream f(name, ios::out | ios::binary);

    f.write((char*)&bfType, sizeof(bfType));
    f.write((char*)&bfSize, sizeof(bfSize));
    f.write((char*)&bfReserved1, sizeof(bfReserved1));
    f.write((char*)&bfReserved2, sizeof(bfReserved2));
    f.write((char*)&bfOffBits, sizeof(bfOffBits));
    f.write((char*)&biSize, sizeof(biSize));
    f.write((char*)&biWidth, sizeof(biWidth));
    f.write((char*)&biHeight, sizeof(biHeight));
    f.write((char*)&biPlanes, sizeof(biPlanes));
    f.write((char*)&biBitCount, sizeof(biBitCount));
    f.write((char*)&biCompression, sizeof(biCompression));
    f.write((char*)&biSizeImage, sizeof(biSizeImage));
    f.write((char*)&biXPelsPerMeter, sizeof(biXPelsPerMeter));
    f.write((char*)&biYPelsPerMeter, sizeof(biYPelsPerMeter));
    f.write((char*)&biClrUsed, sizeof(biClrUsed));
    f.write((char*)&biClrImportant, sizeof(biClrImportant));

    f.write((char*)&pal, sizeof(pal));

    int r = 4 - (biWidth % 4);
    unsigned char temp = 0;

    for (int i = biHeight - 1; i >= 0; i--) {
        for (int j = 0; j < biWidth; j++)
            f.write((char*)&mat[i][j], sizeof(mat[i][j]));

        if (biWidth % 4 != 0)
            for (int k = 0; k < r; k++)
                f.write((char*)&temp, sizeof(temp));
    }
    f.flush();
    f.close();
}

void immagine::salvaImmagine(char* name, int val) {

    ofstream f(name, ios::out | ios::binary);

    f.write((char*)&bfType, sizeof(bfType));
    f.write((char*)&bfSize, sizeof(bfSize));
    f.write((char*)&bfReserved1, sizeof(bfReserved1));
    f.write((char*)&bfReserved2, sizeof(bfReserved2));
    f.write((char*)&bfOffBits, sizeof(bfOffBits));
    f.write((char*)&biSize, sizeof(biSize));
    f.write((char*)&biWidth, sizeof(biWidth));
    f.write((char*)&biHeight, sizeof(biHeight));
    f.write((char*)&biPlanes, sizeof(biPlanes));
    f.write((char*)&biBitCount, sizeof(biBitCount));
    f.write((char*)&biCompression, sizeof(biCompression));
    f.write((char*)&biSizeImage, sizeof(biSizeImage));
    f.write((char*)&biXPelsPerMeter, sizeof(biXPelsPerMeter));
    f.write((char*)&biYPelsPerMeter, sizeof(biYPelsPerMeter));
    f.write((char*)&biClrUsed, sizeof(biClrUsed));
    f.write((char*)&biClrImportant, sizeof(biClrImportant));

    f.write((char*)&pal, sizeof(pal));

    int r = 4 - (biWidth % 4);
    unsigned char temp = 0;

    for (int i = biHeight - 1; i >= 0; i--) {
        for (int j = 0; j < biWidth; j++)
            f.write((char*)&mat[i][j], sizeof(mat[i][j]));

        if (biWidth % 4 != 0)
            for (int k = 0; k < r; k++)
                f.write((char*)&temp, sizeof(temp));
    }
    f.flush();
    f.close();
}

unsigned char** immagine::conv(int dim, int** con) {
    /*
        Convoluzione con 0 - padding: dato un kernel di dimensione k, 
        l'immagine è orlata con (k-1)/2 righe e colonne su tutte le dimensioni con pixel valorizzati a zero.
        Esempio: immagine 6x5, kernel con k=3

                        0000000
        XXXXX           0XXXXX0
        XXXXX           0XXXXX0
        XXXXX   -->     0XXXXX0
        XXXXX           0XXXXX0
        XXXXX           0XXXXX0
        XXXXX           0XXXXX0
                        0000000
    */ 

    // Step 1: calcolare dimensioni nuova matrice (in questo caso, sono le stesse dell'immagine di partenza)
    int h = biHeight;
    int w = biWidth;

    // Spostamento rispetto al centro del kernel
    int offset = (int)(dim / 2);

    // Step 2: dichiarazione e allocazione di una matrice con le dimensioni calcolate
    //          al punto precedente in cui salvare il risultato

    unsigned char** B;
    B = new unsigned char* [h];
    for (int i = 0; i < h; i++) {
        B[i] = new unsigned char[w];
    }

    int sum; // Memorizza temporaneamente la somma pesata relativa al pixel
    for (int i = 0; i < biHeight; i++)
    {
        for (int j = 0; j < biWidth; j++)
        {
            // reset della variabile sum
            sum = 0;
            for (int k = 0; k < dim; k++)
            {
                for (int l = 0; l < dim; l++)
                {
                    // Calcolo degli indici di riga e di colonna di immagine e kernel per applicare la definizione di convoluzione
                    int ik = dim - (1 + k);
                    int jl = dim - (1 + l);
                    int ii = i + (k - offset);
                    int jj = j + (l - offset);

                    if (ii >= 0 && jj >= 0 && ii < biHeight && jj < biWidth)
                        sum += (con[ik][jl]) * (mat[ii][jj]);
                }
            }

            // Gestisco i casi in cui il valore dei pixel sfora il range consentito dalla codifica
            if (sum > 255) {
                B[i][j] = 255;
            }
            else if (sum < 0) {
                B[i][j] = 0;
            }
            else {
                B[i][j] = sum;
            }
        }

    }
    return B;

}

unsigned char** immagine::conv(int dim, float** con, bool padding)
{
    unsigned char** B;
    float sum; // Memorizza temporaneamente la somma pesata relativa al pixel
    int h = biHeight;
    int w = biWidth;

    // Spostamento rispetto al centro del kernel
    int offset = (int)(dim / 2);


    if (padding == true) {
        /*
        *   Convoluzione con 0 - padding: dato un kernel di dimensione k,
        *   l'immagine è orlata con (k-1)/2 righe e colonne su tutte le dimensioni con pixel valorizzati a zero.
        *   Esempio: immagine 6x5, kernel con k=3
        *
        *               0000000
        *   XXXXX       0XXXXX0
        *   XXXXX       0XXXXX0
        *   XXXXX   --> 0XXXXX0
        *   XXXXX       0XXXXX0
        *   XXXXX       0XXXXX0
        *   XXXXX       0XXXXX0
        *               0000000
        */

        // Step 1: calcolare dimensioni nuova matrice (in questo caso, sono le stesse dell'immagine di partenza)
        

        // Step 2: allocazione di una matrice con le dimensioni calcolate
        //          al punto precedente in cui salvare il risultato

        B = new unsigned char* [h];
        for (int i = 0; i < h; i++) {
            B[i] = new unsigned char[w];
        }

        for (int i = 0; i < biHeight; i++)
        {
            for (int j = 0; j < biWidth; j++)
            {
                // reset della variabile sum
                sum = 0.0;
                for (int k = 0; k < dim; k++)
                {
                    for (int l = 0; l < dim; l++)
                    {
                        // Calcolo degli indici di riga e di colonna di immagine e kernel per applicare la definizione di convoluzione
                        int ik = dim - (1 + k);
                        int jl = dim - (1 + l);
                        int ii = i + (k - offset);
                        int jj = j + (l - offset);

                        if (ii >= 0 && jj >= 0 && ii < biHeight && jj < biWidth)
                            sum += (con[ik][jl]) * (mat[ii][jj]);
                    }
                }

                // Gestisco i casi in cui il valore dei pixel sfora il range consentito dalla codifica
                if (sum > 255) {
                    B[i][j] = 255;
                }
                else if (sum < 0) {
                    B[i][j] = 0;
                }
                else {
                    B[i][j] = (unsigned char)sum;
                }
            }
        }
    }
    else {
        /*
        * Convoluzione senza padding -> riduzione delle dimensioni
        * dell'immagine di un fattore legato alla dimensione del kernel.
        * Se dim kernel = k, le dimensioni dell'immagine si riducono di (k-1) su righe e colonne
        */

        h = biHeight - dim + 1;
        w = biWidth - dim + 1;

        // Step 1: calcolare dimensioni nuova matrice
        B = new unsigned char* [h];
        for (int i = 0; i < h; i++) {
            B[i] = new unsigned char[w];
        }

        // Calcolo indici pixel di inizio e fine (devo considerare solo pixel ammissibili)
        for (int i = offset; i < biHeight-offset; i++)
        {
            for (int j = offset; j < biWidth-offset; j++)
            {
                // reset della variabile sum
                sum = 0.0;
                for (int k = 0; k < dim; k++)
                {
                    for (int l = 0; l < dim; l++)
                    {
                        // Calcolo degli indici di riga e di colonna di immagine e kernel per applicare la definizione di convoluzione
                        int ik = dim - (1 + k);
                        int jl = dim - (1 + l);
                        int ii = i + (k - offset);
                        int jj = j + (l - offset);

                        if (ii >= 0 && jj >= 0 && ii < biHeight && jj < biWidth)
                            sum += (con[ik][jl]) * (mat[ii][jj]);
                    }
                }

                // Gestisco i casi in cui il valore dei pixel sfora il range consentito dalla codifica
                if (sum > 255) {
                    B[i][j] = 255;
                }
                else if (sum < 0) {
                    B[i][j] = 0;
                }
                else {
                    B[i][j] = (unsigned char)sum;
                }
            }
        }

    }
    
    return B;
}

void immagine::creaScacchiera(int dim, int off) {

    for (int i = 0; i < biHeight; i++)
        delete mat[i];
    delete mat;

    biWidth = dim * off;
    biHeight = biWidth;

    mat = new unsigned char* [biHeight];
    for (int i = 0; i < biHeight; i++)
        mat[i] = new unsigned char[biWidth];

    for (int i = 0; i < biHeight; i++)
        for (int j = 0; j < biWidth; j++)
            if ((i / off + j / off) % 2 != 0)
                mat[i][j] = 255;
            else
                mat[i][j] = 0;

    char name[100];
    sprintf_s(name, "scacchiera_%d_%d.bmp", dim, off);

    salvaImmagine(name);

}

void compressioneRLE(ofstream& f, immagine im) {

    int h = im.getHeight(), w = im.getWidth();

    f.write((char*)&h, sizeof(h));
    f.write((char*)&w, sizeof(w));

    unsigned char value;
    int count, i;

    for (i = 0, value = im.getMatrix()[0][0], count = 0; i < h; i++)
        for (int j = 0; j < w; j++)
            if (((int)im.getMatrix()[i][j] == (int)value))
                count++;
            else {

                //cout << " " << count << " ";
                //cout << " " << value << " - " << endl;
                f.write((char*)&count, sizeof(count));
                f.write((char*)&value, sizeof(value));
                count = 1;
                value = im.getMatrix()[i][j];

            }
    //cout << " " << count << " ";
    //cout << " " << value << " - " << endl;
    f.write((char*)&count, sizeof(count));
    f.write((char*)&value, sizeof(value));

}

void decompressioneRLE(ifstream& f, immagine& im) {

    unsigned char value, ** mat;
    int count, h, w;

    f.read((char*)&h, sizeof(h));
    f.read((char*)&w, sizeof(w));

    mat = new unsigned char* [h];
    for (int i = 0; i < h; i++)
        mat[i] = new unsigned char[w];

    for (int i = 0, count = 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            if (count == 0) {
                f.read((char*)&count, sizeof(count));
                f.read((char*)&value, sizeof(value));
            }
            mat[i][j] = value;
            count--;
        }

    im.setMatrix(mat);

}
