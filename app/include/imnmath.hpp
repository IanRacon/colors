#ifndef IMNMATH_H
#define IMNMATH_H

#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>

using namespace std;

#define MAX_TIME_DATA_SAMPLES 500

enum gnuplot_splot_type
{
    GNUPLOT_PM3D = 1,
    GNUPLOT_CONTOUR = 2
};

struct imn_gnuplot_params
{
    int image_width;
    int image_height;
    std::string title;
    std::string xlabel;
    std::string ylabel;
    std::string cntrlevels;
    unsigned int stype;
    imn_gnuplot_params()
    {
        image_width = 1024;
        image_height = 1024;
        cntrlevels = "50";
        title = "Wygenerowany przez IMNMATH";
        xlabel = "os X";
        ylabel = "os Y";
        stype = GNUPLOT_PM3D;
    }
};

template <class T = double>
class imn
{
  public:
    static std::vector<double **> time_data2d;
    static std::vector<double **> time_data1d;
    static imn_gnuplot_params plot_params;

    /* Alokuje tablice dwuwymiarowa o wymiarze [nrowa][ncols]
 * Indeksy tablicy zaczynaja sie od 0 */
    static T **matrix(long nrh, long nch)
    {
        T **m;
        m = new T *[nrh];
        for (int i = 0; i < nrh; i++)
            m[i] = new T[nch];
        for (int i = 0; i < nrh; i++)
            for (int j = 0; j < nch; j++)
                m[i][j] = 0;
        return m;
    }
    /* Ustawia wartosci macierzy na wartosc value*/
    static void set_matrix(T **m, long nrh, long nch, T value)
    {
        for (int i = 0; i < nrh; i++)
            for (int j = 0; j < nch; j++)
                m[i][j] = value;
    }

    static void set_matrix_row(T **m, long row, long ncols, T value)
    {
        for (int i = 0; i < ncols; i++)
        {
            m[row][i] = value;
        }
    }

    /* Wypelnia wiersz row macierzy m wartoscia value */
    static void set_matrix_row(T **m, long row, long ncols, T *values)
    {
        for (int i = 0; i < ncols; i++)
        {
            m[row][i] = values[i];
        }
    }
    /* Wypelnia kolumne col macierzy m wartoscia value */
    static void set_matrix_col(T **m, long col, long nrows, T value)
    {
        for (int i = 0; i < nrows; i++)
        {
            m[i][col] = value;
        }
    }

    /* Wypelnia kolumne col macierzy m wartosciami values[i] */
    static void set_matrix_col(T **m, long col, long nrows, T *values)
    {
        for (int i = 0; i < nrows; i++)
        {
            m[i][col] = values[i];
        }
    }

    /* Kopiuje wartosci macierzy srcM do dstM*/
    static void copy_matrix(T **scrM, T **dstM, long nrh, long nch)
    {
        for (int i = 0; i < nrh; i++)
            for (int j = 0; j < nch; j++)
                dstM[i][j] = scrM[i][j];
    }

    /* Wykonuje mnozenie macierzy razy wektor */
    static void matrix_vector(T **mat, T *vec, T *vout, long nrh)
    {
        for (int i = 0; i < nrh; i++)
        {
            vout[i] = 0;
            for (int j = 0; j < nrh; j++)
                vout[i] += mat[i][j] * vec[j];
        }
    }

    /* Wykonuje operacje mnozenia macierzy D =  A * B, przy czym wszystkie macierze sa kwadratowe */
    static void matrix_matrix(T **matA, T **matB, T **matD, long nrh)
    {
        for (int i = 0; i < nrh; i++)
        {
            for (int j = 0; j < nrh; j++)
            {
                matD[i][j] = 0;
                for (int k = 0; k < nrh; k++)
                {
                    matD[i][j] += matA[i][k] * matB[k][j];
                }
            }
        }
    }

    /* Alokuje jednowymiarowy wektor o wymiarze [size] */
    static T *vector(long nrh)
    {
        T *v;
        v = new T[nrh];
        for (int i = 0; i < nrh; i++)
            v[i] = 0;
        return v;
    }

    /* Ustawia wartosci wektora na value */
    static void set_vector(T *v, long size, T value)
    {
        for (int i = 0; i < size; i++)
            v[i] = value;
    }

    /* Kopiuje wartosci wektora scrV do dstV */
    static void copy_vector(T *scrV, T *dstV, long size)
    {
        for (int i = 0; i < size; i++)
            dstV[i] = scrV[i];
    }

    /* Zwalnia pamiec zajeta przez wektor*/
    static void free_vector(T *v)
    {
        delete[] v;
    }
    /* Zwalnia macierz utworzona za pomoca procedury matrix(..)*/
    static void free_matrix(T **m, int nrows)
    {
        for (int i = 0; i < nrows; i++)
            delete[] m[i];
        delete[] m;
    }

    // -------------------------------------------------------------------------
    //                      Funkcje wejscia / wyjscia
    // -------------------------------------------------------------------------

    /**
 * Zapisuje macierz mat do pliku o nazwie filename. Liczba kolumna i wierszy moze
 * byc mniejsza niz rzeczywista rozmiar macierzy.
 * @param forceEnter - dodaje pusta linie co forceEnter liczbe wierszy.
 */
    static void write_matrix(const char *filename, T **mat, int nrows, int ncols, int forceEnter = 0)
    {
        FILE *f;
        int i, j;
        f = fopen(filename, "w");
        for (i = 0; i < nrows; i++)
        {
            for (j = 0; j < ncols; j++)
            {
                if (sizeof(T) == sizeof(double) || sizeof(T) == sizeof(float))
                    fprintf(f, "%20.10e\t", double(mat[i][j]));
                if (sizeof(T) == sizeof(int))
                    fprintf(f, "%16d\t", int(mat[i][j]));
            }
            fprintf(f, "\n");
            if (forceEnter != 0)
                if ((i + 1) % forceEnter == 0)
                    fprintf(f, "\n");
        }
        fclose(f);
    }
    /** 
* Zapisuje do pliku tablice 2d w formacie "x-y-wartosc". Podajemy dx, dy 
* do odpowiedniego przeskalowania osi.  nrows oraz ncols stanowia wymiary
* tablicy data
*/
    static void write_2d_system(const char *filename, T **data, int nx, int ny, double dx, double dy)
    {
        FILE *f;
        f = fopen(filename, "w");
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fprintf(f, "%e\t%e\t%e\n", dx * i, dy * j, double(data[i][j]));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
    static void write_2d_system(const char *filename, T *data, int nx, int ny, double dx, double dy)
    {
        FILE *f;
        f = fopen(filename, "w");
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fprintf(f, "%e\t%e\t%e\n", dx * i, dy * j, double(data[i * ny + j]));
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    /** 
* Generuje wykres 2D w formacie png na podstawie podanej tablicy danych.
* Wymagane jest posiadanie zainstalowanego programu gnuplot.
*/
    static void plot_2d_system(const char *outputname, T **data, int nx, int ny, double dx, double dy)
    {
        write_2d_system("tmp.txt", data, nx, ny, dx, dy);
        ofstream file("tmp.plt");
        file << "set terminal png size " << plot_params.image_width << "," << plot_params.image_height << " enhanced font 'Helvetica,12'\n";
        file << "set o '" << outputname << "'\n";
        file << "set xl '" << plot_params.xlabel << "'\n";
        file << "set yl '" << plot_params.ylabel << "'\n";
        file << "set view map\n";
        file << "set size ratio -1\n";
        file << "set title '" << plot_params.title << "'\n";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << "set contours\n";
            file << "set contour base\n";
            file << "set cntrparam levels " << plot_params.cntrlevels << " \n";
            file << "unset clabel\n";
            if (!(plot_params.stype & GNUPLOT_PM3D))
                file << "unset surface\n";
        }
        file << "splot 'tmp.txt' u 1:2:3 w pm3d";

        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << " lt -1";
        }
        file << " t ''\n";
        file.close();
        int res = system("gnuplot tmp.plt");
        res = system("rm tmp.*");
    }
    static void plot_2d_system(const char *outputname, T *data, int nx, int ny, double dx, double dy)
    {
        write_2d_system("tmp.txt", data, nx, ny, dx, dy);
        ofstream file("tmp.plt");
        file << "set terminal png size " << plot_params.image_width << "," << plot_params.image_height << " enhanced font 'Helvetica,12'\n";
        file << "set o '" << outputname << "'\n";
        file << "set xl '" << plot_params.xlabel << "'\n";
        file << "set yl '" << plot_params.ylabel << "'\n";
        file << "set view map\n";
        file << "set size ratio -1\n";
        file << "set title '" << plot_params.title << "'\n";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << "set contours\n";
            file << "set contour base\n";
            file << "set cntrparam levels " << plot_params.cntrlevels << " \n";
            file << "unset clabel\n";
            if (!(plot_params.stype & GNUPLOT_PM3D))
                file << "unset surface\n";
        }
        file << "splot 'tmp.txt' u 1:2:3 w pm3d";

        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << " lt -1";
        }
        file << " t ''\n";
        file.close();
        int res = system("gnuplot tmp.plt");
        res = system("rm tmp.*");
    }
    static void plot_2d_system(const char *outputname, const std::vector<double> &data, int nx, int ny, double dx, double dy)
    {
        write_2d_system("tmp.txt", data, nx, ny, dx, dy);
        ofstream file("tmp.plt");
        file << "set terminal png size " << plot_params.image_width << "," << plot_params.image_height << " enhanced font 'Helvetica,12'\n";
        file << "set o '" << outputname << "'\n";
        file << "set xl '" << plot_params.xlabel << "'\n";
        file << "set yl '" << plot_params.ylabel << "'\n";
        file << "set view map\n";
        file << "set size ratio -1\n";
        file << "set title '" << plot_params.title << "'\n";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << "set contours\n";
            file << "set contour base\n";
            file << "set cntrparam levels " << plot_params.cntrlevels << " \n";
            file << "unset clabel\n";
            if (!(plot_params.stype & GNUPLOT_PM3D))
                file << "unset surface\n";
        }
        file << "splot 'tmp.txt' u 1:2:3 w pm3d";

        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << " lt -1";
        }
        file << " t ''\n";
        file.close();
        int res = system("gnuplot tmp.plt");
        res = system("rm tmp.*");
    }
    /**
* Generuje wykres 2D w formacie png na podstawie podanej tablicy danych.
* Wymagane jest posiadanie zainstalowanego programu gnuplot.
*/
    static void plot_2d_system(const char *outputname, T **data1, T **data2, int nx, int ny, double dx, double dy)
    {
        write_2d_system("tmp1.txt", data1, nx, ny, dx, dy);
        write_2d_system("tmp2.txt", data2, nx, ny, dx, dy);
        ofstream file("tmp.plt");
        file << "set terminal png size " << plot_params.image_width << "," << plot_params.image_height << " enhanced font 'Helvetica,12'\n";
        file << "set o '" << outputname << "'\n";
        file << "set xl '" << plot_params.xlabel << "'\n";
        file << "set yl '" << plot_params.ylabel << "'\n";
        file << "set view map\n";
        file << "set size ratio -1\n";
        file << "set title '" << plot_params.title << "'\n";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << "set contours\n";
            file << "set contour base\n";
            file << "set cntrparam levels 30\n";
            file << "unset clabel\n";
            if (!(plot_params.stype & GNUPLOT_PM3D))
                file << "unset surface\n";
        }
        file << "splot 'tmp1.txt' u 1:2:3 w pm3d t ''";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << " lt -1";
        }
        file << ", 'tmp2.txt' u 1:2:3 w pm3d t ''";
        if (plot_params.stype & GNUPLOT_CONTOUR)
        {
            file << " lt -1";
        }
        file << "\n";
        file.close();
        int res = system("gnuplot tmp.plt");
        res = system("rm tmp.* tmp1.txt tmp2.txt");
    }

    /**
 * Wypisuje macierz do standardowego wyjscia.
 * @param title - naglowek wypisywania
 * @param format - format zapisu w funkcji fprintf, np. dla double moze byc %16.6e.
 * @param mat   - wskaznik do macierzy
 * @param nrows - liczba wierszy
 * @param ncols - liczba kolumn
 */
    static void print_matrix(const char *title, const char *format, T **mat, int nrows, int ncols)
    {
        printf("%s\n", title);
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                printf(format, double(mat[i][j]));
            }
            printf("\n");
        }
    }

    static void push_data2D(double **data, int nrows, int ncols)
    {

        if (time_data2d.size() > MAX_TIME_DATA_SAMPLES)
        {
            cout << "Za duzo probek... Maksymalna liczba bufora pamieci wynosi:" << MAX_TIME_DATA_SAMPLES << endl;
            return;
        }

        double **in_data = matrix(nrows, ncols);
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                in_data[i][j] = data[i][j];
            }
        }

        time_data2d.push_back(in_data);
        if (time_data2d.size() == 1)
        {
            printf("Dodawanie danych:\n");
        }
        else
        {
            printf("%5d \r", int(time_data2d.size()));
            fflush(stdout);
        }
    }
    static void push_data2D(double *data, int nrows, int ncols)
    {

        if (time_data2d.size() > MAX_TIME_DATA_SAMPLES)
        {
            cout << "Za duzo probek... Maksymalna liczba bufora pamieci wynosi:" << MAX_TIME_DATA_SAMPLES << endl;
            return;
        }

        double **in_data = matrix(nrows, ncols);
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                in_data[i][j] = data[i * ncols + j];
            }
        }

        time_data2d.push_back(in_data);
        if (time_data2d.size() == 1)
        {
            printf("Dodawanie danych:\n");
        }
        else
        {
            printf("%5d \r", int(time_data2d.size()));
            fflush(stdout);
        }
    }

    static void push_vector_data2D(std::vector<double> data, int nrows, int ncols)
    {
        if (time_data2d.size() > MAX_TIME_DATA_SAMPLES)
        {
            cout << "Za duzo probek... Maksymalna liczba bufora pamieci wynosi:" << MAX_TIME_DATA_SAMPLES << endl;
            return;
        }

        double **in_data = matrix(nrows, ncols);
        for (int i = 0; i < nrows; i++)
            for (int j = 0; j < ncols; j++)
                in_data[i][j] = data[i * ncols + j];

        time_data2d.push_back(in_data);

        if (time_data2d.size() == 1)
            printf("Dodawanie danych:\n");
        else
        {
            printf("%5d \r", int(time_data2d.size()));
            fflush(stdout);
        }
    }

    /**
* Zwalniamy pamiec historii.
**/
    static void free_data2D(int nrows, int ncols)
    {
        for (int i = 0; i < int(time_data2d.size()); i++)
        {
            free_matrix(time_data2d[i], nrows);
        }
        time_data2d.clear();
    }

    /**
* Zapisuje sekwencje danych w tablicy data2d do pliku.
**/
    static void write_data2D(const char *filename, int nrows, int ncols, double dx, double dy)
    {
        FILE *f;
        f = fopen(filename, "w");
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                fprintf(f, "%e\t%e\t", dx * i, dy * j);
                for (int t = 0; t < int(time_data2d.size()); t++)
                {
                    fprintf(f, "%e\t", time_data2d[t][i][j]);
                }
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    // -------------------------------------------------------------------------- //
    //
    // Funkcje pomocnicze do wczytywania danych z pliku napisane przez p. Wach.
    //
    // -------------------------------------------------------------------------- //
    /**
 *  Funkcja do wczytywania pola predkosci z laboratorium "Adwekcja w 2D" z pliku o formacie
 *  x y predkoscX predkoscY
 *  Tablice u i v nalezy zaalokowac przed wywolaniem tej funkcji 
 *  nx,ny okresjala wymiary wczytywanego ukladu
*/
    static void uv_load(double **u, double **v, int nx, int ny, const char *filename)
    {
        FILE *f;
        int i, j;
        double x, y, u_temp, v_temp;
        char line[200];
        if (!(f = fopen(filename, "rt")))
        {
            printf("Brak pliku\n");
            exit(0);
        }

        i = 0;
        j = 0;
        while (fgets(line, 200, f) != NULL)
        {
            sscanf(line, "%lf %lf %lf %lf\n", &x, &y, &u_temp, &v_temp);

            const char accept[] = " \t\r\n";
            bool test = (strspn(line, accept) == strlen(line));

            if (!test)
            {
                u[i][j] = u_temp;
                v[i][j] = v_temp;
                j++;
                if (j == (ny))
                {
                    j = 0;
                    i++;
                }
            }
        }
        fclose(f);
    }
};

template <class T>
std::vector<double **> imn<T>::time_data2d;
template <class T>
std::vector<double **> imn<T>::time_data1d;
template <class T>
imn_gnuplot_params imn<T>::plot_params;

typedef imn<float> imnf;
typedef imn<double> imnd;
typedef imn<int> imni;
#endif /* IMNMATH_H */
