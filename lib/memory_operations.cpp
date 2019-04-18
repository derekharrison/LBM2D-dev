/*
 * memory_operations.cpp
 *
 *  Created on: Apr 2, 2019
 *      Author: d-w-h
 */


#include <cstdlib>

bool *allocate_1D_array_bool(const int Nx)
{
    bool *f = new bool[Nx];

    return f;

}


void deallocate_1D_array_bool(bool* f)
{
    delete[] f;

}

double *allocate_1D_array_doubles(const int Nx)
{
    double *f = new double[Nx];

    return f;

}


void deallocate_1D_array_doubles(double* f)
{
    delete[] f;

}

int **allocate_2D_array_int(const int Nx, const int Ny)
{
    int **f = new int*[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new int[Ny];

    return f;

}


void deallocate_2D_array_int(int** f,const int Nx)
{
    for(int i = 0; i < Nx; ++i)
            delete[] f[i];
    delete[] f;

}

double **allocate_2D_array_doubles(const int Nx, const int Ny)
{
    double **f = new double*[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new double[Ny];

    return f;

}


void deallocate_2D_array_doubles(double** f,const int Nx)
{
    for(int i = 0; i < Nx; ++i)
            delete[] f[i];
    delete[] f;

}


double ***allocate_3D_array_doubles(const int Nx,const int Ny,const int Nz)
{
    double ***f = new double**[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new double*[Ny];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            f[i][j] = new double[Nz];

    return f;

}


void deallocate_3D_array_doubles(double*** f,const int Nx,const int Ny)
{
    for(int i = 0; i < Nx; ++i)
    {
        for(int j = 0; j < Ny; ++j)
            delete[] f[i][j];
        delete[] f[i];
    }
    delete[] f;

}

bool **allocate_2D_array_bool(const int Nx,const int Ny)
{
    bool **f = new bool*[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new bool[Ny];

    return f;

}

void deallocate_2D_array_bool(bool*** f,const int Nx)
{
    for(int i = 0; i < Nx; ++i)
            delete[] f[i];
    delete[] f;

}

bool ***allocate_3D_array_bool(const int Nx,const int Ny,const int Nz)
{
    bool ***f = new bool**[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new bool*[Ny];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            f[i][j] = new bool[Nz];

    return f;

}

void deallocate_3D_array_bool(bool*** f,const int Nx,const int Ny)
{
    for(int i = 0; i < Nx; ++i)
    {
        for(int j = 0; j < Ny; ++j)
            delete[] f[i][j];
        delete[] f[i];
    }
    delete[] f;

}

double ****allocate_discrete_velocity_3D(const int Nx,const int Ny,const int Nz,const int Nd)
{
    double ****f = new double***[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new double**[Ny];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            f[i][j] = new double*[Nz];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            for(int k = 0; k < Nz; ++k)
                f[i][j][k] = new double[Nd];

    return f;

}


void deallocate_discrete_velocity_3D(double ****f, const int Nx,const int Ny,const int Nz)
{
    for(int i = 0; i < Nx; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
                delete[] f[i][j][k];
            delete[] f[i][j];
        }
        delete[] f[i];
    }

    delete[] f;

}

double ***allocate_discrete_velocity_2D(const int Nx,const int Ny, const int Nd)
{
    double ***f = new double**[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new double*[Ny];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            f[i][j] = new double[Nd];

    return f;

}


void deallocate_discrete_velocity_2D(double ****f, const int Nx,const int Ny,const int Nz)
{
    for(int i = 0; i < Nx; ++i)
    {
        for(int j = 0; j < Ny; ++j)
            delete[] f[i][j];
        delete[] f[i];
    }
    delete[] f;

}

int ****allocate_4D_array_int(const int Nx,const int Ny,const int Nz,const int Nd)
{
    int ****f = new int***[Nx];

    for(int i = 0; i < Nx; ++i)
        f[i] = new int**[Ny];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            f[i][j] = new int*[Nz];

    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            for(int k = 0; k < Nz; ++k)
                f[i][j][k] = new int[Nd];

    return f;

}


void deallocate_4D_array_int(int ****f, const int Nx,const int Ny,const int Nz)
{
    for(int i = 0; i < Nx; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
                delete[] f[i][j][k];
            delete[] f[i][j];
        }
        delete[] f[i];
    }

    delete[] f;

}
