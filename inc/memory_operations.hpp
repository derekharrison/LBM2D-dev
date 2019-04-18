/*
 * memory_operations.hpp
 *
 *  Created on: Apr 2, 2019
 *      Author: d-w-h
 */

#ifndef MEMORY_OPERATIONS_HPP_
#define MEMORY_OPERATIONS_HPP_

double *allocate_1D_array_doubles(const int Nx);
void deallocate_1D_array_doubles(double* f);
double **allocate_2D_array_doubles(const int Nx,const int Ny);
void deallocate_2D_array_doubles(double** f,const int Nx);
double ***allocate_3D_array_doubles(const int Nx,const int Ny,const int Nz);
void deallocate_3D_array_doubles(double*** f,const int Nx,const int Ny);
double ****allocate_discrete_velocity_3D(const int Nx,const int Ny,const int Nz,const int Nd);
void deallocate_discrete_velocity_3D(double ****f, const int Nx,const int Ny,const int Nz);
double ***allocate_discrete_velocity_2D(const int Nx,const int Ny, const int Nd);
void deallocate_discrete_velocity_2D(double ****f, const int Nx,const int Ny,const int Nz);
int ****allocate_4D_array_int(const int Nx,const int Ny,const int Nz,const int Nd);
void deallocate_4D_array_int(int ****f, const int Nx,const int Ny,const int Nz);
int **allocate_2D_array_int(const int Nx, const int Ny);
void deallocate_2D_array_int(int** f,const int Nx);

bool ***allocate_3D_array_bool(const int Nx,const int Ny,const int Nz);
void deallocate_3D_array_bool(bool*** f,const int Nx,const int Ny);
bool *allocate_1D_array_bool(const int Nx);
void deallocate_1D_array_bool(bool* f);
bool **allocate_2D_array_bool(const int Nx,const int Ny);
void deallocate_2D_array_bool(bool*** f,const int Nx);

#endif /* MEMORY_OPERATIONS_HPP_ */
