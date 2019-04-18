/*
 * main.cpp
 *
 *  Created on: Apr 2, 2019
 *      Author: d-w-h
 */

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <vector>
#include <string>

#include "../inc/memory_operations.hpp"

typedef struct node_coordinates {
    int i;
    int j;
} nd_crdnts;

int main(int argc, char* argv[])
{
    /*Parameter declaration*/
    int Nx, Ny, Nd;
    int max_timesteps, timestep;
    double dt, tau, dt_tau, dx;
    double x_x, x_y, box_nx, box_ny, radius;
    double omega_z, ux, uy, I, m, px, py;
    double Fx, Fy, Fm, Fr, g;
    int size_box;
    std::map<int, int> index_map;


    Nx = 100;
    Ny = 100;
    Nd = 9;

    max_timesteps = 500;
    dt = 1.0;
    dx = 1.0;
    tau = 0.6;

    radius = 6.0;
    m = 1.0;

    x_x = 50;
    x_y = 90;
    ux = 0.0;
    uy = 0.0;
    px = 0.0;
    py = 0.0;
    omega_z = 0.0;
    g = -9.81;
    size_box = radius/dx + 3;
    dt_tau = dt/tau;
    I = 0.5*m*radius*radius;

    /*Memory allocation*/
    double*** fi     = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double*** fieq   = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double*** fiprop = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double**     rho = allocate_2D_array_doubles(Nx, Ny);
    double**      Ux = allocate_2D_array_doubles(Nx, Ny);
    double**      Uy = allocate_2D_array_doubles(Nx, Ny);
    double** vort_field = allocate_2D_array_doubles(Nx, Ny);
    double**       X = allocate_2D_array_doubles(Nx, Ny);
    double**       Y = allocate_2D_array_doubles(Nx, Ny);
    double*       wi = allocate_1D_array_doubles(Nx);
    double*      cix = allocate_1D_array_doubles(Nx);
    double*      ciy = allocate_1D_array_doubles(Nx);
    bool** solid_flags = allocate_2D_array_bool(Nx, Ny);
    int** adj_flags = allocate_2D_array_int(Nx, Ny);
    nd_crdnts* adj_nodes = (nd_crdnts*)malloc(sizeof(nd_crdnts) * 0.5*Nx*Ny);

    /*Initialization*/
    /*Initialize index map*/
    index_map.insert(std::pair <int, int> (0, 0));
    index_map.insert(std::pair <int, int> (1, 3));
    index_map.insert(std::pair <int, int> (2, 4));
    index_map.insert(std::pair <int, int> (3, 1));
    index_map.insert(std::pair <int, int> (4, 2));
    index_map.insert(std::pair <int, int> (5, 7));
    index_map.insert(std::pair <int, int> (6, 8));
    index_map.insert(std::pair <int, int> (7, 5));
    index_map.insert(std::pair <int, int> (8, 6));

    /*Initializing cij and wi*/
    for(int i = 0; i < 27; ++i)
    {
        cix[0] = 0; ciy[0] = 0; wi[0] = 4.0/9.0;
        cix[1] = 1; ciy[1] = 0; wi[1] = 1.0/9.0;
        cix[2] = 0; ciy[2] = 1; wi[2] = 1.0/9.0;
        cix[3] = -1;ciy[3] = 0; wi[3] = 1.0/9.0;
        cix[4] = 0; ciy[4] = -1;wi[4] = 1.0/9.0;
        cix[5] = 1; ciy[5] = 1; wi[5] = 1.0/36.0;
        cix[6] = -1;ciy[6] = 1; wi[6] = 1.0/36.0;
        cix[7] = -1;ciy[7] = -1;wi[7] = 1.0/36.0;
        cix[8] = 1; ciy[8] = -1;wi[8] = 1.0/36.0;
    }

    /*Initializing density and velocity*/
    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
        {
            Ux[i][j] = 0.0;
            Uy[i][j] = 0.0;
            vort_field[i][j] = 0.0;
            X[i][j] = i*dx + 0.5*dx;
            Y[i][j] = j*dx + 0.5*dx;
            rho[i][j] = 1000.0;
            solid_flags[i][j] = false;
            adj_flags[i][j] = 0;
        }

    /*Initializing discrete velocity distributions*/
    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            for(int d = 0; d < Nd; ++d)
            {
                fi[i][j][d] = wi[d]*rho[i][j];
                fieq[i][j][d] = wi[d]*rho[i][j];
                fiprop[i][j][d] = wi[d]*rho[i][j];
            }


    /*LBM algorithm*/
    timestep = 0;
    do
    {
        /*Calculate density and velocity*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j)
            {
                rho[i][j] = 0.0;
                Ux[i][j] = 0.0;
                Uy[i][j] = 0.0;
                for(int d = 0; d < Nd; ++d)
                {
                    rho[i][j] += fiprop[i][j][d];
                    Ux[i][j] += cix[d] * fiprop[i][j][d];
                    Uy[i][j] += ciy[d] * fiprop[i][j][d];
                }
                Ux[i][j] = Ux[i][j]/rho[i][j];
                Uy[i][j] = Uy[i][j]/rho[i][j];
            }

        /*Calculate vorticity*/
        for(int i = 1; i < Nx - 1; ++i)
            for(int j = 1; j < Ny - 1; ++j)
            {
                if(solid_flags[i][j] == 0)
                {
                    vort_field[i][j] = 0.5 * (Uy[i+1][j] - Uy[i-1][j] - Ux[i][j+1] + Ux[i][j-1]) * 1.0 / dx;
                }
            }

        /*Calculate max velocity and vorticity*/
        double max_v_y = 0.0;
        double min_v_y = 1.0;
        double max_velocity = 0.0;
        double max_vorticity = 0.0;
        double dummy_val_vel = 0.0;
        double dummy_val_vort = 0.0;
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j)
            {
                dummy_val_vel  = Ux[i][j] * Ux[i][j] + Uy[i][j] * Uy[i][j];
                dummy_val_vort = vort_field[i][j] * vort_field[i][j];

                if(Uy[i][j] > max_v_y)
                {
                    max_v_y = Uy[i][j];
                }
                if(Uy[i][j] < min_v_y)
                {
                    min_v_y = Uy[i][j];
                }


                if(dummy_val_vel > max_velocity)
                {
                    max_velocity = dummy_val_vel;
                }

                if(dummy_val_vort > max_vorticity)
                {
                    max_vorticity = dummy_val_vort;
                }
            }
        max_velocity  = sqrt(max_velocity);
        max_vorticity = sqrt(max_vorticity);

        /*Calculate equilibrium distribution*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j)
            {
                for(int d = 0; d < Nd; ++d)
                {
                    double uci = cix[d]*Ux[i][j] + ciy[d]*Uy[i][j];
                    double uu = Ux[i][j]*Ux[i][j] + Uy[i][j]*Uy[i][j];
                    fieq[i][j][d] = wi[d]*rho[i][j]*(1.0 + 3.0*uci + 4.5*uci*uci - 1.5*uu);
                }
            }

        /*Collision*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j)
            {
                for(int d = 0; d < Nd; ++d)
                {
                    fi[i][j][d] = fiprop[i][j][d]*(1 - dt_tau) + fieq[i][j][d]*dt_tau;
                }
            }

        /*Streaming internal nodes*/
        for(int i = 1; i < Nx - 1; ++i)
            for(int j = 1; j < Ny - 1; ++j)
            {
                bool is_fluid = solid_flags[i][j] == false;
                if(is_fluid)
                {
                    fiprop[i][j][0] = fi[i][j][0];
                    fiprop[i][j][1] = fi[i-1][j][1];
                    fiprop[i][j][2] = fi[i][j-1][2];
                    fiprop[i][j][3] = fi[i+1][j][3];
                    fiprop[i][j][4] = fi[i][j+1][4];
                    fiprop[i][j][5] = fi[i-1][j-1][5];
                    fiprop[i][j][6] = fi[i+1][j-1][6];
                    fiprop[i][j][7] = fi[i+1][j+1][7];
                    fiprop[i][j][8] = fi[i-1][j+1][8];
                }
            }

        /*Streaming south side*/
        for(int i = 0; i < Nx; ++i)
            for(int d = 0; d < Nd; ++d)
            {
                if(ciy[d] == 1)
                {
                    fiprop[i][0][d] = fi[i][0][index_map[d]];
                }
            }

        /*Streaming north side*/
        for(int i = 0; i < Nx; ++i)
            for(int d = 0; d < Nd; ++d)
            {
                if(ciy[d] == -1)
                {
                    fiprop[i][Ny-1][d] = fi[i][Ny-1][index_map[d]];
                }
            }

        /*Streaming west side*/
        for(int j = 0; j < Ny; ++j)
            for(int d = 0; d < Nd; ++d)
            {
                if(cix[d] == 1)
                {
                    fiprop[0][j][d] = fi[0][j][index_map[d]];
                }
            }

        /*Streaming east side*/
        for(int j = 0; j < Ny; ++j)
            for(int d = 0; d < Nd; ++d)
            {
                if(cix[d] == -1)
                {
                    fiprop[Nx-1][j][d] = fi[Nx-1][j][index_map[d]];
                }
            }

        /*Detecting solid nodes*/
        box_nx = x_x/dx;
        box_ny = x_y/dx;
        int num_solid_nodes = 0;
        for(int i = box_nx - size_box; i <= (box_nx + size_box) && (i >= 0) && (i <= Nx - 1); ++i)
            for(int j = box_ny - size_box; j <= (box_ny + size_box) && (j >= 0) && (j <= Ny - 1); ++j)
            {
                double del_x, del_y, del_r;
                del_x = i*dx + 0.5*dx - x_x;
                del_y = j*dx + 0.5*dx - x_y;
                del_r = sqrt(del_x*del_x + del_y*del_y);
                if(del_r <= radius)
                {
                    solid_flags[i][j] = true;
                    num_solid_nodes++;
                }
                else
                {
                    solid_flags[i][j] = false;
                }
            }

        /*Detecting nodes adjacent to solid regions*/
        int adj_node_counter = 0;
        for(int i = box_nx - size_box; (i <= box_nx + size_box) && (i >= 0) && (i <= Nx - 1); ++i)
            for(int j = box_ny - size_box; (j <= box_ny + size_box) && (j >= 0) && (j <= Ny - 1); ++j)
            {
                bool is_adjacent = false;
                for(int ii = i - 1; (ii <= i + 1) && (ii >= 0) && (ii <= Nx - 1); ++ii)
                    for(int jj = j - 1; (jj <= j + 1) && (jj >= 0) && (jj <= Ny - 1); ++jj)
                    {
                        bool hit_solid = (solid_flags[ii][jj] == true) && (solid_flags[i][j] == false);
                        if(hit_solid)
                        {
                            is_adjacent = true;
                        }
                    }
                if(is_adjacent)
                {
                    adj_nodes[adj_node_counter].i = i;
                    adj_nodes[adj_node_counter].j = j;
                    adj_node_counter++;

                    adj_flags[i][j] = 2; //for testing
                }
                else
                {
                    adj_flags[i][j] = 0; //for testing
                }
            }

        /*Streaming solid boundary nodes*/
        int n_nodes = adj_node_counter;
        double del_omega_z = 0.0;
        double del_px = 0.0;
        double del_py = 0.0;
        for(int n = 0; n < n_nodes; ++n)
        {
            int i, j, is, js;
            bool hit_solid = false;

            i = adj_nodes[n].i;
            j = adj_nodes[n].j;

            //east node solid, 1
            is = i + 1;
            js = j;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] + 0.5*dx - x_x;
                ry = Y[i][j] - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][3] = fi[i][j][1] - 2*wi[1]*rho[i][j]*3.0*(cix[1]*uwx + ciy[1]*uwy);

                delpx = cix[1]*fiprop[i][j][1] - cix[3]*fi[i][j][3];
                delpy = ciy[1]*fiprop[i][j][1] - ciy[3]*fi[i][j][3];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //west node solid, 2
            is = i - 1;
            js = j;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] - 0.5*dx - x_x;
                ry = Y[i][j] - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][1] = fi[i][j][3] - 2*wi[3]*rho[i][j]*3.0*(cix[3]*uwx + ciy[3]*uwy);

                delpx = cix[3]*fiprop[i][j][3] - cix[1]*fi[i][j][1];
                delpy = ciy[3]*fiprop[i][j][3] - ciy[1]*fi[i][j][1];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //north node solid, 3
            is = i;
            js = j + 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] - x_x;
                ry = Y[i][j] + 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][4] = fi[i][j][2] - 2*wi[2]*rho[i][j]*3.0*(cix[2]*uwx + ciy[2]*uwy);

                delpx = cix[2]*fiprop[i][j][2] - cix[4]*fi[i][j][4];
                delpy = ciy[2]*fiprop[i][j][2] - ciy[4]*fi[i][j][4];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //south node solid, 4
            is = i;
            js = j - 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] - x_x;
                ry = Y[i][j] - 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][2] = fi[i][j][4] - 2*wi[4]*rho[i][j]*3.0*(cix[4]*uwx + ciy[4]*uwy);

                delpx = cix[4]*fiprop[i][j][4] - cix[2]*fi[i][j][2];
                delpy = ciy[4]*fiprop[i][j][4] - ciy[2]*fi[i][j][2];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //north-east node solid, 5
            is = i + 1;
            js = j + 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] + 0.5*dx - x_x;
                ry = Y[i][j] + 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][7] = fi[i][j][5] - 2*wi[5]*rho[i][j]*3.0*(cix[5]*uwx + ciy[5]*uwy);

                delpx = cix[5]*fiprop[i][j][5] - cix[7]*fi[i][j][7];
                delpy = ciy[5]*fiprop[i][j][5] - ciy[7]*fi[i][j][7];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //north-west node solid, 6
            is = i - 1;
            js = j + 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] - 0.5*dx - x_x;
                ry = Y[i][j] + 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][8] = fi[i][j][6] - 2*wi[6]*rho[i][j]*3.0*(cix[6]*uwx + ciy[6]*uwy);

                delpx = cix[6]*fiprop[i][j][6] - cix[8]*fi[i][j][8];
                delpy = ciy[6]*fiprop[i][j][6] - ciy[8]*fi[i][j][8];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //south-west node solid, 7
            is = i - 1;
            js = j - 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] - 0.5*dx - x_x;
                ry = Y[i][j] - 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][5] = fi[i][j][7] - 2*wi[7]*rho[i][j]*3.0*(cix[7]*uwx + ciy[7]*uwy);

                delpx = cix[7]*fiprop[i][j][7] - cix[5]*fi[i][j][5];
                delpy = ciy[7]*fiprop[i][j][7] - ciy[5]*fi[i][j][5];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }

            //south-east node solid, 8
            is = i + 1;
            js = j - 1;
            hit_solid = solid_flags[is][js];
            if(hit_solid)
            {
                double uwx, uwy, rx, ry;
                double delpx, delpy;

                rx = X[i][j] + 0.5*dx - x_x;
                ry = Y[i][j] - 0.5*dx - x_y;

                uwx = ux - omega_z*ry;
                uwy = uy + omega_z*rx;

                fiprop[i][j][6] = fi[i][j][8] - 2*wi[8]*rho[i][j]*3.0*(cix[8]*uwx + ciy[8]*uwy);

                delpx = cix[8]*fiprop[i][j][8] - cix[6]*fi[i][j][6];
                delpy = ciy[8]*fiprop[i][j][8] - ciy[6]*fi[i][j][6];

                del_px += delpx;
                del_py += delpy;
                del_omega_z += (1.0/I) * (delpx*rx - delpy*ry);
            }
        }

        /*Write time dependent data*/
        FILE* params = NULL;
        params = fopen("parameters.txt", "w");
        if(params != NULL)
        {
            fprintf(params, "%i\n", Nx);
            fprintf(params, "%i\n", Ny);
            fprintf(params, "%i\n", timestep);
            fprintf(params,"%f\n", max_velocity * 1.05);
            fprintf(params,"%f\n", max_vorticity * 1.05);
            fprintf(params,"%f\n", min_v_y * 1.05);
            fprintf(params,"%f\n", max_v_y * 1.05);
        }
        else
        {
            printf("could not open file params\n");
            exit(2);
        }
        fclose(params);

        FILE* object_location = NULL;
        object_location = fopen("object_location.txt", "a");
        if(object_location != NULL)
        {
            fprintf(object_location,"%f    ", x_x);
            fprintf(object_location,"%f\n", x_y);
        }
        else
        {
            printf("could not open file object_location\n");
            exit(2);
        }
        fclose(object_location);

        FILE* data = NULL;
        char file_name_complete[50];
        char* file_name_suffix = (char*) ".txt";
        char* file_name = (char*) "./simdata/data";

        sprintf(file_name_complete, "%s_%d%s", file_name, timestep, file_name_suffix);
        data = fopen(file_name_complete, "w");
        if(data != NULL)
        {
            for(int j = Ny - 1;j >= 0 ; --j)
                for(int i = 0;i < Nx; ++i)
                {
                    fprintf(data,"%f    ", i*dx + 0.5*dx);
                    fprintf(data,"%f    ", j*dx + 0.5*dx);
                    fprintf(data,"%f    ", Ux[i][j]);
                    fprintf(data,"%f    ", Uy[i][j]);
                    fprintf(data,"%f\n", vort_field[i][j]);
                }
        }
        else
        {
            printf("could not open file params\n");
            exit(2);
        }
        fclose(data);

        px = px + del_px;
        py = py + del_py;

        Fx = del_px/dt;
        Fy = del_py/dt;
        Fm = m*g;

        Fr = Fm + Fy;

        x_x = x_x + 0.5 * Fx/m * dt*dt;
        x_y = x_y + 0.5 * Fr/m * dt*dt;

        omega_z = omega_z + del_omega_z;

//        ux = 0.0;//px/m;
//        uy = -0.25;//py/m;
//        x_y = x_y + uy*dt;


        timestep++;

//        printf("px: %f\tpy: %f\n", px, py);
//        printf("del_px: %f\tdel_py: %f\n", del_px, del_py);
        printf("x_x: %f\tx_y: %f\n", x_x, x_y);
        printf("ux: %f\tuy: %f\n", ux, uy);
        printf("del_px: %f\tdel_py: %f\n", del_px, del_py);
        printf("num solid nodes: %i\n", num_solid_nodes);
        printf("num adj nodes: %i\n", adj_node_counter);
//        printf("omega_z: %f\tdel_omega_z: %f\n", omega_z, del_omega_z);

    } while(timestep < max_timesteps);
    /*End LBM algorithm*/

    double sum = 0.0;
    for(int d = 0; d < Nd; ++d)
    {
        sum += fiprop[2][2][d];
    }
    printf("sum: %f\n", sum);
    /*Processing results*/
    std::cout << "done bro!" << std::endl;

    return 0;
}
