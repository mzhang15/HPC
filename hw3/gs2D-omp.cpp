/******************************************************************************
* FILE: gs2D-omp.cpp
* DESCRIPTION:
*   Implementation of Gauss-Seidel solver 
*   Report timing and error, which is compared to Jacobi's implementation 
* AUTHOR: Mengyang Zhang
* LAST REVISED: 03/28/2021
* COMPILE:
*   g++ -Wall -pedantic -std=c++11 -fopenmp gs2D-omp.cpp && ./a.out N
*   (N is the size of the matrix)
******************************************************************************/
 
#include <iostream>
#include <vector>
#include <math.h>
#include "jacobi2D-omp.cpp"

using namespace std;

class GS2D {
public:
    GS2D(int N) : N(N) {
        h = (double) 1 / (N + 1);
        // printf("h = %7.3f\n", h);
    }

    vector<vector<double>> solution(int max_iters) {
        vector<vector<double>> u(N, vector<double>(N, 0.0)), u_temp;
        vector<vector<double>> f(N, vector<double>(N, 1.0));

        double tol = 1e-5;

        for (int i = 0; i < max_iters; ++i) {
            u_temp = u;
            updateRed(u, f); 
            updateBlack(u, f);

            if (square_loss(u, u_temp) < tol) break;
        }

        return u;
    }
private:
    double square_loss(const vector<vector<double>>& u_next,
                       const vector<vector<double>>& u){
        double sum = 0.0;
        int r = u.size(), c = u[0].size();

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                sum += (u_next[i][j] - u[i][j]) * (u_next[i][j] - u[i][j]);
            }
        } 

        return sqrt(sum);                      
    }

    void updateRed(vector<vector<double>>& u,
                   const vector<vector<double>>& f) {
        int r = u.size(), c = u[0].size();

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                if ((i + j) % 2) continue; // not red point
                double up = i > 0 ? u[i-1][j] : 0.0;
                double down = i < r-1 ? u[i+1][j] : 0.0;
                double left = j > 0 ? u[i][j-1] : 0.0;
                double right = j < c-1 ? u[i][j+1] : 0.0;
                u[i][j] = (h*h*f[i][j] + up + left + down + right) / 4;
            }
        }
    }

    void updateBlack(vector<vector<double>>& u,
                     const vector<vector<double>>& f) {
        int r = u.size(), c = u[0].size();

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                if ((i + j) % 2 == 0) continue; // not black point
                double up = i > 0 ? u[i-1][j] : 0.0;
                double down = i < r-1 ? u[i+1][j] : 0.0;
                double left = j > 0 ? u[i][j-1] : 0.0;
                double right = j < c-1 ? u[i][j+1] : 0.0;
                u[i][j] = (h*h*f[i][j] + up + left + down + right) / 4;
            }
        }
    }
    const int N;
    double h;
};

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("usage: program name N\n");
        exit(1);
    }

    int N = atoi(argv[1]);
    int maxIters = 5000;

    GS2D solver(N);
    Jacobi2D jaco_solver(N);

    // report timing and error 
    #ifdef _OPENMP
    double t_start = omp_get_wtime();
    #endif
    vector<vector<double>> u = solver.solution(maxIters);
    #ifdef _OPENMP
    double t_end = omp_get_wtime();
    printf("Time elapsed is %7.3f\n", t_end - t_start);
    #endif

    vector<vector<double>> u_jaco = jaco_solver.solution(maxIters);

    double square_distance = 0.0;

    for (int i = 0; i < (int) u.size(); ++i) {
        for (int j = 0; j < (int) u[i].size(); ++j) 
            square_distance += (u[i][j] - u_jaco[i][j]) * (u[i][j] - u_jaco[i][j]);
    }

    printf("L2 norm between two solutions: %10.5f\n", square_distance);

    return 0;
}