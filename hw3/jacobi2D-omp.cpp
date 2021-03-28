#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

using namespace std;

class Jacobi2D {
public:
    Jacobi2D(int N) : N(N) {
        h = (double) 1 / (N + 1);
        // printf("h = %7.3f\n", h);
    }

    vector<vector<double>> solution(int maxIters) {
        vector<vector<double>> u(N, vector<double>(N, 0.0)), u_next;
        vector<vector<double>> f(N, vector<double>(N, 1.0));
        double tol = 1e-5;

        u_next = u;
        for (int i = 0; i < maxIters; ++i) {
            update(u_next, u, f); 

            // print_matrix(u_next);

            if (square_loss(u_next, u) < tol) {
                u = u_next;
                break;
            }
            u = u_next;
        }

        return u;
    }
private:
    void print_matrix(const vector<vector<double>>& u) {
        int r = u.size(), c = u[0].size();
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) printf("%7.3f ", u[i][j]);
            cout << "\n";
        }
    }

    void update(vector<vector<double>>& u_next,
                const vector<vector<double>>& u,
                const vector<vector<double>>& f) {
        int r = u.size(), c = u[0].size();

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double up = i > 0 ? u[i - 1][j] : 0.0;
                double down = i < r-1 ? u[i + 1][j]: 0.0;
                double left = j > 0 ? u[i][j-1] : 0.0;
                double right = j < c-1 ? u[i][j + 1]: 0.0;
                u_next[i][j] = (h*h*f[i][j] + up + left + down + right) / 4;
            }
        }
    }

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
    
    const int N;
    double h;
};

// int main(int argc, char *argv[]) {
//     if (argc != 2) {
//         printf("usage: program name N\n");
//         exit(1);
//     }

//     int N = atoi(argv[1]);
//     int maxIters = 5000;

//     Jacobi2D solver(N);

//     #ifdef _OPENMP
//     double t_start = omp_get_wtime();
//     #endif
//     vector<vector<double>> u = solver.solution(maxIters);
//     #ifdef _OPENMP
//     double t_end = omp_get_wtime();
//     printf("Time elapsed is %7.3f\n", t_end - t_start);
//     #endif

//     return 0;
// }
