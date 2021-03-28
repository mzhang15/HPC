#include <iostream>
#include <vector>

using namespace std;

class GS2D {
public:
    GS2D(int N) : N(N) {
        h = 1 / (N + 1);
    }

    vector<vector<double>> solution(int max_iters) {
        vector<vector<double>> u(N, vector<double>(N, 0.0));
        vector<vector<double>> f(N, vector<double>(N, 1.0));
        double res0 = computeResidual(u, f), res, tol = 1e-5;

        for (int i = 0; i < maxIters; ++i) {
            res = computeResidual(u, f);
            printf("Norm of residual at iteration %d = %10f\n", i, res);
            if (res / res0 <= tol) break;

            updateRed(u, f); 
            updateBlack(u, f);
        }

        return u;
    }
private:
    void updateRed(vector<vector<double>>& u,
                   const vector<vector<double>>& f) {
        int r = u.size(), c = u[0].size();

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                if (i + j % 2) continue; // not red point
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

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {

            }
        }
    }
    const int N;
    double h;
};