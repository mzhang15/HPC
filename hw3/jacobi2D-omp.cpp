#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class Jacobi2D {
public:
    Jacobi2D(int N) : N(n) {
        h = 1 / (N + 1);
    }

    vector<vector<double>> solution(int maxIters) {
        vector<vector<double>> u(N, vector<double>(N, 0.0)), u_next;
        vector<vector<double>> f(N, vector<double>(N, 1.0));
        double tol = 1e-5;

        u_next = u;
        for (int i = 0; i < maxIters; ++i) {
            update(u_next, u, f); 
            if (square_loss(u_next, u) < tol) break;
            u = u_next;
        }

        return u;
    }
private:
    void update(vector<vector<double>>& u_next,
                const vector<vector<double>>& u,
                const vector<vector<double>>& f) {
        int r = u.size(), c = u[0].size();
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

    double square_loss(const vector<vector<double>>>& u_next,
                       const vectotr<vector<double>>& u){
        double sum = 0.0;
        int r = u.size(), c = u[0].size();

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                sum += (u_next[i][j] - u[i][j]) * (u_next[i][j] - u[i][j])
            }
        } 

        return sqrt(sum);                      
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

    Jacobi2D solver(N);

    vector<vector<double>> u = solver.solution(maxIters);

    for (int i = 0; i < u.size(); ++i) {
        for (int j = 0; j < u[i].size(); ++j) cout << u[i][j] << ' ';
        cout << endl;
    }

    return 0;
}
