#include <iostream>
#include <vector>
#include <random>
#include <numeric>
#include <Eigen/Dense>

using namespace std;

const double eps = 1e-9;

int d;
vector<vector<double>> p;

double dot_product(const vector<double>& a, const vector<double>& b) {
    double res = 0;
    for (int i = 0; i < d; ++i) {
        res += a[i] * b[i];
    }
    return res;
}


double det(const vector<vector<double>>& m) {
    auto d = m.size();
    if (d == 1) {
        return m[0][0];
    }
    double res = 0;
    auto rows = m.size();
    auto cols = m[0].size();
    Eigen::MatrixXd matrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix(i, j) = m[i][j];
        }
    }
    return matrix.determinant();
}


vector<double> solve_linear_system(const vector<vector<double>>& A, const vector<double>& b) {
    vector<vector<double>> m(d, vector<double>(d + 1));
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            m[i][j] = A[i][j];
        }
        m[i][d] = b[i];
    }

    for (int col = 0; col < d; ++col) {
        int row_with_max_value = col;
        for (int row = col + 1; row < d; ++row) {
            if (abs(m[row][col]) > abs(m[row_with_max_value][col])) {
                row_with_max_value = row;
            }
        }
        swap(m[col], m[row_with_max_value]);

        for (int row = col + 1; row < d; ++row) {
            double c = m[row][col] / m[col][col];
            for (int j = col; j <= d; ++j) {
                m[row][j] -= c * m[col][j];
            }
        }
    }

    vector<double> x(d);
    for (int row = d - 1; row >= 0; --row) {
        x[row] = m[row][d];
        for (int col = row + 1; col < d; ++col) {
            x[row] -= x[col] * m[row][col];
        }
        x[row] /= m[row][row];
    }

    return x;
}

bool is_inside_simplex(const vector<vector<double>>& simplex, const vector<double>& point) {
    vector<vector<double>> A(d, vector<double>(d));
    vector<double> b(d);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            A[i][j] = simplex[j + 1][i] - simplex[0][i];
        }
        b[i] = point[i] - simplex[0][i];
    }

    auto coefs = solve_linear_system(A, b);

    double sum_coefs = accumulate(coefs.begin(), coefs.end(), 0.);

    return all_of(coefs.begin(), coefs.end(), [](double coef)
    {
        return coef >= -eps; }) && sum_coefs <= 1 + eps;
}

vector<vector<double>> find_convex_hull() {
    vector<vector<double>> res;
    while (p.size() >= d + 1) {
        mt19937 gen(random_device{}());
        uniform_real_distribution<> dis(-1, 1);

        vector<vector<double>> simplex(d + 1);
        bool found = false;
        while (!found) {
            vector<vector<double>> dirs(d + 1, vector<double>(d));
            for (int i = 0; i < d + 1; ++i) {
                for (int j = 0; j < d; ++j) {
                    dirs[i][j] = dis(gen);
                }
            }

            for (int i = 0; i < d + 1; ++i) {
                int max_point_index = -1;
                double max_dot_product = -INFINITY;
                for (int j = 0; j < p.size(); ++j) {
                    double dp = dot_product(dirs[i], p[j]);
                    if (dp > max_dot_product) {
                        max_dot_product = dp;
                        max_point_index = j;
                    }
                }
                simplex[i] = p[max_point_index];
            }

            vector<vector<double>> m(d, vector<double>(d));
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    m[i][j] = simplex[j + 1][i] - simplex[0][i];
                }
            }

            if (abs(det(m)) > eps) {
                found = true;
            }
        }

        for (const auto& point : simplex) {
            res.push_back(point);
        }

        p.erase(remove_if(p.begin(), p.end(), [&](const vector<double>& point) {
            return is_inside_simplex(simplex, point);
        }), p.end());
    }

    return res;
}
std::vector<std::vector<double>> generatePointsOnSphere(int numPoints) {
    std::vector<std::vector<double>> points;
    for (int i = 0; i < numPoints; ++i) {
        double y = 1 - (i / static_cast<double>(numPoints - 1)) * 2;
        double radius = 1;

        double theta = 0.1 * i;

        double x = cos(theta) * radius;
        double z = sin(theta) * radius;

        points.push_back({ x, y, z });
    }

    for (int i = 0; i < numPoints; ++i) {
        double y = 1 - (i / static_cast<double>(numPoints - 1)) * 2;
        double radius = 1;

        double theta = 0.1 * i;

        double x = cos(theta) * radius;
        double z = sin(theta) * radius;

        points.push_back({ x, y, z });
    }


    return points;
}




int main() {
    std::vector<std::vector<double>> points;


    /* 4-х Мерный куб( 1 точка осталась )
    d = 4;
    auto n = 18;
   
    std::vector<std::vector<double>> points;
    for (double i = 0; i <= 1; ++i)
        for (double j = 0; j <= 1; ++j)
            for (double k = 0; k <= 1; ++k)
                for (double l = 0; l <= 1; ++l) points.push_back({ i,j,k,l });
    points.push_back({ 0,0.1,0.2,0.3 });
    points.push_back({ 0,0.1,0.3,0.3 });
    points.push_back({ 0,0.1,0.2,0.4 });
    */

    /* Cфера 3-х (Не осталось внутр)
    d = 3;
    auto n = 10;
    points = generatePointsOnSphere(n);
    points.push_back({ 0,0.1,0.2,0.3 });
    points.push_back({ 0,0.1,0.3,0.3 });
    points.push_back({ 0,0.1,0.2,0.4 });
    */

    d = 3;
    auto n = 10;
    
    p.resize(n, vector<double>(d));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            p[i][j] = points[i][j];
        }
    }

    auto convex_hull = find_convex_hull();

    cout << convex_hull.size() << endl;
    for (const auto& point : convex_hull) {
        for (double coord : point) {
            cout << coord << ' ';
        }
        cout << endl;
    }

    return 0;
}