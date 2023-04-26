// Ryabov Aleksandr
// a.ryabov@innopolis.university
// B22-DSAI-01

#include <bits/stdc++.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>

using namespace std;

double limit(double i)
{
    return ( abs(i) < pow(10, -10) ) ? 0.0 : i; 
}

class Matrix
{
    public:
        int row;
        int column;
        vector<vector<double>> table;

        const string error_messange = {"Error: the dimensional problem occurred"};

        Matrix () 
        {
            this->table = vector<vector<double>> {};
        }

       friend ostream& operator<< (ostream& out, const Matrix& mtx)
        {
            for (auto row: mtx.table)
            {
                for (auto elem: row) out << elem << " " ;
                cout << "\n";
            }

            return out;
        }

        friend istream& operator >> (istream& in, Matrix& mtx)
        {
            in >> mtx.row >> mtx.column;

            mtx.column = mtx.row;
            mtx.table.clear();

            auto clm = vector<double> {};
            double temp;
            for (int i = 0; i < mtx.row; i++)
            {
                clm.clear();
                for (int j = 0; j < mtx.column; j++) 
                {
                    in >> temp;
                    clm.push_back(temp);
                }
                mtx.table.push_back(clm);
            }

            return in;
        }

        void operator= (const Matrix& mtx)
        {
            this->row = mtx.row;
            this->column = mtx.column;
            this->table = mtx.table;
        }

        Matrix operator+ (const Matrix& mtx)
        {
            Matrix variantial = *new Matrix();

            if (this->column != mtx.column) throw "Error: the dimensional problem occurred";
            if (this->row != mtx.row) throw "Error: the dimensional problem occurred";

            variantial.row = mtx.row;
            variantial.column = mtx.column;
            variantial.table.clear();

            auto clm = vector<double> {};
            for (int i = 0; i < mtx.row; i++)
            {
                clm.clear();
                for (int j = 0; j < mtx.column; j++)
                { 
                    clm.push_back(this->table[i][j] + mtx.table[i][j]);
                }
                variantial.table.push_back(clm);
            }

            return variantial;
        }


        Matrix operator- (const Matrix& mtx)
        {
            Matrix variantial = *new Matrix();

            if (this->column != mtx.column) throw "Error: the dimensional problem occurred";
            if (this->row != mtx.row) throw "Error: the dimensional problem occurred";
            
            variantial.row = mtx.row;
            variantial.column = mtx.column;
            variantial.table.clear();

            auto clm = vector<double> {};
            for (int i = 0; i < mtx.row; i++)
            {
                clm.clear();
                for (int j = 0; j < mtx.column; j++)
                { 
                    clm.push_back(this->table[i][j] - mtx.table[i][j]);
                }
                variantial.table.push_back(clm);
            }

            return variantial;
        }

        Matrix operator* (const Matrix& mtx)
        {
            Matrix variantial = *new Matrix();

            if (this->column != mtx.row) throw "Error: the dimensional problem occurred";

            variantial.row = this->row;
            variantial.column = mtx.column;
            variantial.table.clear();

            auto clm = vector<double> {}; 
            for (int i = 0 ; i < this->row; i++) 
            {
                clm.clear();

                for (int j = 0; j < mtx.column; j++) 
                {
                    double temp = 0;
                    for (int k = 0; k < mtx.row; k++) 
                    {
                        temp += (this->table[i][k] * mtx.table[k][j]);
                    }
                    clm.push_back(limit(temp));
                }

                variantial.table.push_back(clm);
            }

            return variantial;
        }

        Matrix transpose ()
        {
            auto variantial = *new Matrix();

            variantial.column = this->row;
            variantial.row = this->column;
            variantial.table.clear();

            auto clm = vector<double> {}; 
            for (int j = 0; j < this->column; j++)
            {
                clm.clear();
                for (int i = 0; i < this->row; i++)
                {
                    clm.push_back(this->table[i][j]);
                }
                variantial.table.push_back(clm);
            }

            return variantial;
        }

        void setSize (int r, int c)
        {
            this->column = c;
            this->row = r;
        }
};

class SquareMatrix: public Matrix
{
    public: 
        using Matrix::Matrix;

        friend void determinant(const Matrix& mtx);
        friend Matrix inverse(const Matrix& mtx);

        friend istream& operator >> (istream& in, SquareMatrix& mtx)
        {
            in >> mtx.row; // >> mtx.column;

            mtx.column = mtx.row;
            mtx.table.clear();

            auto clm = vector<double> {};
            double temp;
            for (int i = 0; i < mtx.row; i++)
            {
                clm.clear();
                for (int j = 0; j < mtx.column; j++) 
                {
                    in >> temp;
                    clm.push_back(temp);
                }
                mtx.table.push_back(clm);
            }

            return in;
        }

};

class IdentityMatrix: public SquareMatrix
{
    public:
        IdentityMatrix(int size)
        {
            this->column = size;
            this->row = size;
            this->table = vector<vector<double>>(this->row, vector<double>(this->column, 0));

            for (int i = 0; i < this->row; i++)
            {
                for (int j = 0; j < this->column; j++)
                {
                    if (i == j) this->table[i][j] = 1;
                }
            }

        }

};

class EliminationMatrix: public SquareMatrix
{
    public:
        EliminationMatrix(int p, int r, int rw, const Matrix& mtx)
        {
            //if (mtx.column != mtx.row) throw "NonSquareMatrixException";

            this->column = mtx.row;
            this->row = mtx.row;
            this->table = vector<vector<double>>(this->row, vector<double>(this->column, 0));

            for (int i = 0; i < this->row; i++)
            {
                for (int j = 0; j < this->column; j++)
                {
                    if (i == j) this->table[i][j] = 1;
                }
            }

            double alpha = -1 * (mtx.table[p - 1][r - 1] / mtx.table[rw - 1][r - 1]);

            this->table[p - 1][r - 1] = alpha;
        }
};

class PermutationMatrix: public SquareMatrix
{
    public:
        PermutationMatrix(int p, int r, const Matrix& mtx)
        {
            this->column = mtx.row;
            this->row = mtx.row;
            this->table = vector<vector<double>>(this->row, vector<double>(this->column, 0));

            for (int i = 0; i < this->row; i++)
            {
                for (int j = 0; j < this->column; j++)
                {
                    if (i == j && i != p - 1 && i != r - 1) this->table[i][j] = 1;
                }
            }

            this->table[r - 1][p - 1] = 1;
            this->table[p - 1][r - 1] = 1;

        }
};

class AugmentedMatrix: public Matrix
{
    public:
        AugmentedMatrix(const Matrix& mtx)
        {
            if (mtx.column != mtx.row) throw "NonSquareMatrixException";

            this->column = mtx.row + mtx.column;
            this->row = mtx.row;
            this->table = vector<vector<double>>(this->row, vector<double>(this->column, 0));

            for (int i = 0; i < this->row; i++)
            {
                for (int j = 0; j < this->column; j++)
                {
                    if (j < mtx.column) this->table[i][j] = mtx.table[i][j];
                    else 
                    {
                        if (i == j - mtx.column) this->table[i][j] = 1;
                    }

                }
            }

        }

        AugmentedMatrix(const Matrix& mtx1, const Matrix& mtx2)
        {
            if (mtx1.column != mtx1.row) throw "NonSquareMatrixException";
            
            this->column = mtx1.column + mtx2.column; 
            this->row    = mtx1.row;
            this->table = vector<vector<double>>(this->row, vector<double>(this->column, 0));

            for (int i = 0; i < this->row; i++)
            {
                for (int j = 0; j < this->column; j++)
                {
                    if (j < mtx1.column) this->table[i][j] = mtx1.table[i][j];
                    else this->table[i][j] =  mtx2.table[i][j - mtx1.column];
                }
            }
        }
};

class ColumnVector: public Matrix
{
    public:
        
        friend istream& operator >> (istream& in, ColumnVector& mtx)
        {
            in >> mtx.row; 

            mtx.column = 1;
            mtx.table.clear();

            auto clm = vector<double> {};
            double temp;
            for (int i = 0; i < mtx.row; i++)
            {
                clm.clear();
                for (int j = 0; j < mtx.column; j++) 
                {
                    in >> temp;
                    clm.push_back(temp);
                }
                mtx.table.push_back(clm);
            }

            return in;
        }

};


Matrix inverse(const Matrix& mtx)
{
    if (mtx.column != mtx.row) throw "NonSquareMatrixException";;

    // int step_cnt = 0;

    // cout << "step #" << step_cnt++ << ": Augmented Matrix" << endl;
    Matrix A = AugmentedMatrix(mtx);

    // cout << A << "Direct way:" << endl;

    auto& tb = A.table;

    for (int clm = 0; clm < mtx.column - 1; clm++)
    {
        
        int mx_row = clm;
        for (int i = clm + 1; i < mtx.row; i++)
        {
            if ( abs(tb[mx_row][clm]) < abs(tb[i][clm])) mx_row = i; 
        }    

        if (mx_row != clm) 
        {
            // cout << "step #" << step_cnt++ << ": permutation" << endl;
            auto P { PermutationMatrix(mx_row + 1, clm + 1, mtx) };

            A = P * A;

            // cout << A;
        }

        for (int i = clm + 1; i < mtx.row; i++)
        {
            if (tb[i][clm] == 0) continue;

            // cout << "step #" << step_cnt++ << ": elimination" << endl;
            // cout<< "###    " << i + 1 << "   " << clm + 1 << endl;
            auto E { EliminationMatrix(i + 1, clm + 1, clm + 1, A) };

            A = E * A;

            // cout << A;
        }
    
    }

    // cout << "Way back:" << endl;

    for (int clm = mtx.row - 1; clm >= 0; clm--)
    {
        for (int i = clm - 1; i >= 0; i--)
        {
            if (tb[i][clm] == 0) continue;

            // cout << "step #" << step_cnt++ << ": elimination" << endl;

            auto E { EliminationMatrix(i + 1, clm + 1, clm + 1, A) };
            
            A = E * A;

            // cout << A;
        }
    }

    // cout << "Diagonal normalization:" << endl;

    for (int i = 0; i < A.row; i++)
    {
        double alpha = A.table[i][i];

        for (int j = 0; j < A.column; j++)
        {
            A.table[i][j] = A.table[i][j] / alpha; 
        }
    }

    auto H { Matrix() };
    H.setSize(mtx.row, mtx.column);
    H.table = vector<vector<double>>(H.row, vector<double>(H.column, 0));

    for (int i = 0; i < A.row; i++)
    {
        for (int j = 0; j < A.column; j++)
        {
            if (j > mtx.row - 1) H.table[i][j - mtx.row] = A.table[i][j];
        }
    }

    return H;

}




int main()
{
    cout << std::fixed << std::setprecision(4);

    int n;
    cin >> n;

    vector<pair<int, int>> vec;
    int bf1, bf2;
    for (int i = 0; i < n; i++)
    {
        cin >> bf1 >> bf2;
        vec.push_back(pair<int, int>(bf1, bf2));
    }

    int degree;
    cin >> degree;

    auto A { Matrix() };
    A.setSize(n, degree + 1);
    A.table = vector<vector<double>>(A.row, vector<double>(A.column, 0));

    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < degree + 1; j++)
        {
            A.table[i][j] = pow(vec[i].first, j);
        }
    }

    cout << "A:" << endl;
    cout << A;


    cout << "A_T*A:" << endl;
    auto A_t = A.transpose() * A;
    cout << A_t;


    cout << "(A_T*A)^-1:" << endl;
    auto B { Matrix() };
    B = A_t;
    B = inverse(B);
    cout << B;


    cout << "A_T*b:" << endl;
    auto b { Matrix() };
    b.setSize(vec.size(), 1);
    b.table = vector<vector<double>>(b.row, vector<double>(b.column, 0));
    for (int i = 0; i < n; i++)
    {
        b.table[i][0] = vec[i].second;
    }    
    auto C = A.transpose() * b;
    cout  << C;

    cout << "x~:" << endl;
    auto x = B * C;
    cout << x;


    return 0;
}


