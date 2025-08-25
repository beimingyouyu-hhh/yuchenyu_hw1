#include "../inc/algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = k * a.data[i][j];
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

Matrix get_minor(Matrix a, int row, int col) {
    Matrix minor = create_matrix(a.rows - 1, a.cols - 1);
    int r = 0, c = 0;
    for (int i = 0; i < a.rows; i++) {
        if (i == row) continue;
        c = 0;
        for (int j = 0; j < a.cols; j++) {
            if (j == col) continue;
            minor.data[r][c] = a.data[i][j];
            c++;
        }
        r++;
    }
    return minor;
}

double det_matrix(Matrix a)
{
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    
    // 1x1矩阵直接返回元素值
    if (a.rows == 1) return a.data[0][0];
    
    // 2x2矩阵直接计算
    if (a.rows == 2) {
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    }
    
    // 大于2x2的矩阵使用递归计算
    double det = 0;
    int sign = 1;
    for (int j = 0; j < a.cols; j++) {
        Matrix minor = get_minor(a, 0, j);
        det += sign * a.data[0][j] * det_matrix(minor);
        sign = -sign;
    }
    return det;
}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    
    double det = det_matrix(a);
    if (fabs(det) < 1e-10) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    
    Matrix adj = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            Matrix minor = get_minor(a, i, j);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj.data[j][i] = sign * det_matrix(minor); // 注意这里是j,i而不是i,j
        }
    }
    
    return scale_matrix(adj, 1.0/det);
}

int rank_matrix(Matrix a)
{
    Matrix temp = create_matrix(a.rows, a.cols);
    // 复制原矩阵
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            temp.data[i][j] = a.data[i][j];
        }
    }
    
    int rank = (a.rows < a.cols) ? a.rows : a.cols;
    for (int i = 0; i < rank; i++) {
        // 如果对角线元素为0，寻找非零元素交换行
        if (fabs(temp.data[i][i]) < 1e-10) {
            int j;
            for (j = i + 1; j < a.rows; j++) {
                if (fabs(temp.data[j][i]) > 1e-10) {
                    // 交换行
                    for (int k = 0; k < a.cols; k++) {
                        double t = temp.data[i][k];
                        temp.data[i][k] = temp.data[j][k];
                        temp.data[j][k] = t;
                    }
                    break;
                }
            }
            if (j == a.rows) {
                rank--;
                continue;
            }
        }
        
        // 消元
        for (int j = i + 1; j < a.rows; j++) {
            double ratio = temp.data[j][i] / temp.data[i][i];
            for (int k = i; k < a.cols; k++) {
                temp.data[j][k] -= ratio * temp.data[i][k];
            }
        }
    }
    
    return rank;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    
    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}