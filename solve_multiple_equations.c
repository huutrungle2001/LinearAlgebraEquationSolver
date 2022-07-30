#include <math.h>
#include <stdio.h>
#include <stdlib.h>



// So ky tu ben trai dau phay duoc hien thi
const int INT_PART = 4;

// Luu tru ma tran
typedef struct Matrix {
    int row;
    int col;
    double **data;
}Matrix;

// Doi gia tri 2 so thuc
void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

// In ra ma tran
void printMatrix(Matrix matrix, int num_decimal, FILE *fp) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        fprintf(fp, "\t\t");
        
        for (int j = 0; j < matrix.col; j++) {
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
            fprintf(fp, "%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
        }
        printf("\n");
        fprintf(fp, "\n");
    }
    printf("\n");
    fprintf(fp, "\n");
}

// In ra ma tran bo sung
void printMatrix2(Matrix matrix, int num_decimal, FILE *fp) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        fprintf(fp, "\t\t");
        for (int j = 0; j < matrix.col - 1; j++) {
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
            fprintf(fp, "%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
        }
        printf("| %.*f\n", num_decimal, matrix.data[i][matrix.col - 1] + 0.0);
        fprintf(fp, "| %.*f\n", num_decimal, matrix.data[i][matrix.col - 1] + 0.0);
    }
}

// In ra dang dinh thuc cua ma tran
void printMatrix3(Matrix matrix, int num_decimal, double multiples, FILE *fp) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        fprintf(fp, "\t\t");
        if (i == (matrix.row + 1) / 2 && multiples != 1) {
            printf("%*.*f * ", num_decimal + INT_PART, num_decimal, multiples);
            fprintf(fp, "%*.*f * ", num_decimal + INT_PART, num_decimal, multiples);
        } else {
            printf("       ");
            fprintf(fp, "       ");
            for (int j = 0; j < num_decimal; j++) {
                printf(" ");
                fprintf(fp, " ");
            }
        }
        printf("| ");
        fprintf(fp, "| ");
        for (int j = 0; j < matrix.col; j++) {
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
            fprintf(fp, "%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j] + 0.0);
        }
        printf("|\n");
        fprintf(fp, "|\n");
    }
}

// Tinh hang cua ma tran bac thang
int tinhHang(Matrix A) {
    int r = A.row;
    int c = A.col;
    int rank = 0;
    for (int i = 0; i < r; i++) {
        int check = 0;
        for (int j = 0; j < c; j++) {
            if (A.data[i][j] != 0) {
                check = 1;
                break;
            }
        }
        rank += check;
    }
    return rank;
}

// Tinh dinh thuc bang khu Gauss
void tinhDinhThucBangGauss(Matrix A, int num_decimal, FILE *fp) {
    printf("+-------------------------------------------------------+\n");
    printf("| TINH DINH THUC MA TRAN VUONG A BANG PHUONG PHAP GAUSS |\n");
    printf("+-------------------------------------------------------+\n\n");
    fprintf(fp, "+-------------------------------------------------------+\n");
    fprintf(fp, "| TINH DINH THUC MA TRAN VUONG A BANG PHUONG PHAP GAUSS |\n");
    fprintf(fp, "+-------------------------------------------------------+\n\n");
    // Thuc hien tich det A bang viec dua A ve dang ma tran bac thang
    // va tich cac phan tu tren duong cheo chinh
    double multiples = 1;
    double detA = 1;
    printf("Ta co: \ndet(A) =  \n");
    fprintf(fp, "Ta co: \ndet(A) =  \n");
    printMatrix3(A, num_decimal, multiples, fp);
    // Gan A sang mot ma tran tam thoi A1
    Matrix A1;
    A1.row = A.row;
    A1.col = A.col;
    A1.data = (double **)malloc(sizeof(double *) * A.row);
    for (int i = 0; i < A.row; i++) {
        A1.data[i] = (double *)malloc(sizeof(double) * A.col);
    }
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.col; j++) {
            A1.data[i][j] = A.data[i][j];
        }
    }

    // Thuc hien khu Gauss
    int r = A.row;
    int c = A.col;
    for (int i = 0; i < r - 1; i++) {
        int check = 0;
        // Thuc hien doi hang neu A[i][i] = 0
        for (int j = i; j < r; j++) {
            if (A1.data[j][i] != 0) {
                check = 1;
                if (j != i) {
                    printf("\t=\t(H%d <-> H%d)  \n", i + 1, j + 1);
                    fprintf(fp, "\t=\t(H%d <-> H%d)  \n", i + 1, j + 1);
                    for (int k = 0; k < c; k++) {
                        swap(&A1.data[i][k], &A1.data[j][k]);
                    }
                    multiples *= -1;
                    printMatrix3(A1, num_decimal, multiples, fp);
                }
                break;
            }
        }
        if (!check) {
            detA *= 0;
            continue;
        }
        // Bien doi co ban
        detA *= A1.data[i][i];
        int kt = 0;
        for (int j = i + 1; j < r; j++) {
            if (A1.data[j][i] != 0) {
                kt = 1;
                double t = A1.data[j][i] / A1.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                fprintf(fp, "(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                for (int k = 0; k < c; k++) {
                    A1.data[j][k] -= t * A1.data[i][k];
                }
            }
        }
        if (kt) {
            printf("\n\t=\t\n");
            fprintf(fp, "\n\t=\t\n");
            printMatrix3(A1, num_decimal, multiples, fp);
        }
    }
    detA *= A1.data[r - 1][r - 1];
    printf("\t=\t%.*f", num_decimal, A1.data[0][0] + 0.0);
    fprintf(fp, "\t=\t%.*f", num_decimal, A1.data[0][0] + 0.0);
    for (int i = 1; i < r; i++) {
        printf(" * %.*f", num_decimal, A1.data[i][i] + 0.0);
        fprintf(fp, " * %.*f", num_decimal, A1.data[i][i] + 0.0);
    }
    printf("\n\t=\t%.*f\n", num_decimal, detA + 0.0);
    fprintf(fp, "\n\t=\t%.*f\n", num_decimal, detA + 0.0);
    printf("\nVay det(A) = %.*f\n\n", num_decimal, detA + 0.0);
    fprintf(fp, "\nVay det(A) = %.*f\n\n", num_decimal, detA + 0.0);
}

// Giai he phuong trinh tam giac tren
void giaiHePhuongTrinh(Matrix M, int num_decimal, FILE *fp) {
    // Giai tu duoi len tren
    int r = M.row;
    int c = M.col;
    double *x = (double *)malloc(sizeof(double) * (c - 1));

    x[c - 2] = M.data[r - 1][c - 1] / M.data[r - 1][c - 2];
    for (int i = r - 2; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < c - 1; j++) {
            sum += M.data[i][j] * x[j];
        }
        x[i] = (M.data[i][c - 1] - sum) / M.data[i][i];
    }
    for (int i = 0; i < c - 1; i++) {
        printf("\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, x[i] + 0.0);
        fprintf(fp, "\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, x[i] + 0.0);
    }
}

// Giai he phuong trinh bang phuong phap khu Gauss
void giaiHePhuongTrinhBangGauss(Matrix A, Matrix b, int num_decimal, FILE *fp) {
    printf("+---------------------------------------------+\n");
    printf("| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS |\n");
    printf("+---------------------------------------------+\n\n");
    fprintf(fp, "+---------------------------------------------+\n");
    fprintf(fp, "| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS |\n");
    fprintf(fp, "+---------------------------------------------+\n\n");
    printf("Truoc het ta bien doi ma tran bo sung M ve dang bac thang:\n");
    fprintf(fp, "Truoc het ta bien doi ma tran bo sung M ve dang bac thang:\n");
    // Ta thuc hien bien doi ma tran bo sung ve dang bac thang
    Matrix M;
    int r = A.row;
    int c = A.col;
    M.row = r, M.col = c + 1;
    M.data = (double **)malloc(sizeof(double *) * M.row);
    for (int i = 0; i < M.row; i++) {
        M.data[i] = (double *)malloc(sizeof(double) * M.col);
    }
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            M.data[i][j] = A.data[i][j];
        }
        M.data[i][c] = b.data[i][0];
    }
    printf("Ta co ma tran bo sung: \nM = (A | b) = \n");
    fprintf(fp, "Ta co ma tran bo sung: \nM = (A | b) = \n");
    printMatrix2(M, num_decimal, fp);
    for (int i = 0; i < r - 1; i++) {
        int check = 0;
        for (int j = i; j < r; j++) {
            if (M.data[j][i] != 0) {
                check = 1;
                if (j != i) {
                    printf("(H%d <-> H%d)  ", i + 1, j + 1);
                    fprintf(fp, "(H%d <-> H%d)  ", i + 1, j + 1);
                    for (int k = 0; k <= c; k++) {
                        swap(&M.data[i][k], &M.data[j][k]);
                    }
                    printf("\n\t===>\t\n");
                    fprintf(fp, "\n\t===>\t\n");
                    printMatrix2(M, num_decimal, fp);
                }
                break;
            }
        }
        if (!check) {
            continue;
        }
        int kt = 0;
        for (int j = i + 1; j < r; j++) {
            if (M.data[j][i] != 0) {
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                fprintf(fp, "(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                for (int k = 0; k <= c; k++) {
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if (kt) {
            printf("\n\t===>\t\n");
            fprintf(fp, "\n\t===>\t\n");
            printMatrix2(M, num_decimal, fp);
        }
    }

    printf("\nLuc nay he phuong trinh tro thanh: \n");
    fprintf(fp, "\nLuc nay he phuong trinh tro thanh: \n");
    for (int i = 0; i < M.row; i++) {
        printf("\t");
        fprintf(fp, "\t");
        for (int j = 0; j < i; j++) {
            for(int k = 0; k < 2*INT_PART+7+num_decimal; k++)
                printf(" ");
                fprintf(fp, " ");
        }
        printf("%*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, M.data[i][i] + 0.0, INT_PART, i + 1);
        fprintf(fp, "%*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, M.data[i][i] + 0.0, INT_PART, i + 1);
        for (int j = i + 1; j < M.col - 1; j++) {
            if (M.data[i][j] != 0) {
                printf(" + %*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, M.data[i][j] + 0.0, INT_PART, j + 1);
                fprintf(fp, " + %*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, M.data[i][j] + 0.0, INT_PART, j + 1);
            } else {
                printf("\t");
                fprintf(fp, "\t");
            }
        }
        printf(" = %*.*f\n", num_decimal + INT_PART, num_decimal, M.data[i][M.col - 1] + 0.0);
        fprintf(fp, " = %*.*f\n", num_decimal + INT_PART, num_decimal, M.data[i][M.col - 1] + 0.0);
    }

    Matrix A_bac_thang;
    A_bac_thang.row = A.row;
    A_bac_thang.col = A.col;
    A_bac_thang.data = (double **)malloc(A_bac_thang.row * sizeof(double *));
    for (int i = 0; i < A_bac_thang.row; i++) {
        A_bac_thang.data[i] = (double *)malloc(A_bac_thang.col * sizeof(double));
    }
    for (int i = 0; i < A_bac_thang.row; i++) {
        for (int j = 0; j < A_bac_thang.col; j++) {
            A_bac_thang.data[i][j] = M.data[i][j];
        }
    }
    printf("\n");
    fprintf(fp, "\n");
    // Tinh hang cua ma tran A va ma tran bo sung M
    // Thuc hien kiem tra nghiem cua phuong trinh dua vao dinh ly ve nghiem cua HPTTT

    int rank_A, rank_M;
    rank_A = tinhHang(A_bac_thang);
    rank_M = tinhHang(M);
    printf("Ta co hang cua ma tran A la: r(A) = %d\n", rank_A);
    fprintf(fp, "Ta co hang cua ma tran A la: r(A) = %d\n", rank_A);
    printf("Ta co hang cua ma tran bo sung M la: r(M) = %d\n", rank_M);
    fprintf(fp, "Ta co hang cua ma tran bo sung M la: r(M) = %d\n", rank_M);
    printf("\n");
    fprintf(fp, "\n");
    printf("Ta kiem tra bang Dinh ly ve nghiem cua HPTTT: \n");
    fprintf(fp, "Ta kiem tra bang Dinh ly ve nghiem cua HPTTT: \n");
    if (rank_A == rank_M) {
        if (rank_A == A.row) {
            printf("\tDo r(A) = r(M) = so an, nen he phuong trinh co nghiem duy nhat.\n");
            fprintf(fp, "\tDo r(A) = r(M) = so an, nen he phuong trinh co nghiem duy nhat.\n");
            printf("\tTa giai lan luot cac phuong trinh tu duoi len ta duoc: \n");
            fprintf(fp, "\tTa giai lan luot cac phuong trinh tu duoi len ta duoc: \n");
            // Giai he phuong trinh tam giac tren de rut ra nghiem cua phuong trinh neu no co nghiem duy nhat
            giaiHePhuongTrinh(M, num_decimal, fp);
        } else {
            printf("\tDo r(A) = r(M) < so an, nen he phuong trinh co vo so nghiem:\n");
            fprintf(fp, "\tDo r(A) = r(M) < so an, nen he phuong trinh co vo so nghiem:\n");
        }
    } else {
        printf("\tDo r(A) < r(M), nen he phuong trinh vo nghiem:\n");
        fprintf(fp, "\tDo r(A) < r(M), nen he phuong trinh vo nghiem:\n");
    }
    printf("\n");
}

// Giai he phuong trinh bang Gauss-Joocdan
void giaiHePhuongTrinhBangGaussJoocdan(Matrix A, Matrix b, int num_decimal, FILE *fp) {
    printf("+-----------------------------------------------------+\n");
    printf("| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS-JOOCDAN |\n");
    printf("+-----------------------------------------------------+\n\n");
    fprintf(fp, "+-----------------------------------------------------+\n");
    fprintf(fp, "| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS-JOOCDAN |\n");
    fprintf(fp, "+-----------------------------------------------------+\n\n");
    printf("Truoc het ta bien doi ma tran bo sung M ve dang chinh tac:\n");
    fprintf(fp, "Truoc het ta bien doi ma tran bo sung M ve dang chinh tac:\n");
    // Dua ma tran bo sung M ve dang chinh tac
    Matrix M;
    int r = A.row;
    int c = A.col;
    M.row = r, M.col = c + 1;
    M.data = (double **)malloc(sizeof(double *) * M.row);
    for (int i = 0; i < M.row; i++) {
        M.data[i] = (double *)malloc(sizeof(double) * M.col);
    }
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            M.data[i][j] = A.data[i][j];
        }
        M.data[i][c] = b.data[i][0];
    }
    printf("Ta co ma tran bo sung: \nM = (A | b) = \n");
    fprintf(fp, "Ta co ma tran bo sung: \nM = (A | b) = \n");
    printMatrix2(M, num_decimal, fp);
    // Truoc het dua ma tran M ve dang bac thang
    for (int i = 0; i < r; i++) {
        int check = 0;
        for (int j = i; j < r; j++) {
            if (M.data[j][i] != 0) {
                check = 1;
                if (j != i) {
                    printf("(H%d <-> H%d)  ", i + 1, j + 1);
                    fprintf(fp, "(H%d <-> H%d)  ", i + 1, j + 1);
                    for (int k = 0; k <= c; k++) {
                        swap(&M.data[i][k], &M.data[j][k]);
                    }
                    printf("\n\t===>\t\n");
                    fprintf(fp, "\n\t===>\t\n");
                    printMatrix2(M, num_decimal, fp);
                }
                break;
            }
        }
        if (!check) {
            continue;
        }
        int kt = 0;
        if (M.data[i][i] != 1) {
            kt = 1;
            double t = M.data[i][i];
            printf("(H%d <-> (1/%.*f)*H%d)  ", i + 1, num_decimal, t, i + 1);
            fprintf(fp, "(H%d <-> (1/%.*f)*H%d)  ", i + 1, num_decimal, t, i + 1);
            for (int k = 0; k <= c; k++) {
                M.data[i][k] /= t;
            }
        }
        for (int j = i + 1; j < r; j++) {
            if (M.data[j][i] != 0) {
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                fprintf(fp, "(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                for (int k = 0; k <= c; k++) {
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if (kt) {
            printf("\n\t===>\t\n");
            fprintf(fp, "\n\t===>\t\n");
            printMatrix2(M, num_decimal, fp);
        }
    }
    int rank_M = tinhHang(M);
    int rank_A = 0;
    for (int i = 0; i < r; i++) {
        int check = 0;
        for (int j = 0; j < c; j++) {
            if (M.data[i][j] != 0) {
                check = 1;
                break;
            }
        }
        rank_A += check;
    }
    // Thuc hien tinh hang va kiem tra HPT co nghiem duy nhat khong bang dinh ly ve nghiem cua HPTTT
    if( rank_A < rank_M ){
        printf("\n\tDo r(A) < r(M), nen he phuong trinh vo nghiem.\n");
        fprintf(fp, "\n\tDo r(A) < r(M), nen he phuong trinh vo nghiem.\n");
        return;
    }else if( rank_A == rank_M && rank_A < c ){
        printf("\n\tDo r(A) = r(M) < so an, nen he phuong trinh co vo so nghiem.\n");
        fprintf(fp, "\n\tDo r(A) = r(M) < so an, nen he phuong trinh co vo so nghiem.\n");
        return;
    }
    // Dua A ve dang chinh tac voi ve trai la ma tran don vi, ve phai la vector nghiem
    for (int i = r - 1; i >= 0; i--) {
        int check = 0;
        int kt = 0;
        for (int j = i - 1; j >= 0; j--) {
            if (M.data[j][i] != 0) {
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                fprintf(fp, "(H%d <-> H%d - %.*f*H%d)  ", j + 1, j + 1, num_decimal, t, i + 1);
                for (int k = 0; k <= c; k++) {
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if (kt) {
            printf("\n\t===>\t\n");
            fprintf(fp, "\n\t===>\t\n");
            printMatrix2(M, num_decimal, fp);
        }
    }
    printf("Khi do nghiem duy nhat cua he phuong trinh la: \n");
    fprintf(fp, "Khi do nghiem duy nhat cua he phuong trinh la: \n");
    for (int i = 0; i < r; i++) {
        printf("\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, M.data[i][c] + 0.0);
        fprintf(fp, "\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, M.data[i][c] + 0.0);
    }
    printf("\n");
    fprintf(fp, "\n");
}

// Kiem tra ma tran co phai la ma tran doi xung hay khong
int check_symmetric(Matrix A, FILE *fp) {
    printf("Ta kiem tra ma tran A co la ma tran doi xung hay khong: \n");
    fprintf(fp, "Ta kiem tra ma tran A co la ma tran doi xung hay khong: \n");
    if (A.row == A.col) {
        for (int i = 0; i < A.row; i++) {
            for (int j = 0; j < A.col; j++) {
                if (A.data[i][j] != A.data[j][i]) {
                    return 0;
                }
            }
        }
        return 1;
    } else {
        return 0;
    }
}

// Giai he phuong trinh bang phuong phap Choleski
void giaiHePhuongTrinhBangCholeski(Matrix A, Matrix b, int num_decimal, FILE *fp) {
    // A = L.L^T
    // L.L^T.x = b
    // Giai he phuong trinh
    // L^t.x = y
    // L.y = b
    printf("+-----------------------------------------------------+\n");
    printf("| GIAI HE PHUONG TRINH BANG PHUONG PHAP Choleski      |\n");
    printf("+-----------------------------------------------------+\n\n");
    printf("De A co khai trien Choleski: A = L.L^T, voi T la ma tran tam giac duoi thi A phai la ma tran duong-doi xung.\n");
    fprintf(fp, "+-----------------------------------------------------+\n");
    fprintf(fp, "| GIAI HE PHUONG TRINH BANG PHUONG PHAP Choleski      |\n");
    fprintf(fp, "+-----------------------------------------------------+\n\n");
    fprintf(fp, "De A co khai trien Choleski: A = L.L^T, voi T la ma tran tam giac duoi thi A phai la ma tran duong-doi xung.\n");

    // Truoc het de dua A ve dang khai tren Choleski
    // Ta kiem tra A co phai la ma tran duong-doi xung hay khong.
    // Dong thoi tim ma tran tam giac duoi L trong khai trien (neu co)

    if (!check_symmetric(A, fp)) {
        printf("\tMa tran A khong phai la ma tran duong doi xung. Vui long thu cach khac!\n");
        fprintf(fp, "\tMa tran A khong phai la ma tran duong doi xung. Vui long thu cach khac!\n");
        return;
    }
    int r = A.row;
    int c = A.col;
    Matrix lower_triangular;
    lower_triangular.row = r;
    lower_triangular.col = c;
    lower_triangular.data = (double **)malloc(sizeof(double *) * lower_triangular.row);
    for (int i = 0; i < r; i++) {
        lower_triangular.data[i] = (double *)malloc(sizeof(double) * lower_triangular.col);
        for (int j = 0; j < c; j++) {
            lower_triangular.data[i][j] = 0;
        }
    }

    for (int i = 0; i < r; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;

            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += pow(lower_triangular.data[i][k], 2);
                }
                if (A.data[i][i] - sum < 0) {
                    printf("\tMa tran A khong phai la ma tran duong doi xung. Vui long thu cach khac!\n");
                    fprintf(fp, "\tMa tran A khong phai la ma tran duong doi xung. Vui long thu cach khac!\n");
                    return;
                }
                lower_triangular.data[i][i] = sqrt(A.data[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++) {
                    sum += (lower_triangular.data[i][k] * lower_triangular.data[j][k]);
                }
                lower_triangular.data[i][j] = (A.data[i][j] - sum) / lower_triangular.data[j][j];
            }
        }
    }
    Matrix upper_triangular;
    upper_triangular.row = r;
    upper_triangular.col = c;
    upper_triangular.data = (double **)malloc(sizeof(double *) * upper_triangular.row);
    for (int i = 0; i < r; i++) {
        upper_triangular.data[i] = (double *)malloc(sizeof(double) * upper_triangular.col);
        for (int j = 0; j < c; j++) {
            upper_triangular.data[i][j] = lower_triangular.data[j][i];
        }
    }
    printf("Ta co khai trien Choleski cua A: A = L.L^T = \n");
    fprintf(fp, "Ta co khai trien Choleski cua A: A = L.L^T = \n");
    printMatrix(lower_triangular, num_decimal, fp);
    printf("\t\t*\n");
    printMatrix(upper_triangular, num_decimal, fp);

    // Giai lan luot hai phuong trinh tam giac duoi: Ly = b va tam giac tren L^t.x = y
    // Va dua ra nghiem cua he phuong trinh

    printf("Nghiem cua he phuong trinh L.y = b la: \n");
    fprintf(fp, "Nghiem cua he phuong trinh L.y = b la: \n");
    // Giai he phuong trinh tam giac duoi
    double *y = (double *)malloc(sizeof(double) * c);

    y[0] = b.data[0][0] / lower_triangular.data[0][0];
    for (int i = 1; i < r; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += lower_triangular.data[i][j] * y[j];
        }
        y[i] = (b.data[i][0] - sum) / lower_triangular.data[i][i];
    }
    for (int i = 0; i < c; i++) {
        printf("\t\ty[%d] = %*.*f\n", i + 1, num_decimal + INT_PART, num_decimal, y[i] + 0.0);
        fprintf(fp, "\t\ty[%d] = %*.*f\n", i + 1, num_decimal + INT_PART, num_decimal, y[i] + 0.0);
    }
    printf("\nNghiem cua he phuong trinh L^T.x = y la: \n");
    fprintf(fp, "\nNghiem cua he phuong trinh L^T.x = y la: \n");
    // Giai he phuong trinh tam giac tren
    double *x = (double *)malloc(sizeof(double) * c);

    x[c - 1] = y[r - 1] / upper_triangular.data[r - 1][c - 1];
    for (int i = r - 2; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < c; j++) {
            sum += upper_triangular.data[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / upper_triangular.data[i][i];
    }
    for (int i = 0; i < c; i++) {
        printf("\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, x[i] + 0.0);
        fprintf(fp, "\t\tx[%*d] = %*.*f\n", INT_PART, i + 1, num_decimal + INT_PART, num_decimal, x[i] + 0.0);
    }
}

int main() {
    FILE *fp;
    fp = fopen("output.txt", "w");
    printf("+------------------------------------------------------------------------------------------------------------+\n");
    printf("| CHUONG TRINH GIAI HE PHUONG TRINH DAI SO TUYEN TINH Ax = B bang phuong phap Gauss, Gauss-Joocdan, Choleski |\n");
    printf("+------------------------------------------------------------------------------------------------------------+\n\n");
    fprintf(fp, "+------------------------------------------------------------------------------------------------------------+\n");
    fprintf(fp, "| CHUONG TRINH GIAI HE PHUONG TRINH DAI SO TUYEN TINH Ax = B bang phuong phap Gauss, Gauss-Joocdan, Choleski |\n");
    fprintf(fp, "+------------------------------------------------------------------------------------------------------------+\n\n");

    printf("+-------------------------------------------+\n");
    printf("| NHAP SO CHU SO THAP PHAN SE DUOC HIEN THI |\n");
    printf("+-------------------------------------------+\n\n");
    fprintf(fp, "+-------------------------------------------+\n");
    fprintf(fp, "| NHAP SO CHU SO THAP PHAN SE DUOC HIEN THI |\n");
    fprintf(fp, "+-------------------------------------------+\n\n");
    // Nhap vao so chu so thap phan se duoc hien thi
    int num_decimal;
    printf("Nhap so chu so thap phan se duoc hien thi: ");
    scanf("%d", &num_decimal);
    while (num_decimal < 0) {
        printf("Nhap sai! Vui long nhap lai: ");
        scanf("%d", &num_decimal);
    }
    printf("So chu so thap phan duoc hien thi la: %d\n", num_decimal);
    fprintf(fp, "So chu so thap phan duoc hien thi la: %d\n", num_decimal);
    printf("+-----------------------------------+\n");
    printf("| NHAP VAO A VA b DUOI DANG MA TRAN |\n");
    printf("+-----------------------------------+\n\n");
    fprintf(fp, "+-----------------------------------+\n");
    fprintf(fp, "| NHAP VAO A VA b DUOI DANG MA TRAN |\n");
    fprintf(fp, "+-----------------------------------+\n\n");
    // Nhap vao ma tran A
    Matrix A;
    {
        printf("Nhap vao kich co cua ma tran vuong A: ");
        scanf("%d", &A.row);
        while (A.row <= 0) {
            printf("Nhap sai! Vui long nhap lai: ");
            scanf("%d", &A.row);
        }
        A.col = A.row;

        A.data = (double **)malloc(A.row * sizeof(double *));
        for (int i = 0; i < A.row; i++) {
            A.data[i] = (double *)malloc(A.col * sizeof(double));
        }
        printf("Nhap vao ma tran vuong A kich co %dx%d: \n", A.row, A.col);
        for (int i = 0; i < A.row; i++) {
            for (int j = 0; j < A.col; j++) {
                printf("\ta[%d][%d] = ", i + 1, j + 1);
                scanf("%lf", &A.data[i][j]);
            }
        }
        printf("Ma tran A kich co %dx%d la: \n", A.row, A.col);
        fprintf(fp,"Ma tran A kich co %dx%d la: \n", A.row, A.col);
        printMatrix(A, num_decimal, fp);
    }
    // Nhap vao ma tran b
    Matrix b;
    {
        b.row = A.row;
        b.col = 1;
        b.data = (double **)malloc(b.row * sizeof(double *));
        for (int i = 0; i < b.row; i++) {
            b.data[i] = (double *)malloc(b.col * sizeof(double));
        }
        printf("Nhap vao vecto cot b kich co %dx%d: \n", b.row, b.col);
        for (int i = 0; i < b.row; i++) {
            printf("\tb[%d] = ", i + 1);
            scanf("%lf", &b.data[i][0]);
        }
        printf("Vector cot b kich co %dx%d la: \n", b.row, b.col);
        printMatrix(b, num_decimal, fp);
    }

    
    // Dang cua he phuong trinh
    printf("+-------------------------+\n");
    printf("| He phuong trinh co dang |\n");
    printf("+-------------------------+\n\n");
    fprintf(fp, "+-------------------------+\n");
    fprintf(fp, "| He phuong trinh co dang |\n");
    fprintf(fp, "+-------------------------+\n\n");
    for (int i = 0; i < A.row; i++) {
        printf("\t%*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, A.data[i][0] + 0.0, INT_PART, 1);
        fprintf(fp,"\t%*.*f*x[%*d]", num_decimal + INT_PART, num_decimal, A.data[i][0] + 0.0, INT_PART, 1);
        for (int j = 1; j < A.col; j++) {
            printf(" + %*.*f*x[%*d] ", num_decimal + INT_PART, num_decimal, A.data[i][j] + 0.0, INT_PART, j + 1);
            fprintf(fp," + %*.*f*x[%*d] ", num_decimal + INT_PART, num_decimal, A.data[i][j] + 0.0, INT_PART, j + 1);
        }
        printf("= %*.*f\n", num_decimal + INT_PART, num_decimal, b.data[i][0] + 0.0);
        fprintf(fp, "= %*.*f\n", num_decimal + INT_PART, num_decimal, b.data[i][0] + 0.0);
    }
    printf("\n");
    fprintf(fp, "\n");

    // Thuc hien chuong trinh chinh
    int choice;
    char decision;
    do {
        printf("+------------------------------------------------------------+\n");
        printf("| Xin moi lua chon chuong trinh                              |\n");
        printf("+------------------------------------------------------------+\n");
        printf("| 1. Tinh dinh thuc ma tran A bang phuong phap Gauss         |\n");
        printf("+------------------------------------------------------------+\n");
        printf("| 2. Giai phuong trinh Ax = b bang phuong phap Gauss         |\n");
        printf("+------------------------------------------------------------+\n");
        printf("| 3. Giai phuong trinh Ax = b bang phuong phap Gauss-Joocdan |\n");
        printf("+------------------------------------------------------------+\n");
        printf("| 4. Giai phuong trinh Ax = b bang phuong phap Choleski      |\n");
        printf("+------------------------------------------------------------+\n\n");
        fprintf(fp,"+------------------------------------------------------------+\n");
        fprintf(fp, "| Xin moi lua chon chuong trinh                              |\n");
        fprintf(fp, "+------------------------------------------------------------+\n");
        fprintf(fp, "| 1. Tinh dinh thuc ma tran A bang phuong phap Gauss         |\n");
        fprintf(fp, "+------------------------------------------------------------+\n");
        fprintf(fp, "| 2. Giai phuong trinh Ax = b bang phuong phap Gauss         |\n");
        fprintf(fp,"+------------------------------------------------------------+\n");
        fprintf(fp,"| 3. Giai phuong trinh Ax = b bang phuong phap Gauss-Joocdan |\n");
        fprintf(fp,"+------------------------------------------------------------+\n");
        fprintf(fp,"| 4. Giai phuong trinh Ax = b bang phuong phap Choleski      |\n");
        fprintf(fp,"+------------------------------------------------------------+\n\n");
        printf("Nhap lua chon: ");
        scanf("%d", &choice);
        while (choice < 1 || choice > 4) {
            printf("Nhap sai! Vui long nhap lai: ");
            scanf("%d", &choice);
        }
        switch (choice) {
            case 1:
                tinhDinhThucBangGauss(A, num_decimal, fp);
                break;
            case 2:
                giaiHePhuongTrinhBangGauss(A, b, num_decimal, fp);
                break;
            case 3:
                giaiHePhuongTrinhBangGaussJoocdan(A, b, num_decimal, fp);
                break;
            case 4:
                giaiHePhuongTrinhBangCholeski(A, b, num_decimal, fp);
                break;
        }
        printf("Ban co muon tiep tuc khong? (y/n): ");
        fflush(stdin);
        scanf("%c", &decision);
        while (decision != 'y' && decision != 'n') {
            printf("Nhap sai! Vui long nhap lai: ");
            fflush(stdin);
            scanf("%c", &decision);
        }
    } while (decision == 'y');
    printf("CHUONG TRINH KET THUC! XIN CAM ON!\n");
    fprintf(fp, "CHUONG TRINH KET THUC! XIN CAM ON!\n");
    return 0;
}