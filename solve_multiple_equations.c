#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// Giai he phuong trinh dai so tuyen tinh Ax = b

struct Matrix {
    int row;
    int col;
    double **data;
};

void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

// In ra ma tran
void printMatrix(struct Matrix matrix, int num_decimal) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        for (int j = 0; j < matrix.col; j++) {
            printf("%.*f ", num_decimal, matrix.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// In ra ma tran suy rong
void printMatrix2(struct Matrix matrix, int num_decimal) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        for (int j = 0; j < matrix.col-1; j++) {
            printf("%.*f ", num_decimal, matrix.data[i][j]);
        }
        printf("| %.*f\n", num_decimal, matrix.data[i][matrix.col-1]);
    }
}

// Ma tran bac thang rut gon thu duoc tu phep khu Gauss
struct Matrix duaVeMaTranBacThang(struct Matrix A1, int num_decimal) {
    // Chuyen doi ma tran suy rong A1 ve dang rut gon bac thang
    struct Matrix mat_ans = A1;
    int r = A1.row;
    int c = A1.col;
    for (int i = 0; i < r-1; i++) {
        printf("Ta chuan hoa cot %d: \n", i+1);
        printf("\tTa kiem tra A1[%d][%d]?\n", i+1, i+1);
        if( mat_ans.data[i][i] == 0 ){
            int check = 1;
            printf("\tDo A1[%d][%d] = 0 nen ta di tim hang i>%d ma co A1[i][%d] != 0.\n", i+1, i+1, i+1, i+1);
            for(int j = i+1; j < r; j++){
                if( mat_ans.data[j][i] != 0 ){
                    check = 0;
                    printf("Do A1[%d][%d] != 0 nen ta doi hai hang %d va %d cho nhau.\n", j+1, i+1, i+1, j+1);
                    for(int k = 0; k < c; k++){
                        swap(&mat_ans.data[i][k], &mat_ans.data[j][k]);
                    }
                    printMatrix2(mat_ans, num_decimal);
                    break;
                }
            }
            if( check ){
                printf("\tKhong co hang i>%d nao ma A1[i][%d] != 0 nen tat ca phan tu cot %d da duoc chuan hoa.\n", i+1, i+1, i+1);
                break;
            }
        }
        if( mat_ans.data[i][i] != 1){
            printf("\tDo A[%d][%d] != 1 nen ta tim hang i>%d ma co A1[i][%d] = 1 de giai don gian hon.\n", i+1, i+1, i+1, i+1);
            int check = 1;
            for(int j = i+1; j < r; j++){
                if( mat_ans.data[j][i] == 1 ){
                    check = 0;
                    printf("\tDo A[%d][%d] = 1 nen ta doi hai hang %d va %d cho nhau.\n", j+1, i+1, i+1, j+1);
                    for(int k = 0; k < c; k++){
                        swap(&mat_ans.data[i][k], &mat_ans.data[j][k]);
                    }
                    printMatrix2(mat_ans, num_decimal);
                    break;
                }
            }
            if( check ){
                printf("\tKhong co hang i>%d nao co A1[i][%d] = 1 nen ta tiep tuc chuan hoa nhu binh thuong.\n", i+1, i+1, i+1);
            }
        }
        int check = 0;  // Kiem tra xem cot i+1 da duoc chuan hoa chua
        for (int j = i + 1; j < r; j++) {
            if( mat_ans.data[j][i] != 0 ){
                check = 1;
            }
            float t = mat_ans.data[j][i] / mat_ans.data[i][i];
            for (int k = 0; k < c; k++) {
                mat_ans.data[j][k] -= t * mat_ans.data[i][k];
            }
        }
        if( check ){
            printf("\tTa thuc hien chuan hoa cot %d.\n", i+1, i+1, i+1);
            printf("\tMa tran A1 sau khi chuan hoa cot %d:\n", i + 1);
            printMatrix2(mat_ans, num_decimal);
        }else{
            printf("\tCot %d da duoc chuan hoa.\n", i+1);
        }
    }
    return mat_ans;
}

void tinhDinhThuc(struct Matrix A, struct Matrix A_bac_thang, int num_decimal){
    // Tinh dinh thuc
    printf("Tinh dinh thuc ma tran A:\n");
    printf("\tdet(A) = det(A_bac_thang) = A_bac_thang[1][1]*...*A_bac_thang[%d][%d] = ", A.row, A.row);
    double det_A = 1.0;
    for(int i = 0; i < A.row; i++){
        det_A *= A_bac_thang.data[i][i];
    }
    printf("%.*f\n\n", num_decimal, det_A);
}

// Giai he phuong trinh bang phuong phap Gauss
void giaiHePhuongTrinh(struct Matrix A, struct Matrix b, struct Matrix A_bac_thang, int num_decimal){
    // Chuyen doi ma tran A ve dang rut gon bac thang
    tinhDinhThuc(A, A_bac_thang, num_decimal);
    // Giai he phuong trinh
    struct Matrix x = A_bac_thang;
    for(int i = 0; i < A.row; i++){
        x.data[i][A.col] = b.data[i][0];
    }
    for(int i = A.row-1; i >= 0; i--){
        for(int j = i+1; j < A.col; j++){
            x.data[i][0] -= x.data[i][j] * x.data[j][0];
        }
        x.data[i][0] /= x.data[i][i];
    }
    printf("\nX = \n");
    printMatrix(x, num_decimal);
}

int main() {
    printf("+----------------------------------------------------------------------------------------------------------+\n");
    printf("|CHUONG TRINH GIAI HE PHUONG TRINH DAI SO TUYEN TINH Ax = B bang phuong phap Gauss, Gauss-Joicdan, Solepski|\n");
    printf("+----------------------------------------------------------------------------------------------------------+\n\n");
    int num_decimal;
    printf("Nhap so chu so thap phan se duoc hien thi: ");
    scanf("%d", &num_decimal);
    while (num_decimal < 0) {
        printf("Nhap sai! Vui long nhap lai: ");
        scanf("%d", &num_decimal);
    }
    // Nhap vao ma tran A
    struct Matrix A;
    {
        printf("Nhap vao kich co cua ma tran vuong A: ");
        scanf("%d", &A.row);
        while( A.row <= 0 ){
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
            printf("\t");
            for (int j = 0; j < A.col; j++) {
                scanf("%lf", &A.data[i][j]);
            }
        }
    }
    // Nhap vao ma tran b
    struct Matrix b;
    {
        b.row = A.row;
        b.col = 1;
        b.data = (double **)malloc(b.row * sizeof(double *));
        for (int i = 0; i < b.row; i++) {
            b.data[i] = (double *)malloc(b.col * sizeof(double));
        }
        printf("Nhap vao ma tran b kich co %dx%d: \n", b.row, b.col);
        for (int i = 0; i < b.row; i++) {
            printf("\t");
            for (int j = 0; j < b.col; j++) {
                scanf("%lf", &b.data[i][j]);
            }
        }
    }

    
    struct Matrix A1;
    {
        A1.row = A.row;
        A1.col = A.col+1;
        A1.data = (double **)malloc(A1.row * sizeof(double *));
        for (int i = 0; i < A1.row; i++) {
            A1.data[i] = (double *)malloc(A1.col * sizeof(double));
        }
        for (int i = 0; i < A1.row; i++) {
            for (int j = 0; j < A1.col-1; j++) {
                A1.data[i][j] = A.data[i][j];
            }
            A1.data[i][A1.col-1] = b.data[i][0];
        }
        printf("Ma tran suy rong:\n\t A1 = (A|b) = \n");
        printMatrix2(A1, num_decimal);
        printf("Ta dua ma tran suy rong A1 ve dang ma tran bac thang: \n");
        struct Matrix A_bac_thang = duaVeMaTranBacThang(A1, num_decimal);
        printf("Vay dang bac thang cua ma tran suy rong la: \n");
        printMatrix2(A_bac_thang, num_decimal);
    }

    // Tinh dinh thuc ma tran A bang phuong phap Gauss:
    printf("+------------------------------------------------------+\n");
    printf("|TINH DINH THUC MA TRAN VUONG A BANG PHUONG PHAP GAUSS:|\n");
    printf("+------------------------------------------------------+\n\n");
    // Dua A ve ma tran bac thang
    struct Matrix A_bac_thang;
    {
        A_bac_thang.row = A.row;
        A_bac_thang.col = A.col;
        A_bac_thang.data = (double **)malloc(A_bac_thang.row * sizeof(double *));
        for (int i = 0; i < A_bac_thang.row; i++) {
            A_bac_thang.data[i] = (double *)malloc(A_bac_thang.col * sizeof(double));
        }
        for(int i = 0; i < A_bac_thang.row; i++){
            for(int j = 0; j < A_bac_thang.col; j++){
                A_bac_thang.data[i][j] = A1.data[i][j];
            }
        }
    }
    printf("Dang bac thang cua ma tran A la:\nA_bac_thang = \n");
    printMatrix(A_bac_thang, num_decimal);
    tinhDinhThuc(A, A_bac_thang, num_decimal);
    
    return 0;
    // Giai he phuong trinh dai so bang phuong phap Gauss:
    printf("+---------------------------------------------------+\n");
    printf("|GIAI HE PHUONG TRINH DAI SO BANG PHUONG PHAP GAUSS:|\n");
    printf("+---------------------------------------------------+\n\n");
    printf("Ta su dung ma tran bac thang cua A o cau b de giai he phuong trinh dai so:\nA_bac_thang = \n");
    printMatrix(A_bac_thang, num_decimal);

    return 0;
}
