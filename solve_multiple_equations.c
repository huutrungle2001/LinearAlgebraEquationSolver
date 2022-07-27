#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// So ky tu ben trai dau phay duoc hien thi
const int INT_PART = 4;


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
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// In ra ma tran bo sung
void printMatrix2(struct Matrix matrix, int num_decimal) {
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        for (int j = 0; j < matrix.col - 1; j++) {
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j]);
        }
        printf("| %.*f\n", num_decimal, matrix.data[i][matrix.col - 1]);
    }
}

// In ra man hinh dinh thuc
void printMatrix3(struct Matrix matrix, int num_decimal, double multiples){
    for (int i = 0; i < matrix.row; i++) {
        printf("\t\t");
        if( i == (matrix.row+1)/2 && multiples != 1){
            printf("%*.*f * ", num_decimal + INT_PART, num_decimal, multiples);
        }else{
            printf("       ");
            for(int j = 0; j < num_decimal; j++){
                printf(" ");
            }
        }
        printf("| ");
        for (int j = 0; j < matrix.col; j++) {
            printf("%*.*f ", num_decimal + INT_PART, num_decimal, matrix.data[i][j]);
        }
        printf("|\n");
    }
}


int tinhHang(struct Matrix A) {
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


void tinhDinhThucBangGauss(struct Matrix A, int num_decimal){
    printf("+-------------------------------------------------------+\n");
    printf("| TINH DINH THUC MA TRAN VUONG A BANG PHUONG PHAP GAUSS |\n");
    printf("+-------------------------------------------------------+\n\n");
    double multiples = 1;
    double detA = 1;
    printf("Ta co: \ndet(A) =  \n");
    printMatrix3(A, num_decimal, multiples);
    struct Matrix A1;
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
    int r = A.row;
    int c = A.col;
    for (int i = 0; i < r-1 ; i++) {
        int check = 0;
        for(int j = i; j < r; j++){
            if(A1.data[j][i] != 0){
                check = 1;
                if( j != i ){
                    printf("\t=\t(H%d <-> H%d)  \n", i + 1, j + 1);
                    for(int k = 0; k < c; k++){
                        swap(&A1.data[i][k], &A1.data[j][k]);
                    }
                    multiples *= -1;
                    printMatrix3(A1, num_decimal, multiples);
                }
                break;
            }
        }
        if( !check ){
            detA *= 0;
            continue;
        }
        
        detA *= A1.data[i][i];
        int kt = 0;
        for(int j = i+1; j < r; j++){
            if( A1.data[j][i] != 0 ){
                kt = 1;
                double t = A1.data[j][i] / A1.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j+1, j+1, num_decimal, t, i+1);
                for(int k = 0; k < c; k++){
                    A1.data[j][k] -= t * A1.data[i][k];
                }
            }
        }
        if( kt ){
            printf("\n\t=\t\n");
            printMatrix3(A1, num_decimal, multiples);
        }
    }
    detA *= A1.data[r-1][r-1];
    printf("\t=\t%.*f", num_decimal, A1.data[0][0]);
    for(int i = 1; i < r; i++){
        printf(" * %.*f", num_decimal, A1.data[i][i]);
    }
    if( detA == 0 ){
        detA = abs(detA);
    }
    printf("\n\t=\t%.*f\n", num_decimal, detA);
    printf("\nVay det(A) = %.*f\n\n", num_decimal, detA);
}

// Giai he phuong trinh bang phuong phap Gauss
void giaiHePhuongTrinh(struct Matrix M, int num_decimal) {
    // Giai he phuong trinh tam giac tren
    int r = M.row;
    int c = M.col;
    double *x = (double *)malloc(sizeof(double) * (c-1));
    
    x[c - 2] = M.data[r - 1][c-1] / M.data[r - 1][c-2];
    for (int i = r - 2; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < c-1; j++) {
            sum += M.data[i][j] * x[j];
        }
        x[i] = (M.data[i][c-1] - sum) / M.data[i][i];
    }
    for(int i = 0; i < c-1; i++) {
        printf("\t\tX[%d] = %*.*f\n", i + 1, num_decimal + INT_PART, num_decimal, x[i]);
    }
}
void giaiHePhuongTrinhBangGauss(struct Matrix A, struct Matrix b, int num_decimal){
    printf("+---------------------------------------------+\n");
    printf("| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS |\n");
    printf("+---------------------------------------------+\n\n");
    printf("Truoc het ta bien doi ma tran bo sung M ve dang bac thang:\n");
    struct Matrix M;
    int r = A.row;
    int c = A.col;
    M.row = r, M.col = c+1;
    M.data = (double **)malloc(sizeof(double *) * M.row);
    for(int i = 0; i < M.row; i++){
        M.data[i] = (double *)malloc(sizeof(double) * M.col);
    }
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            M.data[i][j] = A.data[i][j];
        }
        M.data[i][c] = b.data[i][0];
    }
    printf("Ta co ma tran bo sung: \nM = (A | b) = \n");
    printMatrix2(M, num_decimal);
    for (int i = 0; i < r-1 ; i++) {
        int check = 0;
        for(int j = i; j < r; j++){
            if(M.data[j][i] != 0){
                check = 1;
                if( j != i ){
                    printf("(H%d <-> H%d)  ", i + 1, j + 1);
                    for(int k = 0; k <= c; k++){
                        swap(&M.data[i][k], &M.data[j][k]);
                    }
                    printf("\n\t===>\t\n");
                    printMatrix2(M, num_decimal);
                }
                break;
            }
        }
        if( !check ){
            continue;
        }
        int kt = 0;
        for(int j = i+1; j < r; j++){
            if( M.data[j][i] != 0 ){
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j+1, j+1, num_decimal, t, i+1);
                for(int k = 0; k <= c; k++){
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if( kt ){
            printf("\n\t===>\t\n");
            printMatrix2(M, num_decimal);
        }
    }

    printf("\nLuc nay he phuong trinh tro thanh: \n");
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j <= i; j++) {
            printf("\t");
        }
        printf("%*.*f*x[%d]", num_decimal + INT_PART, num_decimal, M.data[i][i], i + 1);
        for (int j = i + 1; j < M.col - 1; j++) {
            if (M.data[i][j] != 0) {
                printf(" + %*.*f*x[%d]", num_decimal + INT_PART, num_decimal, M.data[i][j], j + 1);
            } else {
                printf("\t");
            }
        }
        printf(" = %*.*f\n", num_decimal + INT_PART, num_decimal, M.data[i][M.col - 1]);
    }

    struct Matrix A_bac_thang;
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
    // Tinh hang cua ma tran A va ma tran bo sung M
    int rank_A, rank_M;
    rank_A = tinhHang(A_bac_thang);
    rank_M = tinhHang(M);
    printf("Ta co hang cua ma tran A la: r(A) = %d\n", rank_A);
    printf("Ta co hang cua ma tran bo sung M la: r(M) = %d\n", rank_M);
    printf("\n");
    printf("Ta kiem tra bang Dinh ly ve nghiem cua HPTTT: \n");
    if( rank_A == rank_M ){
        if( rank_A == A.row ){
            printf("\tDo r(A) = r(M) = so an, nen he phuong trinh co nghiem duy nhat.\n");
            printf("\tTa giai lan luot cac phuong trinh tu duoi len ta duoc: \n");
            giaiHePhuongTrinh(M, num_decimal);
        }else{
            printf("\tHe phuong trinh co vo so nghiem:\n");
        }
    }else{
        printf("\tHe phuong trinh vo nghiem:\n");
    }
    printf("\n");
}

void giaiHePhuongTrinhBangGaussJoocdan(struct Matrix A, struct Matrix b, int num_decimal){
    printf("+-----------------------------------------------------+\n");
    printf("| GIAI HE PHUONG TRINH BANG PHUONG PHAP GAUSS-JOOCDAN |\n");
    printf("+-----------------------------------------------------+\n\n");
    printf("Truoc het ta bien doi ma tran bo sung M ve dang chinh tac:\n");
    struct Matrix M;
    int r = A.row;
    int c = A.col;
    M.row = r, M.col = c+1;
    M.data = (double **)malloc(sizeof(double *) * M.row);
    for(int i = 0; i < M.row; i++){
        M.data[i] = (double *)malloc(sizeof(double) * M.col);
    }
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            M.data[i][j] = A.data[i][j];
        }
        M.data[i][c] = b.data[i][0];
    }
    printf("Ta co ma tran bo sung: \nM = (A | b) = \n");
    printMatrix2(M, num_decimal);
    for (int i = 0; i < r ; i++) {
        int check = 0;
        for(int j = i; j < r; j++){
            if(M.data[j][i] != 0){
                check = 1;
                if( j != i ){
                    printf("(H%d <-> H%d)  ", i + 1, j + 1);
                    for(int k = 0; k <= c; k++){
                        swap(&M.data[i][k], &M.data[j][k]);
                    }
                    printf("\n\t===>\t\n");
                    printMatrix2(M, num_decimal);
                }
                break;
            }
        }
        if( !check ){
            continue;
        }
        int kt = 0;
        if( M.data[i][i] != 1){
            kt = 1;
            double t = M.data[i][i];
            printf("(H%d <-> (1/%.*f)*H%d)  ", i + 1, num_decimal, t, i + 1);
            for(int k = 0; k <= c; k++){
                M.data[i][k] /= t;
            }
        }
        for(int j = i+1; j < r; j++){
            if( M.data[j][i] != 0 ){
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j+1, j+1, num_decimal, t, i+1);
                for(int k = 0; k <= c; k++){
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if( kt ){
            printf("\n\t===>\t\n");
            printMatrix2(M, num_decimal);
        }
    }
    for(int i = 0; i < r; i++){
        if( M.data[i][i] == 0 ){
            printf("\n\tHe phuong trinh vo nghiem hoac vo so nghiem.\n");
            return;
        }
    }
    for (int i = r-1; i >= 0 ; i--) {
        int check = 0;
        int kt = 0;
        for(int j = i-1; j >= 0; j--){
            if( M.data[j][i] != 0 ){
                kt = 1;
                double t = M.data[j][i] / M.data[i][i];
                printf("(H%d <-> H%d - %.*f*H%d)  ", j+1, j+1, num_decimal, t, i+1);
                for(int k = 0; k <= c; k++){
                    M.data[j][k] -= t * M.data[i][k];
                }
            }
        }
        if( kt ){
            printf("\n\t===>\t\n");
            printMatrix2(M, num_decimal);
        }
    }
    printf("Khi do nghiem duy nhat cua he phuong trinh la: \n");
    for(int i = 0; i < c-1; i++) {
        printf("\t\tX[%d] = %*.*f\n", i + 1, num_decimal + INT_PART, num_decimal, x[i]);
    }
    printf("\n");

}



int main() {
    printf("+------------------------------------------------------------------------------------------------------------+\n");
    printf("| CHUONG TRINH GIAI HE PHUONG TRINH DAI SO TUYEN TINH Ax = B bang phuong phap Gauss, Gauss-Joocdan, Solepski |\n");
    printf("+------------------------------------------------------------------------------------------------------------+\n\n");
    
    printf("+-------------------------------------------+\n");
    printf("| NHAP SO CHU SO THAP PHAN SE DUOC HIEN THI |\n");
    printf("+-------------------------------------------+\n\n");
    // Nhap vao so chu so thap phan se duoc hien thi
    int num_decimal;
    printf("Nhap so chu so thap phan se duoc hien thi: ");
    scanf("%d", &num_decimal);
    while (num_decimal < 0) {
        printf("Nhap sai! Vui long nhap lai: ");
        scanf("%d", &num_decimal);
    }

    printf("+-----------------------------------+\n");
    printf("| NHAP VAO A VA b DUOI DANG MA TRAN |\n");
    printf("+-----------------------------------+\n\n");
    // Nhap vao ma tran A
    struct Matrix A;
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
        printf("Ma tran A la: \n");
        printMatrix(A, num_decimal);
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
        printf("Nhap vao vecto cot b kich co %dx%d: \n", b.row, b.col);
        for (int i = 0; i < b.row; i++) {
            printf("\tb[%d] = ", i + 1);
            scanf("%lf", &b.data[i][0]);
        }
        printf("Vector cot b la: \n");
        printMatrix(b, num_decimal);
    }

    printf("+-------------------------+\n");
    printf("| He phuong trinh co dang |\n");
    printf("+-------------------------+\n\n");
    for (int i = 0; i < A.row; i++) {
        printf("\t%*.*f*x[%d]", num_decimal + INT_PART, num_decimal, A.data[i][0], 1);
        for (int j = 1; j < A.col; j++) {
            printf(" + %*.*f*x[%d] ", num_decimal + INT_PART, num_decimal, A.data[i][j], j + 1);
        }
        printf("= %*.*f\n", num_decimal + INT_PART, num_decimal, b.data[i][0]);
    }
    printf("\n");
    
    printf("+------------------------------------------------------------+\n");
    printf("| Xin moi lua chon chuong trinh                              |\n");
    printf("+------------------------------------------------------------+\n");
    printf("| 1. Tinh dinh thuc ma tran A bang phuong phap Gauss         |\n");
    printf("+------------------------------------------------------------+\n");
    printf("| 2. Giai phuong trinh Ax = b bang phuong phap Gauss         |\n");
    printf("+------------------------------------------------------------+\n");
    printf("| 3. Giai phuong trinh Ax = b bang phuong phap Gauss-Joocdan |\n");
    printf("+------------------------------------------------------------+\n");
    printf("| 4. Giai phuong trinh Ax = b bang phuong phap Solepski      |\n");
    printf("+------------------------------------------------------------+\n\n");

    int choice;
    char decision;
    do{
        printf("Nhap lua chon: ");
        scanf("%d", &choice);
        while (choice < 1 || choice > 4){
            printf("Nhap sai! Vui long nhap lai: ");
            scanf("%d", &choice);
        }
        switch (choice) {
            case 1:
                tinhDinhThucBangGauss(A, num_decimal);
                break;
            case 2:
                giaiHePhuongTrinhBangGauss(A, b, num_decimal);
                break;
            case 3:
                giaiHePhuongTrinhBangGaussJoocdan(A, b, num_decimal);
                break;
            case 4:
                //giaiHePhuongTrinh(A, num_decimal);
                break;
        }
        printf("Ban co muon tiep tuc khong? (y/n): ");
        fflush(stdin);
        scanf("%c", &decision);
        while(decision != 'y' && decision != 'n'){
            printf("Nhap sai! Vui long nhap lai: ");
            fflush(stdin);
            scanf("%c", &decision);
        }
    }while(decision == 'y');
    
    return 0;
}
