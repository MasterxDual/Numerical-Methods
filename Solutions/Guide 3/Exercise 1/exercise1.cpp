#include <cstdio>
#include <stdlib.h>
#include <cmath>

// Now you can handle up to 50x50 matrices
#define MAX_SIZE 50  

/**
 * Reads an augmented matrix from a file text called data.dat
 * The expected format is:
 *   a11 a12 ... a1n b1
 *   a21 a22 ... a2n b2
 *   ...
 *   an1 an2 ... ann bn
 * @param filename Name of the file to read
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param n Pointer to the size of the system (number of equations)
 * @return true if the file was read successfully, false otherwise
 *  */
bool read_array_file(const char* filename, double a[][MAX_SIZE+1], double b[], int* n);


/* The matrix to solve is:
4x1 - x2 + 2x3 + 3x4 = 20
0x1 -2x2 + 7x3 - 4x4 = -7
0x1 + 0x2 + 6x3 + 5x4 = 4
0x1 + 0x2 + 0x3 + 3x4 = 6

In other terms: 
A = | 4 -1  2  3 |
    | 0 -2  7 -4 |
    | 0  0  6  5 |
    | 0  0  0  3 |
b = | 20 |
    | -7 |
    |  4 |
    |  6 |

We'll apply Gaussian elimination with partial pivoting to solve the system of equations.

The solution is:
x1 = 3.0
x2 = -4.0
x3 = -1.0
x4 = 2.0

Determinant of A = -144.0

*/

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, product, sum, aux;
    
    // Defining arrays using the global MAX_SIZE
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Read array from file using function
    if(!read_array_file("data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Original system of equations:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);       // 10 espacios totales, 6 decimales
    }
    printf("\n");
    // Walk the rows of the matrix (Gaussian elimination)
    for(int i = 1; i <= n-1; i++) {

        // ¿What happens if a[i][i] is near zero?
        // We'll use partial pivoting to avoid division by zero or numerical instability
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l; // We find the row with the largest element in column i
                }
            }
            for(int m = i; m <= n; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux; // Swap rows p and i
            }
            aux = b[p];
            b[p] = b[i];
            b[i] = aux; // Swap the independent term
        }

        // Makes zero the elements below the diagonal in the current column
        for(int j = i+1; j <= n; j++) {
            factor = a[j][i] / a[i][i]; // Sin el signo negativo
            
            // Traverses the columns in row j
            for(int k = i; k <= n; k++) {
                a[j][k] = a[j][k] - factor * a[i][k]; // Corregido: a[j][k] - factor * a[i][k]
            }
            b[j] = b[j] - factor * b[i]; // Corregido: b[j] - factor * b[i]
        }
    }

    // We print the matrix after Gaussian elimination, It's more convenient
    printf("The matrix after Gaussian elimination is:\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("The vector b after Gaussian elimination is:\n");
    for(int i = 1; i <= n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");


    // We verify the determinant of the matrix
    product = 1.0;
    for(int i = 1; i <= n; i++) {
        product = product * a[i][i];
    }

    printf("------------------DETERMINANT------------------\n");
    printf("The determinant of the matrix is: %lf\n\n", product);

    if(product == 0) {
        printf("The determinant of the matrix is zero, the system has no unique solution.\n");
        exit(0); // Exit with error code
    }

    // We perform back substitution to find the solution
    X[n] = b[n] / a[n][n];

    for(int i = n-1; i >= 1; i--) {
        sum = b[i];
        for(int j = i+1; j <= n; j++) {
            sum = sum - a[i][j] * X[j];
        }
        sum = sum / a[i][i];
        X[i] = sum;
    }
    printf("------------------SOLUTION------------------\n");
    printf("The solution of the system is:\n");
    for(int i = 1; i <= n; i++) {
        printf("X[%d] = %lf\n", i, X[i]);
    }
    

    return 0;
}

bool read_array_file(const char* filename, double a[][MAX_SIZE+1], double b[], int* n) {
    FILE *fp;
    char c;
    
    // Open data file
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("❌ Error: Cannot open file '%s'\n", filename);
        return false;
    }
    
    printf("✅ File '%s' opened\n\n", filename);

    // Count rows in file
    /* int rows = 0;
    while((c = fgetc(fp)) != EOF) {
        if(c == '\n') {
            rows++;
        }
    } */

    // This alternative works better than the previous one
    int rows = 0;
    while(!feof(fp)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), fp) != NULL) {
            rows++;
        }
    }
    
    // System must be square, so n = rows
    *n = rows;
    printf("📊 Size of the system: %d x %d\n", *n, *n);

    // Close and reopen the file to reset the pointer
    fclose(fp);
    fp = fopen(filename, "r");
    
    // Check maximum size using global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: System too big (%d). Maximum allowed: %d\n", *n, MAX_SIZE);
        fclose(fp);
        return false;
    }

    // Read augmented matrix from file
    // Expected format: each row contains n coefficients + 1 independent term
    // Example for 3x3: a11 a12 a13 b1
    //                   a21 a22 a23 b2  
    //                   a31 a32 a33 b3
    
    int i, j;
    for(i = 1; i <= *n; i++) {
        // Reading the matrix coefficients
        for(j = 1; j <= *n; j++) {
            if(fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("❌ Error reading element a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        // Read the term independent
        if(fscanf(fp, "%lf", &b[i]) != 1) {
            printf("❌ Error reading independent term b[%d]\n", i);
            fclose(fp);
            return false;
        }
    }
    
    fclose(fp);
    printf("✅ Array read successfully from file\n\n");
    
    return true;
}
