#include <cstdio>
#include <stdlib.h>

// Now you can handle up to 50x50 matrixs
#define MAX_SIZE 50  

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
    int rows = 0;
    while((c = fgetc(fp)) != EOF) {
        if(c == '\n') {
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

int main(int argc, char const *argv[]) {
    int n;
    double factor, product, sum;
    
    // Defining arrays using the global MAX_SIZE
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Read array from file using function
    if(!read_array_file("datos.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Original system of equations:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%6.2lf ", a[i][j]);
        }
        printf("| %6.2lf\n", b[i]);
    }
    printf("\n");
    // Walk the rows of the matrix (Eliminación Gaussiana)
    for(int i = 1; i <= n-1; i++) {
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
    printf("The solution of the system is:\n");
    for(int i = 1; i <= n; i++) {
        printf("X[%d] = %lf\n", i, X[i]);
    }
    

    return 0;
}
