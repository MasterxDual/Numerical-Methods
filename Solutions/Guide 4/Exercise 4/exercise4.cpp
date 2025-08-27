#include <cstdio>
#include <stdlib.h>
#include <cmath>

/* This file have bug fixes done with Github Copilot in Gauss-Siedel and Relaxation Method. */

// Now you can handle up to 50x50 matrixs
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

/**
 * Checks if the matrix is ​​diagonally dominant and that there are no zeros on the diagonal.
 * Returns 0 if it's everything ok, otherwise returns 1
 * @param a Matrix of coefficients
 * @param n Size of the matrix
 * @return 0 if everything is ok, 1 if there is a zero on the diagonal
 *         or a warning if the matrix is not diagonally dominant
 */
int checkDiagonalDominance(double a[][MAX_SIZE+1], int n);

/** Function to initialize the guess vector with zeros
 * @param Xv Vector to initialize == Old X
 * @param n Size of the vector
 * 
 *  */ 
void initializeGuess(double Xv[], int n);

/** Function to compute the error between Xn and Xv
 * @param Xn New X = Solution of the matrix
 * @param Xv Old X = Previous iteration of the solution
 * @param n Size of the vectors
 * @return The computed error (Euclidean norm)
 *  */
double computeError(double Xn[], double Xv[], int n);


/** Function to print the solution
 * @param methodName Name of the method used
 * @param Xn Solution vector
 * @param n Size of the vector
 * @param iterations Number of iterations taken to converge
 * @param error Final error
 *  */ 
void printSolution(const char* methodName, double Xn[], int n, int iterations, double error);


/**
 * Implementation of the Jacobi method for solving linear systems
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void jacobiMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/**
 * Implementation of the Jacobi method for solving linear systems
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void gaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);


/**
 * Implementation of the Gauss-Seidel method optimized for banded matrices
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void gaussSeidelBand(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/** Function to calculate the bandwidth of a matrix
 * The bandwidth is defined as the width of the band around the main diagonal
 * that contains all the non-zero elements of the matrix.
 * @param a Matrix of coefficients
 * @param n Size of the matrix
 * @return The bandwidth of the matrix
 *  */
int computeBandwidth(double a[][MAX_SIZE+1], int n);

/**
 * Implementation of the Relaxation method for solving linear systems
 * It's similar to Gauss-Seidel but with a relaxation factor omega and one additional line of code
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void relaxationMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);
/* 
The matrix to solve is larger, take a look to data.dat file.
We use Gauss-Seidel without bandwidth optimization.
    Tolerance = 1e-11

    The solution is:
        x1 = 0.46
        x2 = 0.53
        x3 = 0.51
        x4 = 0.50
        x5 = 0.50
        x6 = 0.50
        x7 = 0.50
        ...
        x49 = 0.53
        x50 = 0.46

    
    The method converged in 16 iterations with an error of 0.0
    
    We use Gauss-Seidel with bandwidth optimization.
    The bandwidth of the matrix is 2

    The solution is:
        x1 = 0.46
        x2 = 0.53
        x3 = 0.51
        x4 = 0.50
        x5 = 0.50
        x6 = 0.50
        x7 = 0.50
        ...
        x48 = 0.51
        x49 = 0.53
        x50 = 0.46
    The method converged in 16 iterations with an error of 0.0
    
    On larger matrices or with smaller bandwidth, the optimization:
    --> Reduce operations: Avoid unnecessary multiplications by zero
    --> Improved speed: Fewer iterations in the inner loop
    --> In your case: Since bandwidth (2) is very small compared to n (50), both methods do practically the same job.
    */



int main(int argc, char const *argv[]) {
    int n, p;
    double factor, product, sum, aux;
    double Xv[MAX_SIZE+1], Xn[MAX_SIZE+1]; // Old X, New X
    double tolerance, old_error, new_error;
    int iterations;
    double omega;

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
            printf("%10.6lf ", a[i][j]);  // 10 total spaces, 6 decimal places
        }
        printf("| %10.6lf\n", b[i]);       // 10 total spaces, 6 decimal places
    }
    printf("\n");

    // Verificación de diagonal dominante
    if (checkDiagonalDominance(a, n) != 0) {
        printf("The method cannot continue. The matrix has zeros on the diagonal. The program exits.");
        return 1;
    }
    printf("Verificación completada.\n");

    printf("Choose a method to solve the system:\n");
    printf("1. Jacobi Method\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Relaxation Method\n");
    printf("4. Gauss-Seidel with Bandwidth Optimization\n");
    printf("Enter option: ");

    int option;
    scanf("%d", &option);

    switch(option) {
        case 1:
            jacobiMethod(a, b, Xv, Xn, n);
            break;
        case 2:
            gaussSeidel(a, b, Xv, Xn, n);
            break;
        case 3:
            relaxationMethod(a, b, Xv, Xn, n);
            break;
        case 4:
            gaussSeidelBand(a, b, Xv, Xn, n);
            break;
        default:
            printf("❌ Invalid option.\n");
    }


    return 0;
}

int checkDiagonalDominance(double a[][MAX_SIZE+1], int n) {
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;

        // We first check if there is zero on the diagonal
        if (fabs(a[i][i]) == 0.0) {
            printf("❌ Error: Zero element on the diagonal at position a[%d][%d].\n", i, i);
            return 1;
        }

        // Sum of off-diagonal elements (elementos fuera de la diagonal)
        for (int j = 1; j <= n; j++) {
            if (j != i) {
                sum += fabs(a[i][j]);
            }
        }

        // We check the dominance condition
        if (fabs(a[i][i]) < sum) {
            printf("⚠️  Warning: The matrix is not diagonally dominant at row %d.\n", i);
        }
    }

    return 0; // Everything is OK
}

void initializeGuess(double Xv[], int n) {
    for(int i = 1; i <= n; i++) {
        Xv[i] = 0.0;
    }
}

double computeError(double Xn[], double Xv[], int n) {
    double error = 0.0;
    for(int i = 1; i <= n; i++) {
        error += pow(Xn[i] - Xv[i], 2);
    }
    return sqrt(error);
}

void printSolution(const char* methodName, double Xn[], int n, int iterations, double error) {
    printf("------------------SOLUTION OF %s------------------\n", methodName);
    printf("The solution of the system is:\n");
    for(int i = 1; i <= n; i++) {
        printf("Xn[%d] = %10.6lf\n", i, Xn[i]);
    }
    printf("The method converged in %d iterations with an error of %10.6lf\n", iterations, error);
}


void jacobiMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;
            for(int j = 1; j <= n; j++) {
                if(j != i) {
                    sum += a[i][j] * Xv[j];
                }
            }
            Xn[i] = (b[i] - sum) / a[i][i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("JACOBI METHOD", Xn, n, iterations, new_error);
}

void gaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;  // Reset sum for each row
            
            // Sum elements before diagonal (using NEW values Xn)
            for(int j = 1; j <= i-1; j++) {
                sum += a[i][j] * Xn[j];
            }
            
            // Sum elements after diagonal (using OLD values Xv)
            for(int j = i+1; j <= n; j++) {
                sum += a[i][j] * Xv[j];
            }
            
            Xn[i] = (b[i] - sum) / a[i][i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("GAUSS-SEIDEL", Xn, n, iterations, new_error);
}

int computeBandwidth(double a[][MAX_SIZE+1], int n) {
    int bw = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (fabs(a[i][j]) > 1e-12) { // There's a coefficient different from zero
                int dist = abs(i - j);
                if (dist > bw) {
                    bw = dist;
                }
            }
        }
    }
    return bw;
}

void gaussSeidelBand(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    int bw = computeBandwidth(a, n);
    printf("📏 Bandwidth of matrix = %d\n", bw);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for (int i = 1; i <= n; i++) {
            sum = 0.0;

            // Loop through only columns within the band
            int jmin = (i - bw > 1) ? i - bw : 1;
            int jmax = (i + bw < n) ? i + bw : n;

            for (int j = jmin; j <= jmax; j++) {
                if (j != i) {
                    if (j < i) sum += a[i][j] * Xn[j]; // Updated now
                    else       sum += a[i][j] * Xv[j]; // Still old
                }
            }

            Xn[i] = (b[i] - sum) / a[i][i];
        }

        new_error = computeError(Xn, Xv, n);

        if (new_error > old_error) {
            printf("❌ The method does not converge, stopping process.\n");
            return;
        }

        old_error = new_error;

        for (int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while (new_error > tolerance);

    printSolution("GAUSS-SEIDEL with BAND", Xn, n, iterations, new_error);
}


void relaxationMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error, omega;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    printf("Please enter relaxation factor (0 < omega < 2):");
    scanf("%lf", &omega);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;  // Reset sum for each row
            
            // Sum elements before diagonal (using NEW values Xn)
            for(int j = 1; j <= i-1; j++) {
                sum += a[i][j] * Xn[j];
            }
            
            // Sum elements after diagonal (using OLD values Xv)
            for(int j = i+1; j <= n; j++) {
                sum += a[i][j] * Xv[j];
            }
            
            // Calculate Gauss-Seidel step
            double gauss_seidel = (b[i] - sum) / a[i][i];
            
            // Apply relaxation factor (SOR)
            Xn[i] = omega * gauss_seidel + (1.0 - omega) * Xv[i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("RELAXATION METHOD", Xn, n, iterations, new_error);
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
