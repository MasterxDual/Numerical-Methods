#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define MAX_SIZE 100 
#include "gauss.h"

/**
 * Function to read Xi, Yi data pairs from a file
 * Expected format: first line contains number of points, then xi yi pairs
 * @param filename Name of the file to read
 * @param X Array to store X values
 * @param Y Array to store Y values
 * @param n Pointer to store the number of data points read
 * @return 1 if successful, 0 otherwise
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * Function to display the data points
 * @param X Array of X values
 * @param Y Array of Y values
 * @param n Number of data points
 */
void print_data_points(double X[], double Y[], int n);


/**
 * Function to print the interpolating polynomial Pn(x)
 * @param a Array of coefficients
 * @param n Degree of polynomial (n-1 is the highest power)
 */
void print_polynomial(double a[], int n);

int main(int argc, char const *argv[]) {
    // Arrays for data points to read from text file
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of Least Squares polynomial
    double a[MAX_SIZE+1];
    // n = number of data points and degree is the polynomial degree
    int n, degree;
    // Sums for least squares calculations
    double sumxy, sumx, sumy, distance;
    double Sr, St, r, average_y, f_xi;

    /* printf("Enter the degree of the polynomial (degree >= 1): ");
    scanf("%d", &degree); */

    // Read data points from file
    if (!read_data_points("data.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    // Verify if we have enough data points
    /* if(n < degree) {
        printf("Error: Not enough data points for the chosen polynomial degree.\n");
        return 1;
    } */

    // Se construye la matriz A y el vector b segun la funcion planteada en el ejercicio 3 que depende de un seno y un coseno
    double sum_ysincos;
    double sum_sincos;
    for(int l = 0; l <= 1; l++) {
        sum_ysincos= 0.0;
        for(int k = 0; k < n; k++) {
            sum_ysincos += (Y[k] * pow(sin(X[k]), 1-l) * pow(cos(X[k]), l));
        }
        b[l] = sum_ysincos;
        for(int m = 0; m <= 1; m++) {
            sum_sincos = 0.0;
            for(int k = 0; k < n; k++) {
                sum_sincos += (pow(sin(X[k]), 2-l-m) * pow(cos(X[k]), l+m));
            }
            A[l][m] = sum_sincos;
        }
    }

    // We use the function from gauss.h to solve the system with Gaussian elimination
    gauss_elimination(2, A, b, solution);

    // We copy the solution to a[i] to give relevance to our context
    for(int i = 0; i <= 1; i++) {
        a[i] = solution[i];
    }

    // Print the polynomial coefficients
    printf("------------------SOLUTION------------------\n");
    printf("The solution of the system is:\n");
    for(int i = 0; i <= 1; i++) {
        printf("a[%d] = %lf\n", i, a[i]);
    }
    
    printf("\n------------------INTERPOLATING POLYNOMIAL------------------\n");
    // print_polynomial(a, n);
    printf("The polynomial is: f(x) = %lf * sin(x) + %lf * cos(x)\n", a[0], a[1]);


    // Calculate Sr, St, r 
    sumy = 0.0;
    St = 0.0;
    Sr = 0.0;

    // Calculate average of Y
    for(int i = 0; i < n; i++) {
        sumy += Y[i];
    }
    average_y = sumy / n;

    // Calculate St - total sum of squares
    for(int i = 0; i < n; i++) {
        St += pow((average_y - Y[i]), 2);
    }
    
    f_xi = 0.0;

    // Calculate Sr - sum of squared residuals
    for(int i = 0; i < n; i++) {
        /* for(int k = 0; k <= degree; k++) {
            f_xi += a[k] * pow(X[i], k);
        } */
        f_xi = a[0] * sin(X[i]) + a[1] * cos(X[i]);
        Sr += pow((f_xi - Y[i]), 2);
    }
    
    // Calculate r - correlation coefficient
    r = sqrt((St - Sr) / St);
    
    // Print correlation coefficient 
    printf("\nThe correlation coefficient r is: %lf\n", r);

    // Verify goodness of fit
    distance = fabs(r - 1.0);

    // Interpretation of goodness of fit
    if(distance < 0.1) {
        printf("The fit is very good (r is close to 1)\n");
    } else if(distance < 0.25) {
        printf("The fit is good\n");
    } else if(distance < 0.5) {
        printf("The fit is acceptable\n");
    } else {
        printf("The fit is poor\n");
    }


    return 0;
}

int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    printf("File '%s' opened successfully\n", filename);
    
    // Read number of data points
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: Cannot read number of data points\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Invalid number of points (%d)\n", *n);
        fclose(fp);
        return 0;
    }
    
    // Read data points
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: Cannot read data point %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Successfully read %d data points\n\n", *n);
    return 1;
}


void print_data_points(double X[], double Y[], int n) {
    printf("Data Points:\n");
    printf("=============\n");
    printf("   i  |      Xi      |      Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

void print_polynomial(double a[], int n) {
    /* printf("Pn(x) = ");
    
    // Handle the first term (constant term a0)
    if (fabs(a[0]) > 1e-10) {  // Avoid printing very small values as zero
        printf("%.6f", a[0]);
    } else {
        printf("0");
    } */
    
    // Handle the rest of the terms
    /* for (int i = 1; i < n; i++) {
        if (fabs(a[i]) > 1e-10) {  // Only print if coefficient is significant
            // Print sign
            if (a[i] > 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
            
            // Print coefficient (absolute value since sign is already printed)
            double coeff = fabs(a[i]);
            if (coeff != 1.0) {
                printf("%.6f", coeff);
            }
            
            // Print variable part
            if (i == 1) {
                printf("x");
            } else {
                printf("x^%d", i);
            }
        }
    } */
    printf("\n\n");
}