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

/**
 * Function to evaluate the actual function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double valuate_function(double x);

/* In this exercise we use the implementation done in exercise 5, so we have the data points of data.txt
    Then, we use the interpolation of grade 2 and we have the result:
    f(h) = A*h² + B*h + C ==> f(h) = 0.998895*h² - 0.093447*h + 0.002763
    The exercise requires to calculate relative density of air at h = 10.5 km
    So we evaluate f(10.5) and we have the result f(10.5) = 0.322342
    r = 0.999988 --> Excellent fit
    Mean Square Error = 0.001018
*/

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
    double sumxy, sumx, sumy, distance, mean_square_error;
    double Sr, St, r, average_y, f_xi, relative_density_of_air;

    printf("Enter the degree of the polynomial (degree >= 1): ");
    scanf("%d", &degree);

    // Read data points from file
    if (!read_data_points("data.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    // Verify if we have enough data points
    if(n < degree) {
        printf("Error: Not enough data points for the chosen polynomial degree.\n");
        return 1;
    }

    // Construct the matrix A and vector b for the normal equations
    // We use 0,...,n-1 for the points, but in the calculus we use 1,....N

    for(int l = 4; l >= 2; l--) {
        sumxy = 0.0;
        for(int k = 0; k < n; k++) {
            sumxy += pow(X[k], l-2) * Y[k];
        }
        b[4-l] = sumxy;
        for(int m = 0; m <= 2; m++) {
            sumx = 0.0;
            for(int k = 0; k < n; k++) {
                sumx += pow(X[k], l-m);
            }
            if(l == 2 && m == 2) {
                A[4-l][m] = n;
                break;
            }
            A[4-l][m] = sumx;
        }
    }

    // We use the function from gauss.h to solve the system with Gaussian elimination
    gauss_elimination(degree + 1, A, b, solution);

    // We copy the solution to a[i] to give relevance to our context
    for(int i = 2; i >= 0; i--) {
        a[2-i] = solution[i];
    }

    // Print the polynomial coefficients
    printf("------------------SOLUTION------------------\n");
    printf("The solution of the system is:\n");
    for(int i = 0; i <= degree; i++) {
        printf("a[%d] = %lf\n", i, a[i]);
    }
    
    printf("\n------------------INTERPOLATING POLYNOMIAL------------------\n");
    print_polynomial(a, degree + 1);

    // Now we evaluate the polynomial at h = 10.5 km to find the relative density of air
    relative_density_of_air = 0.0;
    for(int i = 0; i <= degree; i++) {
        relative_density_of_air += a[i] * pow(10.5, i);
    }
    printf("The relative density of air at h = 10.5 km is: %lf\n", relative_density_of_air);

    // Calculate mean square error
    for(int i = 0; i < n; i++) {
        f_xi = 0.0;
        for(int k = 0; k <= degree; k++) {
            f_xi += a[k] * pow(X[i], k);
        }
        mean_square_error += pow((f_xi - Y[i]), 2);
    }
    mean_square_error /= n;
    mean_square_error = sqrt(mean_square_error);
    printf("The mean square error is: %lf\n", mean_square_error);

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

    // Calculate Sr - sum of squared residuals
    for(int i = 0; i < n; i++) {
        f_xi = 0.0;
        for(int k = 0; k <= degree; k++) {
            f_xi += a[k] * pow(X[i], k);
        }
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
    printf("Pn(x) = ");
    
    // Handle the first term (constant term a0)
    if (fabs(a[0]) > 1e-10) {  // Avoid printing very small values as zero
        printf("%.6f", a[0]);
    } else {
        printf("0");
    }
    
    // Handle the rest of the terms
    for (int i = 1; i < n; i++) {
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
    }
    printf("\n\n");
}

double valuate_function(double x) {
    // f(x) = e^(-x^2)
    return exp(- (x*x));
}