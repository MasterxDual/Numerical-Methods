#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100 
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


int main(int argc, char const *argv[]) {
    double X_hat, sum, product, error, fx, Pn;
    double X1[MAX_POINTS], Y1[MAX_POINTS];
    double X2[MAX_POINTS], Y2[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of interpolating polynomial
    double a[MAX_SIZE+1];
    int n, option;

    // Read data points from file
    if (!read_data_points("results_gauss.txt", X1, Y1, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    if (!read_data_points("results_seidel.txt", X2, Y2, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X1, Y1, n);
    print_data_points(X2, Y2, n);

    // Calculate error between both methods
    printf("Calculating errors between Gauss and Gauss-Seidel results:\n");
    printf("===============================================\n");
    printf("   i  |      X_Gauss      |   X_Seidel   |    Error    \n");
    printf("------|-------------------|-------------------|-------------\n");
    double numerator = 0, denominator = 0;
    for (int i = 0; i < n; i++) {
        numerator += fabs(X1[i] - X2[i]);
        denominator += fabs(X1[i]);
        error = numerator / denominator;
        printf("%4d  | %15.10f | %15.10f | %12.10f\n", i + 1, X1[i], X2[i], error);
    }
    printf("\n");
    

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