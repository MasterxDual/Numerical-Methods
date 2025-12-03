#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100
#define PI 3.14159265358979323846

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
 * Function to define function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double f(double x);

/**
 * Function to save computed values to a text file
 * @param x Array of x values
 * @param fp Array of first derivative values
 * @param n Number of subintervals
 */
void save_in_txt(double x[], double fp[], int n);

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x);

int main(int argc, char const *argv[]) {
    // Interval extremes
    double a, b;    
    // Number of subintervals
    int n;
    // Distance between points
    double h;
    // First derivative value
    double fp[MAX_SIZE + 1];
    // Second derivative value
    double fpp[MAX_SIZE + 1];
    // Third derivative value
    double fppp[MAX_SIZE + 1];
    // Domain values
    double x[MAX_SIZE + 1];
    // Error type O(h) or O(h²). Remember that centered finite difference operators have error O(h²) for O(h) and have error O(h⁴) for O(h²)
    int error_type;
    int number_of_derivative;
    // Variable to store if it's a function or a data table
    int function_data_table;
    
    printf("Insert the derivative that you want to calculate: 1.First derivative 2.Second derivative 3.Third derivative\n");
    scanf("%d", &number_of_derivative);
    
    printf("Insert the error that you want to use: 1. O(h) 2. O(h²)\n");
    scanf("%d", &error_type);

    printf("Do you have a function or Do you have a data table?\n");
    printf("1. I have a function\n");
    printf("2. I have a data table\n");
    scanf("%d", &function_data_table);

    if(error_type == 2) {
        switch(number_of_derivative) {
            case 1: 
                if(function_data_table == 2) {
                    // Arrays for polynomial coefficients calculation to use for Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Data points to read from text file
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Data points to calculate the integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Number of points
                    int n;

                    // Read data points from file
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("Failed to read data from file. Exiting.\n");
                        return 1;
                    }
                    // Print the data points
                    print_data_points(X, Y, n);

                    // Use code of above if original data is not uniformed spaced
                    h = X[1] - X[0]; // Assuming uniform spacing in the original data
                    
                    // Calculate first derivative with error O(h²)
                    // Forward difference O(h²) for first point (needs 3 points)
                    fp[0] = (-Y[2] + 4*Y[1] - 3*Y[0])/(2*h);
                    
                    // Backward difference O(h²) for last point (needs 3 points)
                    fp[n-1] = (3*Y[n-1] - 4*Y[n-2] + Y[n-3])/(2*h);
                
                    // Calculate derivatives for all inner points with appropriate formulas
                    /* for(int i = 1; i <= n-2; i++) {
                        // if(i >= 2 && i <= n-3) {
                            // Central difference O(h⁴) when we have 5 points available
                        fp[i] = (-Y[i+2] + 8*Y[i+1] - 8*Y[i-1] + Y[i-2])/(12*h);
                        // } else {
                            // Central difference O(h²) when we only have 3 points
                            // fp[i] = (Y[i+1] - Y[i-1])/(2*h);
                        // }
                    } */
                    
                    for(int i = 1; i <= n-2; i++) {
                        // x[i] = a + i*h;
                        // fp[i] = (f(x[i]+h) - f(x[i]-h)) / (2*h); // This is O(h²)
                        fp[i] = (Y[i+1] - Y[i-1]) / (2*h); // This is O(h²)
                    }

                    // Print results
                    printf("x\t\tf'(x)\n");
                    for(int i = 0; i < n; i++) {
                        printf("%lf\t%lf\n", X[i], fp[i]);
                    }
                
                    // Save x[i] and fp[i] in a text file
                    save_in_txt(X, fp, n-1);
                
                    // Finally, we print the results.txt file in a graph using Python to visualize the results
                    system("python3 graph_points.py");
                    break;
                }
        }
    }

    return 0;
}

double f(double t) {
    return 10 * exp(-t/10) * sin(2*t);
}

void save_in_txt(double x[], double fp[], int n) {
    FILE *archivo = fopen("results.txt", "w");
    if (archivo == NULL) {
        printf("Error: Unable to create file.\n");
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", x[i], fp[i]);
    }

    fclose(archivo);
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

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Find the correct spline interval
    if (x <= X[0]) {
        k = 0;  // Use first spline for x smaller than first point
    } else if (x >= X[n-1]) {
        k = n - 2; // Use last spline for x larger than last point
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluate S_k(x) = a*x^3 + b*x^2 + c*x + d
    double y = solution[4*k] * pow(x, 3) +
               solution[4*k+1] * pow(x, 2) +
               solution[4*k+2] * x +
               solution[4*k+3];

    return y;
}

double X(double tita) {
    // x(tita) = R*cos(tita) + (L² - R² * sin²(tita))^(1/2)
    // Where R = 90 mm = 0.09 m and L = 2.5*R = 0.225 m
    double R = 0.09;
    double L = 0.225;
    return (R * cos(tita)) + (sqrt(pow(L, 2) - pow(R, 2) * pow(sin(tita), 2)));
}

void points_generator_deg(const char *file_name) {
    FILE *fp = fopen(file_name, "w");
    if (!fp) {
        perror("Error al abrir el archivo");
        return;
    }

    double tita_deg, tita_rad, x;
    // double h = 2 * PI / (N - 1);

    for (int i = 0; i <= 36; i++) {
        tita_deg = i * 5; // From 0 to 180 degrees in steps of 5 degrees
        tita_rad = tita_deg * (PI / 180.0); // Convert to radians

        x = X(tita_rad);
        fprintf(fp, "%.6f %.6f\n", tita_deg, x);
    }

    fclose(fp);
    printf("Archivo '%s' generado correctamente.\n", file_name);
}

void points_generator_rad(const char *file_name) {
    FILE *fp = fopen(file_name, "w");
    if (!fp) {
        perror("Error al abrir el archivo");
        return;
    }

    double tita_deg, tita_rad, x;
    // double h = 2 * PI / (N - 1);

    for (int i = 0; i <= 36; i++) {
        tita_deg = i * 5; // From 0 to 180 degrees in steps of 5 degrees
        tita_rad = tita_deg * (PI / 180.0); // Convert to radians

        x = X(tita_rad);
        fprintf(fp, "%.6f %.6f\n", tita_rad, x);
    }

    fclose(fp);
    printf("Archivo '%s' generado correctamente.\n", file_name);
}