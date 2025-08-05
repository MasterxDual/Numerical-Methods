#include <iostream>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

/* Both methods are closed methods, meaning they require an interval [a, b] where the function changes sign.
Also, both methods require an initial guess for the root.
They use root location */

/* Watching convergence_analysis.png
Left graph: Global convergence (logarithmic scale)
    The False Rule converges faster than Bisection. (It reaches the desired error (1e-5) in just 7 iterations.)
    Bisection requires about 14 iterations to reach the same level of accuracy.
    The steeper slope of the red curve indicates faster convergence.
Right graph: First 7 iterations (linear scale)
    This graph allows a closer look at the initial behavior of both methods.
    From the second iteration on, the False Rule reduces the error much more efficiently than Bisection.
    By the fourth iteration, the False Rule has almost reached the desired accuracy.
In problems where speed is required and the function is continuous and well-behaved, the False Rule may be preferred over Bisection.
*/

// Function prototypes
double math_function(double x);
int get_error_type();
double calculate_error(double c, double old_c, int error_type);
void update_interval(double *a, double *b, double c);
void print_results(double c, double error, int max_iterations);
void save_convergence_data(const char* filename, double* errors, int iterations);

/* 
To find 'a' and 'b', we used the graphic method, which is a common practice in numerical methods.
Bisection Method applied
Root 1:
    a = -1.05
    b = -0.9
    tolerance = 1e-5
    --> Aproximated root = -1.000003, absolute error = 0.000009 <--> -1.000003 +- 0.000009
    --> Number of iterations = 14
Root 2:
    a = 0.9
    b = 1.2
    tolerance = 1e-5
    --> Aproximated root = 1.000003, absolute error = 0.000009 <--> 1.000003 +- 0.000009
    --> Number of iterations = 15

Regula Falsi Method applied
Root 1:
    a = -1.05
    b = -0.9
    tolerance = 1e-5
    --> Aproximated root = -0.999998, absolute error = 0.000007 <--> -0.999998 +- 0.000007
    --> Number of iterations = 7
Root 2:
    a = 0.9
    b = 1.2
    tolerance = 1e-5
    --> Aproximated root = 0.999988, absolute error = 0.000008 <--> 0.999988 +- 0.000008
    --> Number of iterations = 19

Comparing both methods, the Regula Falsi method took fewer iterations than the bisection method, demonstrating that it is more effective (only in cases where
the functions analyzed do not grow rapidly near their roots).
In general, the false rule method tends to converge faster than the bisection method.
This is because the false rule uses linear interpolation to estimate the root at each iteration,
which can provide more accurate approximations and allow for faster convergence for rapidly changing functions.
At the second root, the Regula Falsi method requires a greater number of iterations than the bisection method because at that root, the function grows rapidly, decreasing the convergence of the Regula Falsi method.
*/
int main(int argc, char const *argv[]) {
    // Variable definitions
    double a, b, tolerance, old_c, c, error;
    int max_iterations = 0;
    int method, error_type;
    
    // Array to store convergence data (max 1000 iterations)
    double errors_history[1000];
    
    // Request user to input values
    printf("Ingrese el valor de a: ");
    scanf("%lf", &a);
    
    printf("Ingrese el valor de b: ");
    scanf("%lf", &b);
    
    printf("Ingrese la tolerancia del error: ");
    scanf("%lf", &tolerance);
    
    // Check if the function has a root in the interval [a, b]
    if(math_function(a) * math_function(b) > 0) {
        printf("No se puede garantizar la existencia de una raíz en el intervalo [%lf, %lf]\n", a, b);
        exit(0);
    }
    old_c = a; // old_c will be equal to a at minus one iteration

    printf("Ingrese el tipo de metodo (1 para bisección, 2 para regula falsi): ");
    scanf("%d", &method);
    
    switch(method) {
        case 1:
            printf("Método de Bisección seleccionado\n");
            error_type = get_error_type();
            do {
                c = (a + b) / 2; // We divide the interval in half
                max_iterations++; // We count the number of iterations
                
                update_interval(&a, &b, c);
                error = calculate_error(c, old_c, error_type);
                
                // Store error for convergence analysis
                errors_history[max_iterations - 1] = error;
                
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            
            // Save convergence data to file
            save_convergence_data("bisection_convergence.txt", errors_history, max_iterations);
            printf("\n✅ Datos de convergencia guardados en 'bisection_convergence.txt'\n");
            break;
        case 2:
            printf("Método de Regula Falsi seleccionado\n");
            error_type = get_error_type();
            do {
                // This formula is derived from the equation of a line
                c = ((a * math_function(b)) - (b * math_function(a))) / (math_function(b) - math_function(a)); // We calculate the point c using the Regula Falsi formula
                max_iterations++; // We count the number of iterations
                
                update_interval(&a, &b, c);
                error = calculate_error(c, old_c, error_type);
                
                // Store error for convergence analysis
                errors_history[max_iterations - 1] = error;
                
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            
            // Save convergence data to file
            save_convergence_data("regula_falsi_convergence.txt", errors_history, max_iterations);
            printf("\n✅ Datos de convergencia guardados en 'regula_falsi_convergence.txt'\n");
            break;
        default:
            printf("Opción inválida. Seleccione 1 o 2.\n");
    }
    
    return 0;
}

// Mathematical function to find its root
double math_function(double x) {
    return (pow(x, 10) - 1); //This is the function from exercise 2
}

// Function to get error type from user
int get_error_type() {
    int error_type;
    printf("¿Solicita error absoluto o porcentual? (1 para absoluto, 2 para porcentual): ");
    scanf("%d", &error_type);
    return error_type;
}

// Function to calculate error (absolute or percentage)
double calculate_error(double c, double old_c, int error_type) {
    if(error_type == 2) {
        return fabs((c - old_c) / c) * 100; // Percentage error
    } else if(error_type == 1) {
        return fabs(c - old_c); // Absolute error
    } else {
        printf("Tipo de error inválido. Seleccione 1 o 2.\n");
        exit(0);
    }
}

// Function to update interval based on sign analysis
void update_interval(double *a, double *b, double c) {
    if(math_function(*a) * math_function(c) > 0) {
        *a = c; // The root is in the interval [c, b]
    } else if(math_function(*a) * math_function(c) < 0) {
        *b = c; // The root is in the interval [a, c]
    } else {
        printf("La raíz exacta es (sin error): %lf\n", c);
        exit(1);
    }
}

// Function to print final results
void print_results(double c, double error, int max_iterations) {
    printf("La raíz aproximada es: %lf, con un error de %lf\n", c, error);
    printf("Número de iteraciones: %d\n", max_iterations);
}

// Function to save convergence data to file for plotting
void save_convergence_data(const char* filename, double* errors, int iterations) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("❌ Error: No se pudo crear el archivo %s\n", filename);
        return;
    }
    
    // Write header with description
    fprintf(file, "# Datos de convergencia - f(x) = x^10 - 1\n");
    fprintf(file, "# Columna 1: Iteración\n");
    fprintf(file, "# Columna 2: Error absoluto\n");
    fprintf(file, "Iteracion\tError\n");
    
    // Write data
    for (int i = 0; i < iterations; i++) {
        fprintf(file, "%d\t\t\t%.12e\n", i + 1, errors[i]);
    }
    
    fclose(file);
    
    // Print summary to console
    printf("📊 Archivo creado: %s\n", filename);
    printf("   - Total de iteraciones: %d\n", iterations);
    printf("   - Error inicial: %.6e\n", errors[0]);
    printf("   - Error final: %.6e\n", errors[iterations-1]);
}
