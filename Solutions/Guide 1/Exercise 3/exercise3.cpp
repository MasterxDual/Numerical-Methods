#include <iostream>
#include <cmath>

using namespace std;

/* Both methods are closed methods, meaning they require an interval [a, b] where the function changes sign.
Also, both methods require an initial guess for the root.
They use root location */

// Function prototypes
double math_function(double x);
int get_error_type();
double calculate_error(double c, double old_c, int error_type);
void update_interval(double *a, double *b, double c);
void print_results(double c, double error, int max_iterations);

/* 
Using basic maths to solve the equation we finally get the following function:
    f(m) = m - (m * exp(-98 / m)) - 49.95 = 0
    We use this function in math_function(m) to find its root.
Bisection Method applied
    a = 61
    b = 64
    tolerance = 1e-3
    --> Aproximated root = 63.539307, absolute error = 0.000732 <--> 63.5393 +- 0.0007
    --> Number of iterations = 12

Regula Falsi applied
    a = 61
    b = 64
    tolerance = 1e-3
    --> Aproximated root = 63.540001, absolute error = 0.000222 <--> 63.5400 +- 0.0002
    --> Number of iterations = 3

Now we use one or another root to confirm if the root is correct.
We did this writing in a paper.
If the equation gives us aproximately 0, then the root is correct.
*/
int main(int argc, char const *argv[]) {
    // Variable definitions
    double a, b, tolerance, old_c, c, error;
    int max_iterations = 0;
    int method, error_type;
    
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
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            break;
        case 2:
            printf("Método de Regula Falsi seleccionado\n");
            error_type = get_error_type();
            do {
                // This formula is derived from the equation of a line
                c = (a * math_function(b) - b * math_function(a)) / (math_function(b) - math_function(a)); // We calculate the point c using the Regula Falsi formula
                max_iterations++; // We count the number of iterations
                
                update_interval(&a, &b, c);
                error = calculate_error(c, old_c, error_type);
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            break;
        default:
            printf("Opción inválida. Seleccione 1 o 2.\n");
    }
    
    return 0;
}

// Mathematical function to find its root
double math_function(double m) {
    return (m - (m * exp(-98 / m)) - 49.95); //This is the function from exercise 3
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


