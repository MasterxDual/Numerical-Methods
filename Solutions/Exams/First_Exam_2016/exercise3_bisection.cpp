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
void print_results(double c, double error, int max_iterations, int error_type);


/*  
Method used: Bisection Method
    Interval is: [0.8, 1.3] 
    Tolerance = 1e-8
    Root = 0.83255460
    Absolute error = 0.00000001 
    Porcentual error = 0.00000089%
    Iterations number = 26
    With Lagrange Polynomial:
    If we use the same interval we can see that there's no any root. So we need to use another extended interval. 
    In the interval [0.5, 2.0]:
    Root = 0.61527059 
    Absolute error = 0.00000001
    Percentual error = 0.00000091%
    Iterations number = 28
 */

int main(int argc, char const *argv[]) {
    // Variable definitions
    double a, b, tolerance, old_c, c, error;
    int max_iterations = 0;
    int method;
    int error_type;
    
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

            print_results(c, error, max_iterations, error_type);
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

            print_results(c, error, max_iterations, error_type);
            break;
        default:
            printf("Opción inválida. Seleccione 1 o 2.\n");
    }
    
    return 0;
}

// Mathematical function to find its root
double math_function(double x) {
    return (-15.33 + 12.6089*x + 28.0355*pow(x, 2) - 13.056*pow(x, 3));
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
        return fabs((c - old_c) / c); // Percentage error
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
void print_results(double c, double error, int max_iterations, int error_type) {
    if(error_type == 2) {
        printf("La raíz aproximada es: %0.8lf, con un error de %0.8lf%%\n", c, error*100);
        printf("Número de iteraciones: %d\n", max_iterations);
    } else {
        printf("La raíz aproximada es: %0.8lf, con un error de %0.8lf\n", c, error);
        printf("Número de iteraciones: %d\n", max_iterations);
    }
}


