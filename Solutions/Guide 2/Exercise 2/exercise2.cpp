#include <iostream>
#include <cmath>

using namespace std;

/*
CABE ACLARAR, QUE SE SUPONE QUE LAS FUNCIONES QUE ESTAN EN EL ENUNCIADO SON g(x), NO HAY QUE HACER NINGUN DESPEJE, SINO QUE SE USAN TAL CUAL.

g(x) = pow(x, 5) - 3*pow(x, 3) - 2*pow(x, 2) + 2, f(x) = pow(x, 5) - 3*pow(x, 3) - 2*pow(x, 2) + 2
    Viendo la gráfica, vemos que no converge en ningún punto fijo
    Nunca converge. Su derivada es mayor a 1.

g(a) = cos(sin(a)), f(a) = cos(sin(a))
    El punto fijo de g(a) o la raiz de f(a) es igual a 0.768189133355
    La raiz aproximada es: 0.768169156737, con un error de 0.000000000001
    Número de iteraciones: 35
    La funcion evaluada en la raiz aproximada es igual a 0.768169156737

g(n) = pow(n, n-cos(n)), f(n) = pow(n, n-cos(n))
    Viendo la gráfica, vemos que converge únicamente en la primer raíz.
    El punto fijo de g(n) o la raiz de f(n) es igual a 0.999999999999
    La raiz aproximada es: 0.999999999999, con un error de 0.000000000001
    Número de iteraciones: 32
    La funcion evaluada en la raiz aproximada es igual a 1.000000000000
*/

// Function prototypes
double g(double x);
double gprima(double x);
double calculate_error(double new_x, double previous_x, int error_type);
void show_results(const char* method, double root, double error, int iterations);
void verify_convergence(const char* method, int iterations, double root, int max_iter);

int main(int argc, char const *argv[]) {
    double x0, x1, x2, error, tolerance;
    int max_iterations = 0, error_type, method;

    // X0 is the initial value from which the root is searched on the x axis, approaching with a straight line x increasingly closer to the root
    printf("Ingrese el valor inicial x0 desde donde empezará a buscar la raíz");
    scanf("%lf", &x0);
    printf("Usted necesita el error porcentual o el absoluto? (1 para absoluto, 0 para porcentual): ");
    scanf("%d", &error_type);
    printf("Ingrese la tolerancia");
    scanf("%lf", &tolerance);
    printf("Ingrese el tipo de metodo (1 para punto fijo, 2 para Newton-Raphson, 3 para secante): ");
    scanf("%d", &method);

    switch(method) {
        case 1: 
            // Fixed point method
            do {
                max_iterations++;

                // If the slope of the curve at the point where it intersects the x line with the function is greater --> it is not possible to find the root in a reasonable time
                if(fabs(gprima(x0)) >= 1) {
                    printf("El método no converge en la iteracion %d\n", max_iterations);
                    exit(0);
                }
            
                x1 = g(x0);
                error = calculate_error(x1, x0, error_type);
                x0 = x1; // We update x0 for the next iteration
            } while(error > tolerance);

            show_results("punto fijo", x1, error, max_iterations);
            break;

        case 2:
            // Newton-Raphson method
            do {
                max_iterations++;
        
                // If the derivative is too small, it may lead to division by zero or slow convergence
                if(gprima(x0) < 10e-4) {
                    printf("La derivada es muy pequeña en la iteración %d\n", max_iterations);
                    exit(0);
                }
        
                // We apply the Newton-Raphson formula
                x1 = x0 - (g(x0) / gprima(x0));
                error = calculate_error(x1, x0, error_type);
                x0 = x1; // We update x0 for the next iteration
            } while(error > tolerance && max_iterations < 10000);
        
            show_results("Newton-Raphson", x1, error, max_iterations);
            verify_convergence("Newton-Raphson", max_iterations, x1, 10000);
            break;

        case 3: 
            // Secant method - needs a second initial point (must be close to the root, like x0)
            printf("Ingrese un segundo valor inicial x1 para el método de la secante: ");
            scanf("%lf", &x1);
            
            do {
                max_iterations++;
        
                // We apply the secant method formula
                x2 = x1 - ((g(x1) * (x1 - x0)) / (g(x1) - g(x0)));
                error = calculate_error(x2, x1, error_type);
        
                x0 = x1; // We update x0 for the next iteration
                x1 = x2; // We update x1 for the next iteration
            } while(error > tolerance && max_iterations < 10000);
        
            show_results("la secante", x2, error, max_iterations);
            verify_convergence("la secante", max_iterations, x2, 10000);
            break;

        default:
            printf("Opción inválida. Seleccione 1, 2 o 3.\n");
    }
    
    return 0;
}

// Function definitions

// Function g(x) in the fixed-point method
// In the Netwon-Rapson method, f(x) is used; see the conceptual difference in a bibliography.
double g(double x) {
    // return cos(sin(x)); //Para segunda función g(n)
    return pow(x, x - cos(x)); //Para tercera función g(n)
}

// Numerical derivative of g(x) using finite differences
// In the Netwon-Rapson method, f'(x) is used; see the conceptual difference in a bibliography.
double gprima(double x) {
    return (g(x + 0.001) - g(x)) / 0.001;
}

// Function to calculate the error
double calculate_error(double new_x, double previous_x, int error_type) {
    if(error_type == 1) {
        // Absolute error
        return fabs(new_x - previous_x);
    } else {
        // Percentage error
        return fabs((new_x - previous_x) / new_x) * 100;
    }
}

// Function to display results
void show_results(const char* method, double root, double error, int iterations) {
    printf("La raiz aproximada es: %.12lf, con un error de %.12lf\n", root, error);
    printf("Número de iteraciones: %d\n", iterations);
    printf("La funcion evaluada en la raiz aproximada es igual a %.12lf\n", g(root));
}

// Function to verify convergence
// This function checks if the method converges and if the function value at the root is close to zero
void verify_convergence(const char* method, int iterations, double root, int max_iter) {
    double valor_funcion = g(root);
    
    if(iterations < max_iter) {
        if(fabs(valor_funcion) < 0.01) {
            printf("✓ El método de %s CONVERGE correctamente.\n", method);
            printf("  - Convergió en %d iteraciones (< %d)\n", iterations, max_iter);
            printf("  - g(raiz) = %.12lf está cerca de 0\n", valor_funcion);
        } else {
            printf("⚠ El método converge en iteraciones pero g(raiz) = %.12lf NO está cerca de 0.\n", valor_funcion);
            printf("  Posible problema: la raíz encontrada no es precisa.\n");
        }
    } else {
        printf("✗ El método de %s NO CONVERGE.\n", method);
        printf("  - Alcanzó el máximo de iteraciones (%d)\n", max_iter);
        printf("  - g(raiz) = %.12lf\n", valor_funcion);
    }
}
