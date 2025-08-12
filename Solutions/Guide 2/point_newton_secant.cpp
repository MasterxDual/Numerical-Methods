#include <iostream>
#include <cmath>

using namespace std;

// Function g(x) in the fixed-point method
// In the Netwon-Rapson method, f(x) is used; see the conceptual difference in a bibliography.
double g(double x) {
    return x * x - 2 * x + 1;  // Ejemplo: x² - 2x + 1
}

// Numerical derivative of g(x) using finite differences
// In the Netwon-Rapson method, f'(x) is used; see the conceptual difference in a bibliography.
double gprima(double x) {
    return (g(x + 0.001) - g(x)) / 0.001;
}

// Function to calculate the error
double calcular_error(double x_nuevo, double x_anterior, int error_type) {
    if(error_type == 1) {
        // Absolute error
        return fabs(x_nuevo - x_anterior);
    } else {
        // Percentage error
        return fabs((x_nuevo - x_anterior) / x_nuevo) * 100;
    }
}

// Function to display results
void mostrar_resultados(const char* metodo, double raiz, double error, int iteraciones) {
    printf("La raiz aproximada es: %lf, con un error de %lf\n", raiz, error);
    printf("Número de iteraciones: %d\n", iteraciones);
    printf("La funcion evaluada en la raiz aproximada es igual a %lf\n", g(raiz));
}

// Function to verify convergence
// This function checks if the method converges and if the function value at the root is close to zero
void verificar_convergencia(const char* metodo, int iteraciones, double raiz, int max_iter) {
    double valor_funcion = g(raiz);
    
    if(iteraciones < max_iter) {
        if(fabs(valor_funcion) < 0.01) {
            printf("✓ El método de %s CONVERGE correctamente.\n", metodo);
            printf("  - Convergió en %d iteraciones (< %d)\n", iteraciones, max_iter);
            printf("  - g(raiz) = %lf está cerca de 0\n", valor_funcion);
        } else {
            printf("⚠ El método converge en iteraciones pero g(raiz) = %lf NO está cerca de 0.\n", valor_funcion);
            printf("  Posible problema: la raíz encontrada no es precisa.\n");
        }
    } else {
        printf("✗ El método de %s NO CONVERGE.\n", metodo);
        printf("  - Alcanzó el máximo de iteraciones (%d)\n", max_iter);
        printf("  - g(raiz) = %lf\n", valor_funcion);
    }
}

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
                error = calcular_error(x1, x0, error_type);
                x0 = x1; // We update x0 for the next iteration
            } while(error > tolerance);

            mostrar_resultados("punto fijo", x1, error, max_iterations);
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
                error = calcular_error(x1, x0, error_type);
                x0 = x1; // We update x0 for the next iteration
            } while(error > tolerance && max_iterations < 10000);
        
            mostrar_resultados("Newton-Raphson", x1, error, max_iterations);
            verificar_convergencia("Newton-Raphson", max_iterations, x1, 10000);
            break;

        case 3: 
            // Secant method - needs a second initial point (must be close to the root, like x0)
            printf("Ingrese un segundo valor inicial x1 para el método de la secante: ");
            scanf("%lf", &x1);
            
            do {
                max_iterations++;
        
                // We apply the secant method formula
                x2 = x1 - ((g(x1) * (x1 - x0)) / (g(x1) - g(x0)));
                error = calcular_error(x2, x1, error_type);
        
                x0 = x1; // We update x0 for the next iteration
                x1 = x2; // We update x1 for the next iteration
            } while(error > tolerance && max_iterations < 10000);
        
            mostrar_resultados("la secante", x2, error, max_iterations);
            verificar_convergencia("la secante", max_iterations, x2, 10000);
            break;

        default:
            printf("Opción inválida. Seleccione 1, 2 o 3.\n");
    }
    
    return 0;
}
