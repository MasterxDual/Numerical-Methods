#include <iostream>
#include <cmath>

using namespace std;

/* Looking for the graph of f(x) = 30809.6 * pow(x, 2) - 1026.99 * pow(x, 3) - 2620861.9
    We can see that there's three roots. We discard negative root because we need to find a depth (it's positive)
    We'll use the Newton-Raphson method:
    Aproximately root is: 11.861, with an error of 0.000
    Number of iterations: 3
    Function evalutaed in the aproximately root is equal to 0.000
    ✓ Method converges Newton-Raphson correctly.
      - Converged in 3 iterations (< 10000)
      - f(root) = 0.000 is close to 0
    
    Aproximately root is: 26.314, with an error of 0.000
    Number of iterations: 3
    Function evalutaed in the aproximately root is equal to -0.001
    ✓ Method converges Newton-Raphson correctly.
      - Converged in 3 iterations (< 10000)
      - f(root) = -0.001 is close to 0

    The first root is the one we need to find the depth of the sphere because 0 < d < 2 * r = 20 cm
    So d = 11.861 cm 
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
    printf("Enter the initial value x0 from where you will start searching for the root");
    scanf("%lf", &x0);
    printf("Do you need the percentage or absolute error? (1 for absolute, 0 for percentage): ");
    scanf("%d", &error_type);
    printf("Enter the tolerance");
    scanf("%lf", &tolerance);
    printf("Enter the method type (1 for fixed point, 2 for Newton-Raphson, 3 for secant): ");
    scanf("%d", &method);

    switch(method) {
        case 1: 
            // Fixed point method
            do {
                max_iterations++;

                // If the slope of the curve at the point where it intersects the x line with the function is greater --> it is not possible to find the root in a reasonable time
                if(fabs(gprima(x0)) >= 1) {
                    printf("The method does not converge in the iteration %d\n", max_iterations);
                    exit(0);
                }
            
                x1 = g(x0);
                error = calculate_error(x1, x0, error_type);
                x0 = x1; // We update x0 for the next iterationx_anterior
            } while(error > tolerance);

            show_results("punto fijo", x1, error, max_iterations);
            break;

        case 2:
            // Newton-Raphson method
            do {
                max_iterations++;
        
                // If the derivative is too small, it may lead to division by zero or slow convergence
                if(fabs(gprima(x0)) < 10e-4) {
                    printf("The derivative is very small in the iteration %d\n", max_iterations);
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
            printf("Enter a second initial value x1 for the secant method: ");
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
            printf("Invalid option. Select 1, 2 o 3.\n");
    }
    
    return 0;
}

// Function definitions

// Function g(x) in the fixed-point method
// In the Netwon-Rapson method, f(x) is used; see the conceptual difference in a bibliography.
double g(double x) {
    return 30809.6 * pow(x, 2) - 1026.99 * pow(x, 3) - 2620861.9;
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
    printf("Aproximately root is: %.3lf, with an error of %.3lf\n", root, error);
    printf("Number of iterations: %d\n", iterations);
    if(method != "punto fijo") {
        printf("Function evalutaed in the aproximately root is equal to %.3lf\n", g(root));
    }
}

// Function to verify convergence
// This function checks if the method converges and if the function value at the root is close to zero
void verify_convergence(const char* method, int iterations, double root, int max_iter) {
    double valor_funcion = g(root);
    
    if(iterations < max_iter) {
        if(fabs(valor_funcion) < 0.01) {
            printf("✓ Method converges %s correctly.\n", method);
            printf("  - Converged in %d iterations (< %d)\n", iterations, max_iter);
            printf("  - f(root) = %.3lf is close to 0\n", valor_funcion);
        } else {
            printf("⚠ Method converges in iterations but g(root) = %.3lf is NOT close to 0.\n", valor_funcion);
            printf("  Possible issue: the root found is not accurate.\n");
        }
    } else {
        printf("✗ Method %s don't converges.\n", method);
        printf("  - Reached maximum number of iterations (%d).\n", max_iter);
        printf("  - f(root) = %.3lf\n", valor_funcion);
    }
}
