#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100

#include "gauss.h"

/**
 * Function to define function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double f(double x);

int main(int argc, char const *argv[]) {
    int choice, subintervals;
    // Limits of integration
    double a, b;
    // Integral of Composition Trapeze
    double sum = 0.0;
    // Values to divide the interval of integration [a,b]
    double x = 0.0;
    // Distance between two consecutive points
    double h = 0.0;
    printf("Insert the limits of integration:\n");
    scanf("%lf %lf", &a, &b);
    /* printf("Please enter the number of subintervals:\n");
    scanf("%d", &subintervals); */
    double tolerance = 0.1;
    double error = 1000;
    double Igauss = 1.4626517536;
    subintervals = 2;
    do {
        sum = 0.0;
        // Calculate I
        sum = f(a) + f(b);
        h = (b - a) / subintervals;        
        for(int i = 1; i <= subintervals-1; i++) {
            x = a + i * h;
            sum += 2 * f(x);
        }
        sum = ((b-a)/(2*subintervals)) * sum;
        error = ((fabs(Igauss - sum))/fabs(Igauss)) * 100.0;
        if(error > tolerance) {
            subintervals++;
        }
    } while(error > tolerance);
    
    // Print integral using Trapeze
    printf("The integral is: %lf\n", sum);
    printf("The number of subintervals is: %d\n", subintervals);
    printf("The aproximated error is: %lf%%\n", error);

    return 0;
}
double f(double x) {
    return (exp(x*x));
}

