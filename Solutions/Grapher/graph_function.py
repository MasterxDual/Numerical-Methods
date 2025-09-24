#!/usr/bin/env python3
"""
Simple Function Grapher for Root Finding
========================================================

This program allows you to graph simple mathematical functions to aid in the visual
finding of roots in numerical methods.

Basic Features:
- Function graphs on Cartesian axes
- Automatic detection of approximate roots
- Simple and straightforward interface

Author: TOBIAS FUNES and BLACKBOXAI
Version: 2.0
"""

import matplotlib.pyplot as plt
import numpy as np

class SimpleGrapher:
    """Simplified class for graphing mathematical functions."""

    def __init__(self):
        """Initializes the grapher."""
        self.basic_functions = {
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'exp': np.exp, 'log': np.log, 'sqrt': np.sqrt,
            'pi': np.pi, 'e': np.e, 'abs': np.abs
        }

    def evaluate_secure_function(self, expr: str, x_vals: np.ndarray) -> np.ndarray:
        """
        Evaluates a function safely using numpy.

        Args:
            expr: String with the mathematical expression
            x_vals: Array of x values

        Returns:
            Array with the function values  
        """
        # Create a safe environment
        safe_dict = {
            'x': x_vals,
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'exp': np.exp, 'log': np.log, 'sqrt': np.sqrt,
            'pi': np.pi, 'e': np.e, 'abs': np.abs
        }

        # Evaluate using numpy for vectorization
        # Use safer evaluation with specific exceptions
        try:
            # Create a safe environment with allowed names
            allowed_names = {
                k: v for k, v in safe_dict.items() 
                if k not in ['__builtins__', 'eval', 'exec', 'import', 'open']
            }
            
            # Evaluate the expression safely
            # pylint: disable=eval-used
            y_vals = eval(expr, {"__builtins__": {}}, allowed_names)
            
            # Verify that the result is numeric
            if not isinstance(y_vals, np.ndarray):
                y_vals = np.array(y_vals)

            return y_vals
            
        except (ValueError, SyntaxError, NameError, ZeroDivisionError) as e:
            print(f"Error evaluating function '{expr}': {str(e)}")
            return np.full_like(x_vals, np.nan)

    def graficar(self, expr: str, x_min: float = -10, x_max: float = 10,
                 points_number: int = 1000) -> None:
        """
        Graphs a function and shows the approximate roots.
        Args:
        expr: String containing the mathematical expression
        x_min: Minimum value of x
        x_max: Maximum value of x
        num_points: Number of points for the graph
                """
        try:
            # Generate x points
            x_vals = np.linspace(x_min, x_max, points_number)

            # Evaluate the function
            y_vals = self.evaluate_secure_function(expr, x_vals)

            # Check for valid data
            if np.all(np.isnan(y_vals)):
                print("Error: No se pudo evaluar la función en el rango dado.")
                return

            # Create the figure
            plt.figure(figsize=(10, 6))

            # Plot the function
            plt.plot(x_vals, y_vals, 'b-', linewidth=2, label=f'f(x) = {expr}')

            # Configure axes
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)

            # Find approximate roots (where it crosses the x-axis)
            # Only consider non-NaN values
            valids = ~np.isnan(y_vals)
            if np.any(valids):
                y_valids = y_vals[valids]
                x_valids = x_vals[valids]

                # Find sign changes
                signos = np.sign(y_valids)
                cambios_signo = np.where(signos[:-1] * signos[1:] < 0)[0]

                if len(cambios_signo) > 0:
                    roots = []
                    for idx in cambios_signo:
                        # Calculate root using simple linear interpolation
                        x1, x2 = x_valids[idx], x_valids[idx+1]
                        y1, y2 = y_valids[idx], y_valids[idx+1]
                        if abs(y2 - y1) > 1e-10:  # Avoid division by zero
                            x_root = x1 - y1 * (x2 - x1) / (y2 - y1)
                            roots.append(x_root)

                    # Mark found roots
                    if roots:
                        plt.scatter(roots, [0]*len(roots), color='red', s=100, zorder=5, marker='o')
                        for _, root in enumerate(roots):
                            plt.annotate(f'Raíz ≈ {root:.4f}',
                                       (root, 0),
                                       xytext=(0, 10),
                                       textcoords='offset points',
                                       fontsize=10,
                                       ha='center',
                                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))

            # Configure title and labels
            plt.title(f'Graph of f(x) = {expr}', fontsize=14, pad=20)
            plt.xlabel('x', fontsize=12)
            plt.ylabel('f(x)', fontsize=12)
            plt.grid(True, alpha=0.3)
            plt.legend(fontsize=10)

            # Adjust limits 
            y_min, y_max = np.nanmin(y_vals), np.nanmax(y_vals)
            if np.isfinite(y_min) and np.isfinite(y_max):
                margen = (y_max - y_min) * 0.1
                plt.ylim(y_min - margen, y_max + margen)

            plt.tight_layout()

            # Save the plot as an image
            archive_name = f"grafica_{expr.replace('**', 'pow').replace('/', 'div').replace('-', 'menos')}.png"
            plt.savefig(archive_name, dpi=150, bbox_inches='tight')
            print(f"Saved graph as: {archive_name}")
            print("Note: In this environment, graphs are saved as PNG files.")
            print("To view the graph, open the generated file.")

            # Close the plot to free memory
            plt.close()

        except (ValueError, TypeError, RuntimeError, OverflowError) as e:
            print(f"Error al graficar: {str(e)}")

def show_instructions():
    """Shows basic usage instructions."""
    print("\n" + "="*50)
    print("     SIMPLE FUNCTION GRAPHER")
    print("="*50)
    print("Function examples:")
    print("  x**2 - 4          (root in x = 2)")
    print("  x**3 - x - 2      (root in x ≈ 1.52)")
    print("  exp(-x) - x       (root in x ≈ 0.57)")
    print("  cos(x)            (root in x ≈ 0.74)")
    print("  sin(x)            (root in x = 0)")
    print("="*50)

def main():
    """Main function of the program."""
    grapher = SimpleGrapher()

    print("Simple Function Grapher for Root Localization")
    print("This program helps you visualize functions and find roots.")

    show_instructions()

    while True:
        try:
            # Obtain function from user
            expr = input("\nEnter the function f(x): ").strip()

            if not expr:
                print("Error: You must enter a function.")
                continue

            # Obtain range 
            try:
                x_min_str = input("Minimum x value (default -10): ").strip()
                x_min = float(x_min_str) if x_min_str else -10.0

                x_max_str = input("Maximum x value (default 10): ").strip()
                x_max = float(x_max_str) if x_max_str else 10.0

                if x_min >= x_max:
                    print("Error: The minimum value must be less than the maximum value.")
                    continue

            except ValueError:
                print("Error: Please enter valid numbers for the range.")
                continue

            print("Generating graph...")
            grapher.graph(expr, x_min, x_max)

            # Ask if the user wants to continue
            continuar = input("\nWould you like to graph another function? (y/n): ").strip().lower()
            if continuar not in ['y', 'yes']:
                print("Goodbye!")
                break

        except KeyboardInterrupt:
            print("\n\nProgram interrupted by the user.")
            break
        except (ValueError, TypeError) as e:
            print(f"Unexpected error: {e}")
            print("Try with a simple function.")

if __name__ == "__main__":
    main()