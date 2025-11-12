import numpy as np
import matplotlib.pyplot as plt
import os

def load_data(filename):
    """Carga datos desde un archivo si existe"""
    if os.path.exists(filename):
        data = np.loadtxt(filename)
        return data[:,0], data[:,1]
    return None, None

# --- Solución exacta ---
# x_exact = np.linspace(0, 1, 200)
# y_exact = np.sqrt(np.exp(x_exact*x_exact))

# --- Cargar soluciones numéricas ---
methods = {
    "Método del Ejercicio": ("method_of_exercise_results.txt", "r", "o"),
    "Exacta": ("exact_results.txt", "g", "s"),
    "Runge-Kutta 4": ("runge_kutta_4_results.txt", "c", "D")
}

plt.figure(figsize=(9, 7))

# Graficar cada método que exista
for label, (filename, color, marker) in methods.items():
    x, y = load_data(filename)
    if x is not None:
        plt.plot(x, y, color=color, marker=marker, linestyle='-', label=label)

# Graficar solución exacta
# plt.plot(x_exact, y_exact, 'b-', linewidth=2, label="Solución Exacta")

# --- Estética del gráfico ---
plt.title("Comparación de Métodos Numéricos vs Solución Exacta", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("y(x)", fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0, 1)
plt.tight_layout()

plt.margins(x=0.05, y=0.05)  # Márgenes del 5% en X y 5% en Y

plt.show()
