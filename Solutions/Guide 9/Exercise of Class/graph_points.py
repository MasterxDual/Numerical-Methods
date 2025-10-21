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
x_exact = np.linspace(0, 1, 200)
y_exact = (1/4.0) * ((x_exact**2) + x_exact + 2)**2

# --- Cargar soluciones numéricas ---
methods = {
    "Euler": ("euler_results.txt", "r", "o"),
    "Heun": ("heun_results.txt", "g", "s"),
    "Punto Medio": ("midpoint_results.txt", "m", "^"),
    "Runge-Kutta 4": ("runge_kutta_4_results.txt", "c", "D")
}

plt.figure(figsize=(9, 7))

# Graficar cada método que exista
for label, (filename, color, marker) in methods.items():
    x, y = load_data(filename)
    if x is not None:
        plt.plot(x, y, color=color, marker=marker, linestyle='-', label=label)

# Graficar solución exacta
plt.plot(x_exact, y_exact, 'b-', linewidth=2, label="Solución Exacta")

# --- Estética del gráfico ---
plt.title("Comparación de Métodos Numéricos vs Solución Exacta", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("y(x)", fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0, 1)
plt.tight_layout()

plt.show()
