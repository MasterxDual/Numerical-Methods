import matplotlib.pyplot as plt
import numpy as np

# Leer datos
bisection = np.loadtxt('bisection_convergence.txt', skiprows=4)
regula = np.loadtxt('regula_falsi_convergence.txt', skiprows=4)

# Crear figura con 2 subgráficos
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Gráfico 1: Comparativo
ax1.semilogy(bisection[:,0], bisection[:,1], 'b-o', label='Bisección', linewidth=2, markersize=4)
ax1.semilogy(regula[:,0], regula[:,1], 'r-s', label='Regula Falsi', linewidth=2, markersize=4)
ax1.set_xlabel('Número de Iteraciones')
ax1.set_ylabel('Error Absoluto (escala log)')
ax1.set_title('Convergencia Comparativa\nf(x) = x¹⁰ - 1')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Gráfico 2: Solo primeras iteraciones (escala lineal)
max_iter = min(15, len(bisection), len(regula))
ax2.plot(bisection[:max_iter,0], bisection[:max_iter,1], 'b-o', label='Bisección', linewidth=2, markersize=4)
ax2.plot(regula[:max_iter,0], regula[:max_iter,1], 'r-s', label='Regula Falsi', linewidth=2, markersize=4)
ax2.set_xlabel('Número de Iteraciones')
ax2.set_ylabel('Error Absoluto')
ax2.set_title(f'Primeras {max_iter} Iteraciones\n(escala lineal)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Ajustar diseño y mostrar
plt.tight_layout()
plt.savefig('convergence_analysis.png', dpi=300, bbox_inches='tight')
print(f"✅ Gráfico guardado como 'convergence_analysis.png'")
print(f"📊 Bisección: {len(bisection)} iteraciones")
print(f"📊 Regula Falsi: {len(regula)} iteraciones")
plt.show()