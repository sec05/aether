import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import re

mesh_sizes = [2**k for k in range(3, 10)]  # Mesh sizes from 8 to 512
l2_errors = []
linf_errors = []

PATH_TO_EXE = os.path.join(os.path.dirname(__file__), "../build/src/aether")
PATH_TO_INPUT = os.path.join(os.path.dirname(__file__), "../input.conf")

for mesh_size in mesh_sizes:
    with open(PATH_TO_INPUT, "r+") as f:
        content = f.read().splitlines()
        for i, line in enumerate(content):
            if line[0:2] == "Nx":
                content[i] = f"Nx: {mesh_size}"
            elif line[0:2] == "Ny":
                content[i] = f"Ny: {mesh_size}"
        content = "\n".join(content)
        f.seek(0)
        f.write(content)
        f.truncate()

    result = subprocess.run([PATH_TO_EXE], capture_output=True, text=True)
    output = result.stdout

    m = re.search(r"L2\s*=\s*([\d.eE+-]+)\s+Linf\s*=\s*([\d.eE+-]+)", output)
    if m:
        l2_errors.append(float(m.group(1)))
        linf_errors.append(float(m.group(2)))

mesh_sizes = np.array(mesh_sizes, dtype=float)
l2_errors = np.array(l2_errors)
linf_errors = np.array(linf_errors)

# error ~ C * N^(-p)  =>  log(err) = log(C) - p*log(N); slope = -p
l2_slope, _ = np.polyfit(np.log(mesh_sizes), np.log(l2_errors), 1)
linf_slope, _ = np.polyfit(np.log(mesh_sizes), np.log(linf_errors), 1)
l2_rate, linf_rate = -l2_slope, -linf_slope
print(f"L2 convergence rate:   {l2_rate:.3f}")
print(f"Linf convergence rate: {linf_rate:.3f}")

# Theoretical O(h^2) reference (h ~ 1/N), anchored to the coarsest mesh
ref_order = 4/3
l2_ref = l2_errors[0] * (mesh_sizes / mesh_sizes[0]) ** (-ref_order)
linf_ref = linf_errors[0] * (mesh_sizes / mesh_sizes[0]) ** (-ref_order)

plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.loglog(mesh_sizes, l2_errors, marker='o', label=f'L2 error (rate {l2_rate:.2f})')
plt.loglog(mesh_sizes, l2_ref, 'k--', label=f'$O(h^{ref_order})$ reference')
plt.xlabel('Mesh Size (Nx = Ny)')
plt.ylabel('Error')
plt.title('L2 Error Convergence')
plt.legend()
plt.grid(True, which="both", ls="--")

plt.subplot(1, 2, 2)
plt.loglog(mesh_sizes, linf_errors, marker='o', label=f'Linf error (rate {linf_rate:.2f})')
plt.loglog(mesh_sizes, linf_ref, 'k--', label=f'$O(h^{ref_order})$ reference')
plt.xlabel('Mesh Size (Nx = Ny)')
plt.ylabel('Error')
plt.title('Linf Error Convergence')
plt.legend()
plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('convergence_plot.png', dpi=150)