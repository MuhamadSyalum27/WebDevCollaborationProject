import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tkinter import Tk, Label, Entry, Button, StringVar

# Fungsi persamaan diferensial phugoid
def phugoid(y, t, omega_n, zeta, delta_alpha):
    u, u_dot = y
    u_dot_dot = -2 * zeta * omega_n * u_dot - omega_n**2 * u + delta_alpha
    return [u_dot, u_dot_dot]

# Fungsi untuk menggambar grafik phugoid
def gambar_grafik(u, u_dot, t, stable_point_index):
    plt.figure(figsize=(12, 6))

    plt.subplot(2, 1, 1)
    plt.plot(t, u)
    plt.title('Phugoid Mode: Kecepatan terhadap Waktu')
    plt.xlabel('Waktu (s)')
    plt.ylabel('Kecepatan (u)')

    plt.subplot(2, 1, 2)
    plt.plot(t, u_dot)
    plt.title('Phugoid Mode: Kecepatan Sudut terhadap Waktu')
    plt.xlabel('Waktu (s)')
    plt.ylabel('Kecepatan Sudut (u_dot)')

    # Menampilkan grafik output
    plt.legend()
    plt.tight_layout()
    plt.show()

# Fungsi untuk memulai simulasi phugoid
def start_phugoid_simulation(omega_n, zeta, delta_alpha, u0, u_dot0, simulation_time):
    t = np.linspace(0, simulation_time, 3000)
    solution = odeint(phugoid, [u0, u_dot0], t, args=(omega_n, zeta, delta_alpha))
    u, u_dot = solution.T
    gambar_grafik(u, u_dot, t, 0)

# Fungsi yang akan dipanggil saat tombol "Simulasi" ditekan
def simulasi_phugoid():
    omega_n = float(omega_n_entry.get())
    zeta = float(zeta_entry.get())
    delta_alpha = float(delta_alpha_entry.get())
    u0 = float(u0_entry.get())
    u_dot0 = float(u_dot0_entry.get())
    simulation_time = float(simulation_time_entry.get())
    start_phugoid_simulation(omega_n, zeta, delta_alpha, u0, u_dot0, simulation_time)

# Membuat GUI menggunakan Tkinter
root = Tk()
root.title("Simulasi Phugoid")

# Label dan Entry untuk omega_n
Label(root, text="Natural frequency (rad/s):").grid(row=0, column=0)
omega_n_var = StringVar(value="")
omega_n_entry = Entry(root, textvariable=omega_n_var)
omega_n_entry.grid(row=0, column=1)

# Label dan Entry untuk zeta
Label(root, text="Damping ratio:").grid(row=1, column=0)
zeta_var = StringVar(value="")
zeta_entry = Entry(root, textvariable=zeta_var)
zeta_entry.grid(row=1, column=1)

# Label dan Entry untuk delta_alpha
Label(root, text="Gangguan sudut serang:").grid(row=2, column=0)
delta_alpha_var = StringVar(value="")
delta_alpha_entry = Entry(root, textvariable=delta_alpha_var)
delta_alpha_entry.grid(row=2, column=1)

# Label dan Entry untuk u0
Label(root, text="Kecepatan Awal:").grid(row=3, column=0)
u0_var = StringVar(value="")
u0_entry = Entry(root, textvariable=u0_var)
u0_entry.grid(row=3, column=1)

# Label dan Entry untuk u_dot0
Label(root, text="Kecepatan Sudut Awal:").grid(row=4, column=0)
u_dot0_var = StringVar(value="")
u_dot0_entry = Entry(root, textvariable=u_dot0_var)
u_dot0_entry.grid(row=4, column=1)

# Label dan Entry untuk waktu simulasi
Label(root, text="Waktu Simulasi (s):").grid(row=5, column=0)
simulation_time_var = StringVar(value="")
simulation_time_entry = Entry(root, textvariable=simulation_time_var)
simulation_time_entry.grid(row=5, column=1)

# Tombol untuk memulai simulasi
Tombol_simulasi = Button(root, text="Simulasi", command=simulasi_phugoid)
Tombol_simulasi.grid(row=6, column=0, columnspan=2)

# Menampilkan GUI
root.mainloop()
