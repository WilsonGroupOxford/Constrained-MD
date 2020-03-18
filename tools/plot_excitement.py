import numpy as np

import matplotlib.pyplot as plt


emass = 5.48579909065e-4
masses = {"Au": 196.9667,
          "C": 12.011,
          "N": 14.007,
          "O": 15.999,}
data = np.genfromtxt("bonds.csv", skip_header=1, dtype=float, delimiter=",")

with open("./bonds.csv") as fi:
    names = fi.readline().split(",")

fig, ax = plt.subplots()
ax.plot(range(data.shape[0]), data[:, 1], label=names[1])
ax.plot(range(data.shape[0]), data[:, 4], label=names[4], linestyle="dashed")
ax.legend()
fig.savefig("./bond-excitement.pdf")

linestyles = ["solid", "dashed", "dotted", "dashdot", "solid", "dashed", "dotted", "dashdot"]
AMPS_TABLE = np.empty([7, data.shape[0]//2 + 1])
for i in range(7):
    FFT = np.fft.rfft(data[:,i] - 1.0)
    FREQS = np.fft.rfftfreq(data.shape[0], d=1)
    AMPS = np.abs(FFT)
    name_a, name_b = names[i].split("->")
    mass_a, mass_b = masses[name_a.split("%")[0]], masses[name_b.split("%")[0]]
    print(mass_a, mass_b)
    reduced_mass = mass_a * mass_b / (mass_a + mass_b)
    FORCE_CONSTANTS = (2 * np.pi * FREQS)**2 * reduced_mass
    print(name_a, name_b, FORCE_CONSTANTS[AMPS > 100])
    plt.plot(FREQS, AMPS, label=names[i], linestyle=linestyles[i])
    AMPS_TABLE[i, :] = AMPS

np.savetxt("./frequencies.txt", AMPS_TABLE.T, header=", ".join(names))

plt.xlabel("Frequency (arb. units)")
plt.ylabel("Amplitude")
plt.xlim(0, 0.5)
plt.ylim(0, 50)
plt.legend()
plt.savefig("./fft.pdf")
