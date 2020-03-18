import scipy.signal
import numpy as np
import matplotlib.pyplot as plt

colors = ["red", "blue", "green", "yellow", "orange", "brown", "purple"]

freqs_dict = dict()
amps_dict = dict()

bond_names = ["Au-C1" ,
              "C1-O1",
"Au-C2",
"C2-O2 ",
"Au-N2",
"N2-N1",
"N1-O3"]

data_by_bond = [[] for _ in range(7)]
angles = np.linspace(0, float(f"{np.pi:.5f}"), 251)

for angle in angles:
    try:
        filename = f"frequencies_{angle:0.3f}.txt"
        data = np.genfromtxt(filename, skip_header=1)
    except IOError:
        filename = f"frequencies_{angle-0.001:0.3f}.txt"
        data = np.genfromtxt(filename, skip_header=1)
    freqs = np.fft.rfftfreq(2*data.shape[0]-1, d=1).T
    freqs_dict[angle] = [[] for _ in range(7)]
    amps_dict[angle] = [[] for _ in range(7)]
    for i in range(7):
        data_by_bond[i].append(data[:, i])
        threshold = np.max(data[:, i])
        over_threshold = data[:, i] > np.max(data[:, 1]) / 100
        if (np.all(~over_threshold)):
            continue
        arg_rels = scipy.signal.argrelmax(data[:, i], axis=0, order=1)
        arg_rels_bool = np.zeros([data[:, i].shape[0]], dtype=bool)
        arg_rels_bool[arg_rels] = True
        max_and_over_threshold = np.logical_and(over_threshold, arg_rels_bool)
        print(freqs[arg_rels], data[:, i][arg_rels])
        freqs_dict[angle][i] = freqs[max_and_over_threshold]
        amps_dict[angle][i] = data[:, i][max_and_over_threshold]
        # plt.plot(freqs, data[:, i], color=colors[i])
        # plt.stem(freqs[max_and_over_threshold], data[:, i][max_and_over_threshold], colors[i], use_line_collection=True)

no_xlabels = 7 # how many labels to see on axis x
step_x = int(250 / (no_xlabels - 1)) # step between consecutive labels
x_positions = np.arange(0,250,step_x) # pixel count at label position
x_labels = freqs[::step_x] # labels you want to see

no_ylabels = 8 # how many labels to see on axis x
step_y = int(250 / (no_ylabels - 1)) # step between consecutive labels
y_positions = np.arange(0,250,step_y) # pixel count at label position
y_labels = angles[::step_y] / np.pi # labels you want to see




for index, bond in enumerate(data_by_bond):
    fig, ax = plt.subplots()
    bond = np.vstack(bond)
    ax.imshow(bond, aspect="auto", interpolation="nearest")
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{item:.2f}" for item in x_labels])
    ax.set_xlabel("Frequency")

    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"{item:.2f} π" for item in y_labels])
    ax.set_ylabel("N2-Au-C1 Angle")
    fig.savefig(f"./spectrum-{bond_names[index]}.pdf")
    plt.close(fig)
