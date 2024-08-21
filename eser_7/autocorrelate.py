import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.signal as signal
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def autocorrelation(data: np.ndarray) -> np.ndarray:
    mean = data.mean()
    v_autocorrelation = signal.correlate(data - mean, data - mean, mode="same")
    v_autocorrelation /= v_autocorrelation.max()
    return v_autocorrelation


def exp_func(x, tau, a):
    return a * np.exp(-x / tau)


def fit_exp_autocorrelation(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(data[data.size // 2 : -1].size)
    popt, pcov = optimize.curve_fit(exp_func, x, data[data.size // 2 : -1])
    return popt, pcov


def print_autocorr(autocorr: np.ndarray, popt: np.ndarray) -> tuple[Figure, Axes]:
    fig: Figure
    ax: Axes
    fig, ax = plt.subplots()
    time_step = np.arange(autocorr[autocorr.size // 2 : -1].size)
    ax.plot(time_step, autocorr[autocorr.size // 2 : -1], label="Autocorrelation")
    ax.plot(
        time_step,
        exp_func(time_step, *popt),
        label="Fitted curve," + r" $\tau$ :" + f" {popt[0]:.2f}",
    )
    ax.set_ylabel("Autocorrelation")
    ax.set_xlabel("steps")
    ax.set_title("Autocorrelation vs time steps")
    return fig, ax


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file", required=True, help="Name of the file to be analyzed, can be only csv"
    )
    parser.add_argument(
        "idx", required=False, default=1, help="Column index for autocorrelation"
    )
    args = parser.parse_args()
    file = args.file
    conv_idx = args.idx

    df = pd.read_csv(file)

    data = df.iloc[:, conv_idx].to_numpy()

    autocorr = autocorrelation(data)
    popt, _ = fit_exp_autocorrelation(autocorr)

    fig: Figure
    fig, _ = print_autocorr(autocorr, popt)
    fig.savefig("autocorrelation.svg")


if __name__ == "__main__":
    main()
