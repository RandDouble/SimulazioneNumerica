import shutil
from pathlib import Path
from scripts.preparation import Phase

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Any


def save_results_starting_equilibrated(phase: Phase):
    phase_name = phase.name.lower()
    file_list_to_save = [
        "potential_energy.dat",
        "acceptance.dat",
    ]
    res_dir = Path(f"eser_7/from_molecular_dynamics/{phase_name}")
    res_dir.mkdir(parents=True, exist_ok=True)

    for file in file_list_to_save:
        shutil.copy(f"NSL_SIMULATOR/OUTPUT/{file}", res_dir / file)
        print(f"Copying {file} to {res_dir / file}")
    print(f"Results for phase {phase_name} saved")


def load_energy_data(data_dir: Path) -> pd.DataFrame:
    common_config = {
        "header": None,
        "sep": "\s+",
        "skipinitialspace": True,
        "keep_default_na": False,
        "skiprows": 1,
    }
    potential_energy_df = pd.read_csv(
        data_dir / "potential_energy.dat",
        names=["block", "actual_pe", "pe_ave", "error"],
        **common_config,
    )
    return potential_energy_df


def print_energy(df: pd.DataFrame) -> list[Figure, Any]:
    name_conversion = {
        "actual_pe": "Potential Energy [J]",
    }

    title_conversion = {
        "actual_pe": "Potential Energy",
    }

    fig: Figure
    ax: Axes
    fig, ax = plt.subplots()

    ax.plot(
        df.iloc[:, 0], df.iloc[:, 1], label="Istantaneous", marker=".", linestyle=""
    )
    ax.errorbar(df.iloc[:, 0], df.iloc[:, 2], yerr=df.iloc[:, 3], label="Average")
    ax.set_xlabel("Block")
    ax.set_ylabel(name_conversion[df.columns[1]])
    ax.set_title(title_conversion[df.columns[1]])

    return fig, ax
