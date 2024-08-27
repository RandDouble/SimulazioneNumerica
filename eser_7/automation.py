import shutil
from pathlib import Path
from scripts.preparation import Phase

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Any


def save_results(
    phase: Phase,
    file_list: list[str],
    save_dir: Path,
    output_dir: Path = Path("NSL_SIMULATOR/OUTPUT"),
) -> None:
    phase_name = phase.name.lower()

    res_dir = save_dir / phase_name
    res_dir.mkdir(parents=True, exist_ok=True)

    for file in file_list:
        res_file = file.removeprefix("CONFIG/").replace(".xyz", f"_{phase_name}.xyz")
        shutil.copy(output_dir / file, res_dir / res_file)
        print(f"Copying {file} to {res_dir / file}")
    print(f"Results for phase {phase_name} saved")


def save_results_starting_equilibrated(phase: Phase):
    file_list_to_save = [
        "potential_energy.dat",
        "acceptance.dat",
        "gofr.dat",
        "partial_gofr.dat",
    ]
    save_dir = Path("eser_7/from_molecular_dynamics/")
    save_results(phase, file_list_to_save, save_dir)


def save_results_equilibration(phase: Phase):
    file_list_to_save = ["potential_energy.dat", "acceptance.dat", "CONFIG/config.xyz"]
    save_dir = Path("eser_7/from_monte_carlo")
    save_results(phase, file_list_to_save, save_dir)


def save_results_production(phase: Phase):
    file_list_to_save = [
        "potential_energy.dat",
        "acceptance.dat",
        "gofr.dat",
        "partial_gofr.dat",
        "pressure.dat",
    ]
    save_dir = Path("eser_7/measurement")
    save_results(phase, file_list_to_save, save_dir)


def save_results_production_MD(phase: Phase):
    file_list_to_save = [
        "potential_energy.dat",
        "acceptance.dat",
        "gofr.dat",
        "partial_gofr.dat",
        "pressure.dat",
    ]
    save_dir = Path("eser_7/measurement_MD")
    save_results(phase, file_list_to_save, save_dir)


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


def load_pressure_data(data_dir: Path) -> pd.DataFrame:
    common_config = {
        "header": None,
        "sep": "\s+",
        "skipinitialspace": True,
        "keep_default_na": False,
        "skiprows": 1,
    }
    pressure_df = pd.read_csv(
        data_dir / "pressure.dat",
        names=["block", "actual_p", "p_ave", "error"],
        **common_config,
    )
    return pressure_df


def load_gofr_data(data_dir: Path) -> pd.DataFrame:
    common_config = {
        "header": None,
        "sep": "\s+",
        "skipinitialspace": True,
        "keep_default_na": False,
        "skiprows": 1,
    }
    gofr_df = pd.read_csv(
        data_dir / "gofr.dat",
        names=["distance", "average", "error"],
        **common_config,
    )
    return gofr_df


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
