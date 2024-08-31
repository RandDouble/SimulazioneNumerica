import shutil
from pathlib import Path
from typing import List, Any

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from scripts.preparation import Phase


def save_result_eq(phase: Phase) -> None:
    phase_name = phase.name.lower()
    file_list_to_copy = [
        "kinetic_energy.dat",
        "output.dat",
        "potential_energy.dat",
        "pressure.dat",
        "seed.out",
        "temperature.dat",
        "total_energy.dat",
        "CONFIG/velocities.out",
        "CONFIG/config.xyz",
    ]

    starting_dir = Path("NSL_SIMULATOR/OUTPUT")

    res_dir = Path(f"eser_4/equilibration_result/{phase_name}")
    res_dir.mkdir(parents=True, exist_ok=True)

    for file in file_list_to_copy:
        res_file = file.removeprefix("CONFIG/")
        shutil.copy(
            starting_dir / file,
            res_dir / res_file,
        )

        print(f"Copying {file} to eser_4/equilibration_result/{phase_name}/{res_file}")

        if "CONFIG" in file:
            res_file = res_file.replace(".xyz", f"_{phase_name}.xyz").replace(
                ".out", f"_{phase_name}.out"
            )
            shutil.copy(starting_dir / file, f"eser_4/measure_conf/{res_file}")
            print(f"Copying {file} to eser_4/measure_conf/{res_file}")

    print(f"Results for phase {phase_name} saved")


def save_result_meas(phase: Phase) -> None:
    phase_name = phase.name.lower()
    file_list_to_copy = [
        "kinetic_energy.dat",
        "output.dat",
        "potential_energy.dat",
        "pressure.dat",
        "seed.out",
        "temperature.dat",
        "total_energy.dat",
        "CONFIG/velocities.out",
        "CONFIG/config.xyz",
    ]

    starting_dir = Path("NSL_SIMULATOR/OUTPUT")

    res_dir = Path(f"eser_4/measure_result/{phase_name}")
    res_dir.mkdir(parents=True, exist_ok=True)

    for file in file_list_to_copy:
        res_file = file.removeprefix("CONFIG/")
        shutil.copy(
            starting_dir / file,
            res_dir / res_file,
        )

        print(f"Copying {file} to eser_4/measure_result/{phase_name}/{res_file}")

    print(f"Results for phase {phase_name} saved")


def load_resulting_data(data_dir: Path) -> List[pd.DataFrame]:
    common_configs = {
        "header": None,
        "sep": "\s+",
        "skipinitialspace": True,
        "keep_default_na": False,
        "skiprows": 1,
    }
    temperature = pd.read_csv(
        data_dir / "temperature.dat",
        names=["block", "actual_t", "t_ave", "error"],
        **common_configs,
    )
    kinetic_energy = pd.read_csv(
        data_dir / "kinetic_energy.dat",
        names=["block", "actual_ke", "ke_ave", "error"],
        **common_configs,
    )
    potential_energy = pd.read_csv(
        data_dir / "potential_energy.dat",
        names=["block", "actual_pe", "pe_ave", "error"],
        **common_configs,
    )
    total_energy = pd.read_csv(
        data_dir / "total_energy.dat",
        names=["block", "actual_te", "te_ave", "error"],
        **common_configs,
    )
    pressure = pd.read_csv(
        data_dir / "pressure.dat",
        names=["block", "actual_p", "p_ave", "error"],
        **common_configs,
    )

    return temperature, kinetic_energy, potential_energy, total_energy, pressure


def print_resulting_data(*dataframes) -> list[Figure, Any]:
    name_conversion = {
        "actual_t": "Temperature [K]",
        "actual_ke": "Kinetic Energy [J]",
        "actual_pe": "Potential Energy [J]",
        "actual_te": "Total Energy [J]",
        "actual_p": "Pressure [Pa]",
    }

    title_conversion = {
        "actual_t": "Temperature",
        "actual_ke": "Kinetic Energy",
        "actual_pe": "Potential Energy",
        "actual_te": "Total Energy",
        "actual_p": "Pressure",
    }

    fig, axs = plt.subplots(
        1, len(dataframes), figsize=(16, 5), sharex=True, layout="constrained"
    )

    for ax, df in zip(axs, dataframes):
        ax.plot(
            df.iloc[:, 0],
            df.iloc[:, 1],
            label="Istantaneous",
            marker="o",
            linestyle="",
            c="C0",
            zorder=0,
        )
        ax.plot(df.iloc[:, 0], df.iloc[:, 2], c="C1", label="Average", zorder=2)
        ax.fill_between(
            df.iloc[:, 0],
            df.iloc[:, 2] - df.iloc[:, 3],
            df.iloc[:, 2] + df.iloc[:, 3],
            alpha=0.6,
            color="C1",
            zorder=1,
        )
        ax.set_xlabel("Block")
        # ax.set_xscale("log")
        ax.set_ylabel(name_conversion[df.columns[1]])
        ax.set_title(title_conversion[df.columns[1]])
        ax.legend()

    return fig, axs


def save_configurations(phase: Phase) -> None:
    phase_name = phase.name.lower()
    dir_to_copy = Path("NSL_SIMULATOR/OUTPUT/CONFIG")

    eq_res_dir = Path(f"eser_4/equilibration_result/{phase_name}_config")

    eq_res_dir.mkdir(parents=True, exist_ok=True)

    shutil.copytree(dir_to_copy, eq_res_dir)
    print(f"Copying CONFIG to eser_4/equilibration_result/{phase_name}_config")

    print(f"Results for phase {phase_name} saved")
