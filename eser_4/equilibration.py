from subprocess import Popen, PIPE, STDOUT
import shutil
from pathlib import Path
from enum import Enum, auto
import io
import os
from typing import List, Any

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


class Phase(Enum):
    SOLID = auto()
    LIQUID = auto()
    GAS = auto()


def compile_program() -> None:
    with Popen(
        "make -C NSL_SIMULATOR/SOURCE/".split(), stdout=PIPE, stderr=STDOUT
    ) as proc:
        for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
            print(line)
        # print(proc.stdout.read())
    print("Program compiled")


def load_config(
    phase: Phase,
    src_dir: Path = Path("eser_4/equilibration_conf"),
    reset_flag: bool = False,
) -> None:
    phase_name = phase.name.lower()
    print(f"Actually working on phase : {phase.name}")
    dest_dir = Path("NSL_SIMULATOR/INPUT")

    shutil.copy(src_dir / f"input.{phase_name}", dest_dir / "input.dat")
    print(f"Configuration for phase {phase_name} loaded")
    shutil.copy(src_dir / "properties.dat", dest_dir / "properties.dat")
    print("Properties loaded")

    if reset_flag:
        shutil.copy(
            src_dir / f"config_{phase_name}.xyz", dest_dir / "CONFIG/config.xyz"
        )
        shutil.copy(
            src_dir / f"velocities_{phase_name}.out", dest_dir / "CONFIG/velocities.in"
        )
    else:
        shutil.copy(src_dir / "config.xyz", dest_dir / "CONFIG/config.xyz")


def run_program(exe_dir: Path, save_config: bool = False) -> None:
    cur_dir = Path(os.getcwd())
    # Changing directory to the executable directory because NSL_SIMULATOR depends on relative paths.
    os.chdir(exe_dir)
    print(os.getcwd())
    command = ["./simulator.exe"]
    if save_config:
        command.append("--write-config")

    print(command)

    with Popen(
        command,
        stdout=PIPE,
        stderr=STDOUT,  # , cwd=os.getcwd()
    ) as proc:
        for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
            print(line)

    os.chdir(cur_dir)
    print("Program executed")


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

    os.makedirs(f"eser_4/equilibration_result/{phase_name}", exist_ok=True)

    for file in file_list_to_copy:
        res_file = file.removeprefix("CONFIG/")
        shutil.copy(
            f"NSL_SIMULATOR/OUTPUT/{file}",
            f"eser_4/equilibration_result/{phase_name}/{res_file}",
        )

        print(f"Copying {file} to eser_4/equilibration_result/{phase_name}/{res_file}")

        if "CONFIG" in file:
            res_file = res_file.replace(".xyz", f"_{phase_name}.xyz").replace(
                ".out", f"_{phase_name}.out"
            )
            shutil.copy(
                f"NSL_SIMULATOR/OUTPUT/{file}", f"eser_4/measure_conf/{res_file}"
            )
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

    os.makedirs(f"eser_4/measure_result/{phase_name}", exist_ok=True)

    for file in file_list_to_copy:
        res_file = file.removeprefix("CONFIG/")
        shutil.copy(
            f"NSL_SIMULATOR/OUTPUT/{file}",
            f"eser_4/measure_result/{phase_name}/{res_file}",
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
            df.iloc[:, 0], df.iloc[:, 1], label="Istantaneous", marker="o", linestyle=""
        )
        ax.errorbar(df.iloc[:, 0], df.iloc[:, 2], yerr=df.iloc[:, 3], label="Average")
        ax.set_xlabel("Block")
        # ax.set_xscale("log")
        ax.set_ylabel(name_conversion[df.columns[1]])
        ax.set_title(title_conversion[df.columns[1]])
        ax.legend()

    return fig, axs


def save_configurations(phase: Phase) -> None:
    phase_name = phase.name.lower()
    dir_to_copy = "CONFIG/"

    os.makedirs(f"eser_4/equilibration_result/{phase_name}_config", exist_ok=True)

    shutil.copytree(
        f"NSL_SIMULATOR/OUTPUT/{dir_to_copy}",
        f"eser_4/equilibration_result/{phase_name}_config",
    )
    print(f"Copying CONFIG to eser_4/equilibration_result/{phase_name}_config")

    print(f"Results for phase {phase_name} saved")
