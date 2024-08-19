from subprocess import Popen, PIPE, STDOUT
import shutil
from pathlib import Path
from enum import Enum, auto
import io
import os


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


def load_config(phase: Phase) -> None:
    phase_name = phase.name.lower()
    print(f"Actually working on phase : {phase.name}")
    src_dir = Path("eser_4/equilibration_conf")
    dest_dir = Path("NSL_SIMULATOR/INPUT")

    shutil.copy(src_dir / f"input.{phase_name}", dest_dir / "input.dat")
    print(f"Configuration for phase {phase_name} loaded")
    shutil.copy(src_dir / "properties.dat", dest_dir / "properties.dat")
    print("Properties loaded")


def run_program(exe_dir: Path) -> None:
    cur_dir = Path(os.getcwd())
    # Changing directory to the executable directory because NSL_SIMULATOR depends on relative paths.
    os.chdir(exe_dir)
    print(os.getcwd())
    with Popen(
        "./simulator.exe",
        stdout=PIPE,
        stderr=STDOUT,  # , cwd=os.getcwd()
    ) as proc:
        for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
            print(line)

    os.chdir(cur_dir)
    print("Program executed")


def save_result(phase: Phase) -> None:
    phase_name = phase.name.lower()
    file_list_to_copy = [
        "kinetic_energy.dat",
        "output.dat",
        "potential_energy.dat",
        "pressure.dat",
        "seed.out",
        "temperature.dat",
        "total_energy.dat",
    ]

    os.makedirs(f"eser_4/equilibration_result/{phase_name}", exist_ok=True)

    for file in file_list_to_copy:
        shutil.copy(
            f"NSL_SIMULATOR/OUTPUT/{file}",
            f"eser_4/equilibration_result/{phase_name}/{file}",
        )
        print(f"Copying {file} to eser_4/equilibration_result/{phase_name}/{file}")

    print(f"Results for phase {phase_name} saved")
