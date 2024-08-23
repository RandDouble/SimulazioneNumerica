"""
Summary
-------
This module contains functions that are used to compile and run the NSL_SIMULATOR program.
"""

import os
import shutil
import subprocess
from enum import Enum, auto
from pathlib import Path


class Phase(Enum):
    SOLID = auto()
    LIQUID = auto()
    GAS = auto()


def compile_program(*args) -> None:
    """
    Parameters
    ----------
    *args : tuple
        Variable number of arguments to be passed to the make command.

    Notes
    -----
    This function compiles a NSL_SIMULATOR by executing the make command in the NSL_SIMULATOR/SOURCE directory.
    The specified arguments are appended to the make command.

    Examples
    --------
    >>> compile_program("arg1", "arg2")
    Program compiled
    """
    command = "make -C NSL_SIMULATOR/SOURCE/".split()
    command.extend(args)
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
    ) as proc:
        for line in proc.stdout:
            print(line, end="")
        # print(proc.stdout.read())
    print("Program compiled")


def run_program(exe_dir: Path, save_config: bool = False) -> None:
    """

    Parameters
    ----------
    exe_dir : (Path)
        The path to the executable directory.
    save_config : (bool, optional)
        Whether to save the configuration. Defaults to False.
    Notes
    -----
    Run NSL_SIMULATOR located in the specified executable directory.
    """
    cur_dir = Path(os.getcwd())
    # Changing directory to the executable directory because NSL_SIMULATOR depends on relative paths.
    os.chdir(exe_dir)
    print(os.getcwd())
    command = ["./simulator.exe"]
    if save_config:
        command.append("--write-config")

    print(command)

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # , cwd=os.getcwd()
        bufsize=1,
        universal_newlines=True,
    ) as proc:
        for line in proc.stdout:
            print(line)

    os.chdir(cur_dir)
    print("Program executed")


def load_config(
    phase: Phase,
    src_dir: Path = Path("eser_4/equilibration_conf"),
    *,
    reset_flag: bool = False,
    is_montecarlo: bool = False,
) -> None:
    """
    Load configuration files for a given phase.

    Parameters
    ----------
    phase : (Phase)
        The phase for which the configuration files are loaded.
    src_dir : Path, optional
        The source directory where the configuration files are located. Default is "eser_4/equilibration_conf".
    reset_flag : bool, optional
        Flag indicating whether to reset the configuration files. Default is False.
    is_montecarlo : bool, optional
        Flag indicating whether the simulation is a Monte Carlo simulation. Default is False.

    Notes
    -----
    This function copies the necessary configuration files from the source directory to the destination directory.
    The destination directory is "NSL_SIMULATOR/INPUT".
    The function also prints messages indicating the progress of the configuration loading process.
    If the reset_flag is True, the function copies specific files based on the phase and simulation type.
    If the reset_flag is False, the function copies the default configuration files.

    Examples
    --------
    >>> load_config(Phase.SOLID, reset_flag=True)
    Actually working on phase : SOLID
    Configuration for phase a loaded
    Properties loaded
    """

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
        if not is_montecarlo:
            shutil.copy(
                src_dir / f"velocities_{phase_name}.out",
                dest_dir / "CONFIG/velocities.in",
            )
    else:
        shutil.copy(src_dir / "config.xyz", dest_dir / "CONFIG/config.xyz")
