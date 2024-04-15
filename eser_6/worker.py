import subprocess
from pathlib import Path
from os import chdir
from dataclasses import dataclass
import shutil
from time import perf_counter
import numpy as np

# Program Constants

TEMP_MIN: float = 0.5
TEMP_MAX: float = 2.0
TEMP_STEP: int = 50


SIMULATIONS: list[tuple[int, float, float]] = [
    # (2, 1.0, 0.0),
    # (2, 1.0, 0.2),
    (3, 1.0, 0.0),
    (3, 1.0, 0.2),
]

SIMULATION_TYPE = {2: "METRO", 3: "GIBBS"}

WORKING_PATH = Path("NSL_SIMULATOR/SOURCE")
PROGRAM_PATH = Path("./simulator.exe")
INPUT_CONFIG_PATH = Path("../INPUT")
OUTPUT_CONFIG_PATH = Path("../OUTPUT")
OUTPUT_RESULT_PATH = Path("../../eser_6/output")

INPUT_CONFIG_PATH_FILE = INPUT_CONFIG_PATH / Path("input.dat")


@dataclass
class Config:
    simulation_type: tuple[int, float, float]
    restart: int
    temp: float
    n_part: int
    rho: float
    r_cut: float
    delta: float
    n_blocks: int
    n_steps: int

    def __str__(self) -> None:
        return_string = f"""SIMULATION_TYPE        {self.simulation_type[0]}    {self.simulation_type[1]}    {self.simulation_type[2]}
RESTART                {self.restart}
TEMP                   {self.temp}
NPART                  {self.n_part}
RHO                    {self.rho}
R_CUT                  {self.r_cut}
DELTA                  {self.delta}
NBLOCKS                {self.n_blocks}
NSTEPS                 {self.n_steps}

ENDINPUT"""
        return return_string


def check_directories() -> bool:
    can_run = True

    print(f"Actual Directory : {Path().cwd()}")
    if WORKING_PATH.exists():
        chdir(WORKING_PATH)
        print(f"New directory : {Path().cwd()}")
    else:
        print("{WORKING_PATH}: Directory Not Found")
        can_run = False

    if PROGRAM_PATH.exists():
        print(f"Found Executable : {PROGRAM_PATH}")
    else:
        print(f"{PROGRAM_PATH}: File Not Found")
        can_run = False

    if INPUT_CONFIG_PATH.exists():
        print(f"Found Input Config File Directory : {INPUT_CONFIG_PATH}")
    else:
        print(f"{INPUT_CONFIG_PATH}: Directory Not Found")
        can_run = False

    if OUTPUT_CONFIG_PATH.exists():
        print(f"Found Input Config File Directory : {OUTPUT_CONFIG_PATH}")
    else:
        print(f"{OUTPUT_CONFIG_PATH}: Directory Not Found")
        can_run = False

    return can_run


def launch_program() -> subprocess.CompletedProcess:
    print("Launching Program")
    result = subprocess.run(PROGRAM_PATH.absolute())
    print(f"Result Code : {result.returncode}")
    return result


def write_config(conf: Config) -> None:
    with open(INPUT_CONFIG_PATH_FILE, "w") as file:
        file.write(str(conf))


def move_output(conf: Config, input_dir: Path, output_dir: Path) -> None:
    conf_out = Path(
        f"output_temp_{conf.temp:.3f}_{SIMULATION_TYPE[conf.simulation_type[0]]}_ext_field_{int(conf.simulation_type[2] * 10)}"
    )
    for file in input_dir.iterdir():
        if file.is_file():
            (output_dir / conf_out).mkdir(exist_ok=True)
            shutil.copy(file, output_dir / conf_out)


def main():
    if not check_directories():
        exit(-1)

    temperatures_range = np.linspace(TEMP_MIN, TEMP_MAX, num=TEMP_STEP)
    actual_config = Config(
        simulation_type=[2, 1.0, 1.0],
        restart=0,
        temp=1.0,
        n_part=50,
        rho=1.0,
        r_cut=0.0,
        delta=0.0,
        n_blocks=20,
        n_steps=20000,
    )

    start_time = perf_counter()
    counter_success = 0
    for sim in SIMULATIONS:
        actual_config.simulation_type = sim
        for temp in np.nditer(temperatures_range):
            actual_config.temp = temp
            print(
                f"Actual temperature : {actual_config.temp}\n"
                + f"Actual Simulation : {SIMULATION_TYPE[actual_config.simulation_type[0]]}"
            )
            write_config(actual_config)
            res = launch_program()
            move_output(actual_config, OUTPUT_CONFIG_PATH, OUTPUT_RESULT_PATH)
            if res.returncode == 0:
                counter_success += 1
    end_time = perf_counter()

    print(
        f"Process completed, success : {counter_success} / {len(SIMULATIONS) * TEMP_STEP}"
    )
    print(f"Elapsed time : {end_time - start_time} s")

    print("Ending")
    exit(0)


if __name__ == "__main__":
    main()
