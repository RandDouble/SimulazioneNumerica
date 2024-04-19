import sqlite3
import numpy as np
from pathlib import Path
import os


def initialize_db(connection: sqlite3.Connection, db_exist : bool):
    cursor = connection.cursor()
    if not db_exist:
        cursor.execute(
            """CREATE TABLE SIMULAZIONI
            (
                ROWID INTEGER PRIMARY KEY AUTOINCREMENT,
                TEMPERATURE FLOAT,
                TIPOLOGIA INTEGER,
                EXTERNAL_FIELD FLOAT,
                COUPLING FLOAT
            )"""
        )

        cursor.execute(
            """CREATE TABLE ACCEPTANCE
            (
            ROWID INTEGER,
            N_BLOCK INTEGER,
            ACCEPTANCE FLOAT,
            FOREIGN KEY ( ROWID )
                REFERENCES SIMULAZIONI (ROWID) ON DELETE CASCADE
            )"""
        )

        cursor.execute(
            """CREATE TABLE MAGNETIZATION 
            (
                ROWID INTEGER,
                BLOCK INTEGER,
                ACTUAL_M FLOAT,
                M_AVE FLOAT,
                ERROR FLOAT,
                FOREIGN KEY ( ROWID )
                    REFERENCES SIMULAZIONI (ROWID) ON DELETE CASCADE

            )"""
        )

        cursor.execute(
            """CREATE TABLE SPECIFIC_HEAT
            (
                ROWID INTEGER,
                BLOCK INTEGER,
                ACTUAL_CV FLOAT,
                CV_AVE FLOAT,
                ERROR FLOAT,
                FOREIGN KEY ( ROWID )
                    REFERENCES SIMULAZIONI (ROWID) ON DELETE CASCADE
            )"""
        )

        cursor.execute(
            """CREATE TABLE  SUSCEPTIBILITY
            (
                ROWID INTEGER,
                BLOCK INTEGER,
                ACTUAL_CHI FLOAT,
                CHI_AVE FLOAT,
                ERROR FLOAT,
                FOREIGN KEY ( ROWID )
                    REFERENCES SIMULAZIONI (ROWID) ON DELETE CASCADE
            )"""
        )

        cursor.execute(
            """CREATE TABLE  TOTAL_ENERGY
            (
                ROWID INTEGER,
                BLOCK INTEGER,
                ACTUAL_TE FLOAT,
                TE_AVE FLOAT,
                ERROR FLOAT,
                FOREIGN KEY ( ROWID )
                    REFERENCES SIMULAZIONI (ROWID) ON DELETE CASCADE
            )"""
        )
        connection.commit()
    
    cursor.execute("UPDATE SQLITE_SEQUENCE SET SEQ = 1 WHERE NAME = 'SIMULAZIONI';")
    cursor.execute("DELETE FROM SIMULAZIONI")
    cursor.execute("DELETE FROM ACCEPTANCE")
    cursor.execute("DELETE FROM MAGNETIZATION")
    cursor.execute("DELETE FROM SPECIFIC_HEAT")
    cursor.execute("DELETE FROM SUSCEPTIBILITY")
    cursor.execute("DELETE FROM TOTAL_ENERGY")
    connection.commit()


def output_handler(connection: sqlite3.Connection, folder: Path):
    cursor = connection.cursor()
    temperature: float = 0.0
    simulation_type: int = 0
    coupling: float = 0.0
    external_field: float = 0.0

    with open(folder / Path("output.dat"), "r") as file_in:
        for line in file_in:
            if line.strip().lower() == "reading input completed!":
                break
            elif "TEMPERATURE" in line.strip().upper():
                temperature = float(line.strip().split()[1])
            elif "SIM_TYPE" in line.strip().upper():
                simulation_type = int(line.strip().split()[1])
                coupling = float(line.strip().split()[2])
                external_field = float(line.strip().split()[3])

            values = ",".join(str(el) for el in [temperature, simulation_type, external_field, coupling])
        cursor.execute(
            f"INSERT INTO SIMULAZIONI (TEMPERATURE, TIPOLOGIA, EXTERNAL_FIELD, COUPLING) VALUES ({values})"
        )
        print(
            f"INSERT INTO SIMULAZIONI (TEMPERATURE, TIPOLOGIA, EXTERNAL_FIELD, COUPLING) VALUES ({values})"
        )

    return cursor.lastrowid


def acceptance_handler(connection: sqlite3.Connection, file: Path, rowid: int) -> None:
    cursor = connection.cursor()
    data = np.loadtxt(file, dtype=np.int32, skiprows=1)
    for row in data:
        cursor.execute(
            f"INSERT INTO ACCEPTANCE (ROWID, N_BLOCK, ACCEPTANCE) VALUES ({rowid}, {row[0]}, {row[1]})"
        )
    connection.commit()


def magnetization_handler(
    connection: sqlite3.Connection, file: Path, rowid: int
) -> None:
    cursor = connection.cursor()
    data = np.loadtxt(file, dtype=np.float64, skiprows=1)
    for row in data:
        s_row = ",".join(row.astype(str).tolist())
        cursor.execute(
            f"INSERT INTO MAGNETIZATION (ROWID, BLOCK, ACTUAL_M, M_AVE, ERROR) VALUES ({rowid}, {s_row})"
        )
    connection.commit()


def specific_heat_handler(
    connection: sqlite3.Connection, file: Path, rowid: int
) -> None:
    cursor = connection.cursor()
    data = np.loadtxt(file, dtype=np.float64, skiprows=1)
    for row in data:
        s_row = ",".join(row.astype(str).tolist())

        cursor.execute(
            f"INSERT INTO SPECIFIC_HEAT (ROWID, BLOCK, ACTUAL_CV, CV_AVE, ERROR) VALUES ({rowid}, {s_row})"
        )
    connection.commit()


def susceptibility_handler(
    connection: sqlite3.Connection, file: Path, rowid: int
) -> None:
    cursor = connection.cursor()
    data = np.loadtxt(file, dtype=np.float64, skiprows=1)
    for row in data:
        s_row = ",".join(row.astype(str).tolist())
        cursor.execute(
            f"INSERT INTO SUSCEPTIBILITY (ROWID, BLOCK, ACTUAL_CHI, CHI_AVE, ERROR) VALUES ({rowid}, {s_row})"
        )
    connection.commit()


def total_energy_handler(
    connection: sqlite3.Connection, file: Path, rowid: int
) -> None:
    cursor = connection.cursor()
    data = np.loadtxt(file, dtype=np.float64, skiprows=1)
    for row in data:
        s_row = ",".join(row.astype(str).tolist())
        cursor.execute(
            f"INSERT INTO TOTAL_ENERGY (ROWID, BLOCK, ACTUAL_TE, TE_AVE, ERROR) VALUES ({rowid}, {s_row})"
        )
    connection.commit()


def navigator(connection: sqlite3.Connection):
    for folder in Path("output").iterdir():
        rowid: int = output_handler(connection, folder)

        for file in folder.iterdir():
            match file.name:
                case "acceptance.dat":
                    acceptance_handler(connection=connection, file=file, rowid=rowid)
                case "magnetization.dat":
                    magnetization_handler(connection=connection, file=file, rowid=rowid)
                case "specific_heat.dat":
                    specific_heat_handler(connection=connection, file=file, rowid=rowid)
                case "susceptibility.dat":
                    susceptibility_handler(
                        connection=connection, file=file, rowid=rowid
                    )
                case "total_energy.dat":
                    total_energy_handler(connection=connection, file=file, rowid=rowid)
                case _:
                    continue


def main():
    os.chdir("eser_6")
    db_exist = False
    if Path('output.db').exists():
        db_exist = True
    connection = sqlite3.connect("output.db")
    initialize_db(connection=connection, db_exist = db_exist)
    navigator(connection=connection)
    connection.close()


if __name__ == "__main__":
    main()
