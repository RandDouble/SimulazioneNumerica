# LSN 2024, esercizi di Stefano Pilosio

Per alcune esigenze personali il mio PC è dotato di [Clang-17](https://clang.llvm.org) e non ho g++, ho verificato con altre macchine un filo meno portatili che tutto il codice compila anche con g++.
Altra nota abbastanza rilevante è che tutto il codice usa lo standard [C++17](https://www.iso.org/standard/68564.html), in quanto tale non compila soltanto con C++11.
Il codice per quanto possibile ha le funzioni commentate seguendo lo standard di [Doxygen](https://www.doxygen.nl).

Per compilare con g++ è necessario usare [Threading Building Block Library](https://www.threadingbuildingblocks.org)

Usando Clang invece è necessario assicurarsi di star usando libc++ come libreria standard. Se si usa invece libstdc++ invece valgono le stesse considerazioni di g++.

Clang fornito da MacOs manca completamente il supporto a una libreria standard parallela... Magari in un futuro scoprirò la bellezza di openmp e migliorerò il codice di conseguenza.

Pertanto, per compilare tutto il codice basta modificare il makefile nella root directory di questo progetto inserendo il compilatore di propria preferenza e scrivere `make`.

## Tabella Esercizi Svolti

I titoli dell'esercizio sono dei link ai notebook, invece i link negli elementi puntati portano al source file dell'esercizio.

### [Esercitazione 1](LSN_Exercises_01.ipynb)

- [x] [Esercizio 1.1](eser_1/1_1/main_1_1.cpp)
- [x] [Esercizio 1.2](eser_1/1_2/main_1_2.cpp)
- [x] [Esercizio 1.3](eser_1/1_3/main_1_3.cpp)

### [Esercitazione 2](LSN_Exercises_02.ipynb)

- [x] [Esercizio 2.1](eser_2/2_1/main_2_1.cpp)
- [ ] [Esercizio 2.2](eser_2/2_1/main_2_1.cpp)

### [Esercitazione 3](LSN_Exercises_03.ipynb)

- [x] [Esercizio 3.1](eser_3/3_1/main_3_1.cpp)

### [Esercitazione 4](LSN_Exercises_04.ipynb)

- [ ] [Esercizio 4.1](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)
- [ ] [Esercizio 4.2](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)

### [Esercitazione 5](LSN_Exercises_05.ipynb)

- [x] [Esercizio 5.1](eser_5/5_1/main_5_1.cpp)

### [Esercitazione 6](LSN_Exercises_06.ipynb)

- [ ] [Esercizio 6.1](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)

### [Esercitazione 7](LSN_Exercises_07.ipynb)
- [x] [Esercizio 7.1](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)
- [x] [Esercizio 7.2](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)
- [x] [Esercizio 7.3](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)
- [ ] [Esercizio 7.4](NSL_SIMULATOR/SOURCE/NSL_SIMULATOR.cpp)

### [Esercitazione 8](LSN_Exercises_08.ipynb)

### [Esercitazione 9](LSN_Exercises_09.ipynb)
### [Esercitazione 10](LSN_Exercises_10.ipynb)
### [Esercitazione 11](LSN_Exercises_11.ipynb)
### [Esercitazione 12](LSN_Exercises_12.ipynb)

---

Per eventuali problematiche potete contattarmi al seguente indirizzo di posta elettronica: stefano.pilosio@studenti.unimi.it
