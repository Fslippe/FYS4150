# Project 3 - Numerical study of the Penning Trap
Authors:
Alessio Canclini and Filip Severin von der Lippe

## Structure of repository
This folder of the repository FYS4150 contains all files used for Project 3 in FYS4150 Autumn 2022.
- figures contains both figures used in the report of this project and other relevant figures produced by the included code.
- latex Contains .tex and .pdf of the report for this project
- src contains all code used to prouduce the reults discussed in the report aswell as the figures found in the figures folder

## Code structure and how to run
The programming languages python and c++ have been used. To run the code make sure to use python 3 and c++11, and have the following packages for python which can be installed with pip by using the following commands:
```
pip3 install pyarma
pip3 install matplotlib
pip3 install seaborn
```
The python file "plot.py" reproduces all figures and results by compiling, linking, and running the necessary c++ code. All this can be done by
```
python3 plot.py
```
It is also possible to run the individual c++ codes by using the following commands:
```
g++ -std=c++11 main.cpp penningtrap.cpp particle.cpp -o main -larmadillo -O2
./main N T n interaction method time_dependency f omega
g++ -std=c++11 frequency_scan.cpp penningtrap.cpp particle.cpp -o frequency_scan -larmadillo -O2
./frequency_scan N T n interaction omega_min omega_max omega_step
```
Here we have the chosen parameters given as command line arguments:
- *int* **N:**                  Number of timesteps
- *double* T:**               Time
- *int* **n:**                  Number of particles in the Penning Trap 
- *bool* **interaction:**       true to run with particle interactions, false otherwise
- *string* **method:**       Method to run, Analyic, Euler or RK4
- *bool* **time_dependency:**   true to run with time dependent $V_0$, false otherwise
- *double* **f:**  Amplitude of time dependent $V_0$
- *double* **omega:** Frequency of time dependent $V_0$
- *double* **omega_min:** Smallest frequency for $V_0$ when performing frequency scan
- *double* **omega_max:** Largest frequency for $V_0$ when performing frequency scan
- *double* **omega_step:** Frequency stepsize for $V_0$ when performing frequency scan

