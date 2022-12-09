# Project 5 - Numerical simulation of the 2+1 dimensional Schrödinger equation 
Authors:
Alessio Canclini and Filip Severin von der Lippe

![](https://github.com/Fslippe/FYS4150/blob/main/project5/figures/animation.gif) ![](https://github.com/Fslippe/FYS4150/blob/main/project5/figures/animation_real.gif))
## Structure of repository
This folder of the repository FYS4150 contains all files used for Project 5 in FYS4150 Autumn 2022.
- figures contains both figures used in the report of this project and other relevant figures produced by the included code.
- latex Contains .tex and .pdf of the report for this project
- code contains all code used to prouduce the reults discussed in the report as well as the figures found in the figures folder

## Code structure and how to run
The programming languages python and c++ have been used. To run the code make sure to use python 3 and c++11, and have the following packages for python which can be installed with pip by using the following commands:

```
pip3 install pyarma
pip3 install matplotlib
pip3 install seaborn
```

Everything is built around one class slit_box with its own .cpp and .hpp files. 
All computations use the main.cpp file which includes the class to compute the results. 
The python file "plot.py" reproduces all figures and results by compiling, linking, and running the necessary c++ code. All this can be done by
```
python3 plot.py
```
It is also possible to run the individual c++ codes by using the following commands:
```
g++ -std=c++11 main.cpp slit_box.cpp -o main -larmadillo -O2
./main h_in dt_in T_in xc_in sigma_x_in px_in yc_in sigma_y_in py_in v0_in savename slits
```
Here we have the chosen parameters given as command line arguments:
- *double* **h_in**                  h constant and position step size      
- *double* **dt_in**                 time step size of simulation
- *double* **T_in **                 time to simulate
- *double* **xc_in**                 x center of initial wave packet
- *double* **sigma_x_in**            x width of initial wave packet
- *double* **px_in**                 wave packet x momentum
- *double* **yc_in**                 y center of initial wave packet
- *double* **sigma_y_in**            y width of initial wave packet
- *double* **py_in**                 wave packet y momentum
- *int* **v0_in**                 potential value inside walls
- *string* **savename**              savename of data files
- *int* **slits**                 number of slits 
