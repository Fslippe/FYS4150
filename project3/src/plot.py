import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 
import os
plt.rcParams.update({"lines.linewidth": 2})
plt.rcParams.update({"font.size": 14})


def run(N, T, n, interaction, method , compile=False):
    """

    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ main.cpp penningtrap.cpp particle.cpp -o main -larmadillo"),
    run_cpp = os.system("./main %s %s %s %s %s" %(N, T, n, interaction, method))   

    if method == "Analytic":
        data = pa.mat()
        data.load("data/r_a.dat")
        r = np.array(data)
        return r
    else:
        data = pa.cube()
        data.load("data/r_%s.dat" %(method))
        r = np.array(data)
        data.load("data/v_%s.dat" %(method))
        v = np.array(data)
        return r, v



def plot_X(X, x_axis, y_axis, label="",show=False):
    """
    Function to plot either position or velocities of particles
    Takes in:
    - X         3D array ()
    - x_axis    Axis of X to plot on x axis (0, 1 or 2)
    - y_axis    Axis of X to plot on y axis (0, 1 or 2)
    - show      Show plot if True, default False
    """

    if np.size(X, axis=1) == 1:
        plt.plot(X[x_axis, 0, :], X[y_axis, 0, :], label="%s" %(label))
    else:
        for i in range(np.size(X, axis=1)):
            plt.plot(X[x_axis, i, :], X[y_axis, i, :], label="Particle %i %s" %(label))
        
    if x_axis == 0:
        x_label = "x"
    elif x_axis == 1:
        x_label = "y"
    elif x_axis == 2:
        x_label = "z"

    if y_axis == 0:
        y_label = "x"
    elif y_axis == 1:
        y_label = "y"
    elif y_axis == 2:
        y_label = "z" 

    plt.xlabel(r"%s [$\mu$s]" %(x_label))
    plt.ylabel(r"%s [$\mu$m]" %(y_label))
    plt.legend()
    if show:
        plt.show()

def plot_Xt(X, t, y_axis, label=""):
    if np.size(X, axis=1) == 1:
        plt.plot(t, X[y_axis, 0, :], label=label)
    else:
            
        for i in range(np.size(X, axis=1)):
            plt.plot(t, X[y_axis, i, :], label = "Particle %i %s" %(i+1, label))

    if y_axis == 0:
        y_label = "x"
    elif y_axis == 1:
        y_label = "y"
    elif y_axis == 2:
        y_label = "z" 

    plt.xlabel(r"$t$ [$\mu$s]")
    plt.ylabel(r"%s [$\mu$m]" %(y_label))
    plt.legend()

def plot_phase_space(r, v, axis, save=False):
    plt.title("Phase space for axis %i" %(axis))
    for i in range(np.size(r, axis=1)):
        plt.plot(r[axis,i,:], v[axis,i,:], label="Particle %i" %(i+1))
    plt.scatter(r[axis, :,0], v[axis, :,0], color="r", label="Start")
    plt.scatter(r[axis, :,-1], v[axis, :,-1], color="b", label="Finish")
    
    if axis==0:
        plt.xlabel(r"$x$ [$\mu$m]")
        plt.ylabel(r"$v_x$ [$\mu$m/$\mu$s]")
    elif axis==1:
        plt.xlabel(r"$y$ [$\mu$m]")
        plt.ylabel(r"$v_y$ [$\mu$m/$\mu$s]") 
    elif axis==2:
        plt.xlabel(r"$z$ [$\mu$m]")
        plt.ylabel(r"$v_z$ [$\mu$m/$\mu$s]")   
    plt.legend()
    plt.axis("equal")

    if save != False:
        plt.savefig("../figures/%s.pdf" %(save), dpi=300, bbox_inches="tight")
    plt.show()

def plot_3D(r, save=False):
    plt.rcParams.update({"font.size": 10})
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    for i in range(np.size(r, axis=1)):
        ax.plot3D(r[0,i,:], r[1,i,:], r[2,i,:], label="Particle %i" %(i+1))
    ax.scatter3D(r[0,:,0], r[1,:,0], r[2,:,0], color="r",s=20, label="Start")
    ax.scatter3D(r[0,:,-1], r[1,:,-1], r[2,:,-1], color="b", s=20, label="Finish")

    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")
    ax.set_zlabel(r"z [$\mu$m]")
    plt.legend()
    if save != False:
        plt.savefig("../figures/%s.pdf" %(save), dpi=300, bbox_inches="tight")

    plt.show()

def compare_analytic(N, T, x_axis, y_axis, save=False):
    r_A = run(N, T, 1, "false", "Analytic")
    r_E, v_E = run(N, T, 1, "false", "Euler")
    r_RK, v_RK = run(N, T, 1, "false", "RK4")
    plt.title("Compare Numerical and Analytic solutions")

    plt.plot(r_A[:,x_axis], r_A[:,y_axis], label="Analytic")
    plot_X(r_RK, x_axis, y_axis, label="RK4")
    plot_X(r_E, x_axis, y_axis, label="Euler")
    plt.axis("equal")
    if save != False:
        plt.savefig("../figures/%s_axis_%i_%i.pdf" %(x_axis, y_axis, save), dpi=300, bbox_inches="tight")
    plt.show()

    t = np.linspace(0, T, N)
    plt.plot(t, r_A[:,y_axis], label="Analytic")
    plot_Xt(r_RK, t, y_axis, label="RK4")
    plot_Xt(r_E, t, y_axis, label="Euler")
    if save != False:
        plt.savefig("../figures/%s_t_axis_%i.pdf" %(x_axis, y_axis, save), dpi=300, bbox_inches="tight")
    plt.show()

def compare_relative_error(N, T, method, save=False, norm=True):
    
    for i in N:
        print(i)
        t = np.linspace(0, T, i)
        r_a = run(i, T, 1, interaction="false", method="Analytic")
        r, v = run(i, T, 1, interaction="false", method= method)
        
        if norm:
            print(r[:,0,:] - r_a.T)
            rel_error = np.linalg.norm(r[:,0,:] - r_a.T, axis=0)  / np.linalg.norm(r_a.T)
            
            plt.plot(t, rel_error, label="N=%i" %(i))
            plt.legend()
            plt.yscale("log")

        else:
            rel_error = (r - r_a.T)
            plt.plot(t, rel_error[0,0,:], label="x")
            plt.plot(t, rel_error[1,0,:], label="y")
            plt.plot(t, rel_error[2,0,:], label="z")
            plt.legend()
            if save != False:
                plt.savefig("../figures/%s_%s_%i.pdf" %(save, method, i), dpi=300, bbox_inches="tight")
            plt.show()

    if norm == True and save != False:
        plt.savefig("../figures/%s_%s_%i_norm.pdf" %(save, method, i), dpi=300, bbox_inches="tight")
    plt.show()


def main():
    N = 50000
    T = 50 
    print (T/N)
    #run(1, 1, 1, "false", "Euler", compile=True) # Compile c++ file

    compare_analytic(N, T, x_axis=1, y_axis=1)
    """
    r, v = run(N, T, 2, interaction="false", method="RK4")
    plot_phase_space(r, v, 0, save="phase_space_x_no_interaction")
    plot_3D(r, save="3D_2_particles_no_interaction")

    r, v = run(N, T, 2, interaction="true", method="RK4")
    plot_phase_space(r, v, 0, save="phase_space_x_with_interaction")
    plot_3D(r, save="3D_2_particles_with_interaction")
    """
    N_array = np.array([4000, 8000, 16000, 32000])
    compare_relative_error(N_array, T, "Euler", save="relative_error", norm=True)

if __name__ == "__main__":
    main()
