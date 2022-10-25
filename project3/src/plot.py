import numpy as np 
import pyarma as pa 
import matplotlib.pyplot as plt 
import os
import seaborn as sns
# See github repo for more information on how to run and reproduce results:
# https://github.com/Fslippe/FYS4150/tree/main/project3

plt.rcParams.update({"lines.linewidth": 2})
plt.rcParams.update({"font.size": 10})
sns.set_style("whitegrid")


def run(N, T, n, interaction, method, time_dependency="false", f=0.4, omega=0.685, compile=False):
    """
    Compile and run c++ code for different parameters
    - N             Number of timepoints
    - T             Time to compute
    - n             Number of particles
    - interaction   Run with particle interactions (true or false)
    - method        Method to run (Analytic, Euler or RK4)
    - compile (opt) Compile and linking code before running (True or False)
    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ -std=c++11 main.cpp penningtrap.cpp particle.cpp -o main -larmadillo -O2"),
    print("Running")
    run_cpp = os.system("./main %s %s %s %s %s %s %s %s" %(N, T, n, interaction, method, time_dependency, f, omega))   

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


def run_frequency_scan(N, T, n, interaction, omega_min, omega_max, omega_step, compile=False):
    """

    """
    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system("g++ -std=c++11 frequency_scan.cpp penningtrap.cpp particle.cpp -o frequency_scan -larmadillo -O2"),
    run_cpp = os.system("./frequency_scan %s %s %s %s %s %s %s" %(N, T, n, interaction, omega_min, omega_max, omega_step))   

  
    data = pa.mat()
    data.load("data/frac_p_left.dat")
    r = np.array(data)
    return r


def plot_X(X, x_axis, y_axis, linestyle="-", label="",show=False):
    """
    Function to plot either position or velocities of particles
    Takes in:
    - X                 3D array containing ([x,y,z], particle, r(t))
    - x_axis            Axis of X to plot on x axis (0, 1 or 2)
    - y_axis            Axis of X to plot on y axis (0, 1 or 2)
    - linestyle (opt)   Linestyle for plot default "-"
    - label (opt)       label for plot when using one particle
    - show  (opt)       Show plot if True, default False
    """
   
    if np.size(X, axis=1) == 1:
        plt.plot(X[x_axis, 0, :], X[y_axis, 0, :], linestyle=linestyle, label="%s" %(label))
    else:
        for i in range(np.size(X, axis=1)):
            plt.plot(X[x_axis, i, :], X[y_axis, i, :], zorder=-0, linestyle=linestyle, label="Particle %i %s" %(i+1, label))
        plt.scatter(X[x_axis, :, 0], X[y_axis, :, 0], s=50, marker="x", zorder=1, color="k", label="%s" %("Start"))
        plt.scatter(X[x_axis, :, -1], X[y_axis, :, -1], s=50, marker="o", zorder=1, color="r", label="%s" %("Finish"))
        

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
    plt.axis("equal")

    plt.xlabel(r"%s ($\mu$s)" %(x_label))
    plt.ylabel(r"%s ($\mu$m)" %(y_label))
    plt.legend()

    if show:
        plt.show()

def plot_Xt(X, t, y_axis, linestyle="-", label="", linewidth=2, ylim=False, first_index =0):
    """
    Function to plot either position or velocities of particles as a function of time
    Takes in:
    - X                 3D array containing ([x,y,z], particle, r(t))
    - t                 time array matching length of X
    - y_axis            Axis of X to plot on y axis (0, 1 or 2 or "norm" for length of vector)
    - linestyle (opt)   Linestyle for plot default "-"
    - label (opt)       label for plot when using one particle
    - ylim (opt)        if True adds y limits [-600,600] to x and y plots
    - first_index (opt) First index to plot for
    """
    
    if y_axis == "norm":
        X = np.linalg.norm(X, axis=0)
        
        plt.plot(t[first_index:], X[0, first_index:], linestyle=linestyle, linewidth=linewidth, label="%s" %(label))
        plt.ylim(0, 600)

    else:
        if np.size(X, axis=1) == 1: 
            plt.plot(t[first_index:], X[y_axis, 0, first_index:], linestyle=linestyle, label="%s" %(label))
            if (y_axis==0) or (y_axis==1) and ylim==True:
                plt.ylim(-600, 600)
        else:
            for i in range(np.size(X, axis=1)):
                plt.plot(t, X[y_axis, i, :], linestyle=linestyle, label = "P_%i %s" %(i+1, label))

    if y_axis == 0:
        y_label = "x"
    elif y_axis == 1:
        y_label = "y"
    elif y_axis == 2:
        y_label = "z" 
    elif y_axis == "norm":
        y_label = "|r|"

    plt.xlabel(r"$t$ ($\mu$s)")
    plt.ylabel(r"%s ($\mu$m)" %(y_label))
    plt.legend()

def plot_phase_space(r, N, v, axis, save=False):
    """
    Function to plot either phase space of particles
    Takes in:
    - r                 3D array containing ([x,y,z], particle, r(t))
    - v                 3D array containing ([x,y,z], particle, r(t))
    - axis              Axis to plot (0, 1 or 2) equaling (x, y or z)
    - save (opt)        Savename of plot, default False 
    """
   

    for i in range(np.size(r, axis=1)):
        plt.plot(r[axis,i,:], v[axis,i,:], zorder=-0, label="Particle %i" %(i+1))
    plt.scatter(r[axis, :,0], v[axis, :,0], s=50, zorder=1, color="k", marker="x", label="Start")
    plt.scatter(r[axis, :,-1], v[axis, :,-1], s=50, zorder=1, color="r", label="Finish")
    plt.axis("equal")
    
    if axis==0:
        plt.xlabel(r"$x$ ($\mu$m)")
        plt.ylabel(r"$v_x$ [$\mu$m/$\mu$s]")
    elif axis==1:
        plt.xlabel(r"$y$ ($\mu$m)")
        plt.ylabel(r"$v_y$ [$\mu$m/$\mu$s]") 
    elif axis==2:
        plt.xlabel(r"$z$ ($\mu$m)")
        plt.ylabel(r"$v_z$ [$\mu$m/$\mu$s]")   
    plt.legend(loc="upper right")

    if save != False:
        plt.savefig("../figures/%s_N%i.pdf" %(save, N), dpi=300, bbox_inches="tight")
    plt.show()

def plot_3D(r, N, save=False):
    """
    Function to plot 3D position of particles
    Takes in:
    - r                 3D array containing ([x,y,z], particle, r(t))
    - save (opt)        Savename of plot, default False 
    """

    plt.rcParams.update({"font.size": 10})
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    for i in range(np.size(r, axis=1)):
        ax.plot3D(r[0,i,:], r[1,i,:], r[2,i,:], zorder=-0, label="Particle %i" %(i+1))
    ax.scatter3D(r[0,:,0], r[1,:,0], r[2,:,0], zorder=1, color="k", marker="x",s=50, label="Start")
    ax.scatter3D(r[0,:,-1], r[1,:,-1], r[2,:,-1], zorder=1, color="r", s=50, label="Finish")

    ax.set_xlabel(r"x ($\mu$m)")
    ax.set_ylabel(r"y ($\mu$m)")
    ax.set_zlabel(r"z ($\mu$m)")
    plt.legend()
    if save != False:
        plt.savefig("../figures/%s_N%i.pdf" %(save, N), dpi=300, bbox_inches="tight")

    plt.show()

def compare_analytic(N, T, x_axis, y_axis, Euler=True, save=False):
    """
    Function to compare numerical methods to analytic solution.
    Only works for a single particle without particle interaction
    Plots both a 2D position plot and a time-position plot
    - N                 Number of timepoints
    - T                 Time to compute
    - x_axis            Axis of position to plot on x-axis in 2D position plot
    - y_axis            Axis of position to plot on y-axis in 2D position plot and time dependency plot
    - Euler (opt)       if False only compares RK4 to analytic
    - save (opt)        Savename of plot, default False 
    """
   
    
    r_A = run(N, T, 1, "false", "Analytic")
    r_E, v_E = run(N, T, 1, "false", "Euler")
    r_RK, v_RK = run(N, T, 1, "false", "RK4")
    plt.title("Compare Numerical and Analytic solutions")

    plt.plot(r_A[:,x_axis], r_A[:,y_axis], "k-", label="Analytic")
    plot_X(r_RK, x_axis, y_axis, "--", label="RK4")
    if Euler:
        plot_X(r_E, x_axis, y_axis, "dotted", label="Euler")

    plt.axis("equal")
    if save != False:
        plt.savefig("../figures/%s_axis_%i_%i_N%i.pdf" %(save, x_axis, y_axis, N), dpi=300, bbox_inches="tight")
    plt.show()

    t = np.linspace(0, T, N+1)
    plt.plot(t, r_A[:,y_axis], "k-", label="Analytic")
    plot_Xt(r_RK, t, y_axis, linestyle="--", label="RK4")
    if Euler:
        plot_Xt(r_E, t, y_axis, linestyle="dotted", label="Euler")
    if save != False:
        plt.savefig("../figures/%s_t_axis_%i_N%i.pdf" %(save, y_axis, N), dpi=300, bbox_inches="tight")

    plt.show()

def compare_error(N, T, method_in, save=False, norm=True):
    """
    Function to compare relative error for different number of timepoints
    Plotting either error for each dimension or using total length of r vector
    Only works for a single particle without particle interaction
    - N                 Array containing different number of timepoints to test for
    - T                 Time to compute
    - method_in         Method to use (Euler or RK4)
    - save (opt)        Savename of plot, default False 
    - norm (opt)        Using length of vector (True or False) default True
    """
   
    delta_max = np.zeros(len(N))
    for i in range(len(N)):
        t = np.linspace(0, T, N[i]+1)
        r_a = run(N[i], T, 1, interaction="false", method="Analytic")
        r, v = run(N[i], T, 1, interaction="false", method= method_in)
        if norm:
            delta = np.linalg.norm(r[:,0,:] - r_a.T, axis=0)
            rel_error = delta / np.linalg.norm(r_a, axis=1)
            delta_max[i] = np.max(delta)
            
            plt.plot(t, rel_error, label="N=%i" %(N[i]))
            plt.legend()
            plt.yscale("log")
            plt.ylabel("Relative Error")
            plt.xlabel(r"Time ($\mu$s)")

        else:
            rel_error = (r - r_a.T) / r_a
            plt.plot(t, rel_error[0,0,:], label="x")
            plt.plot(t, rel_error[1,0,:], label="y")
            plt.plot(t, rel_error[2,0,:], label="z")
            plt.legend()
            if save != False:
                plt.savefig("../figures/%s_%s_%i.pdf" %(save, method_in, N[i]), dpi=300, bbox_inches="tight")
            plt.show()
            
    r_err = 0

    for k in range(1, len(N)):
        print(np.log(delta_max[k]/delta_max[k-1]) / np.log(50/N[k]/(50/N[k-1])))
        r_err += np.log(delta_max[k]/delta_max[k-1]) / np.log(50/N[k]/(50/N[k-1]))

    r_err = r_err/(len(N)-1)

    plt.title(r"%s, Convergence rate $r_{err}=$%.3f" %(method_in, r_err))

    if norm == True and save != False:
        plt.savefig("../figures/%s_%s_norm.pdf" %(save, method_in), dpi=300, bbox_inches="tight")
    plt.show()

def plot_p_fraction_frequency(r, label="",save = False):
    """
    Function to plot fraction of particles left in Penning trap for a range of frequencies and different amplitudes f.
    Takes in:
    - r                 2D array containing 
                        ([frequencies], [fractions for f=0.1], [fractions for f=0.4[fractions for f=0.7])
    - save (opt)        Savename of plot, default False 
    """
    if np.size(r, axis=1) == 2:
        plt.plot(r[:,0], r[:,1], label = label)
    else:

        plt.plot(r[:,0], r[:,1], label = "$f = 0.1$")
        plt.plot(r[:,0], r[:,2], label = "$f = 0.4$")
        plt.plot(r[:,0], r[:,3], label = "$f = 0.7$")
    plt.ylabel("Fraction of Particles left in trap")
    plt.xlabel("Frequency $\omega_v$ [$MHz$]") 
    plt.legend()

    if save != False:
        plt.savefig("../figures/%s.pdf" %(save), dpi=300, bbox_inches="tight")


def main():
    """Compile c++ file main.cpp"""
    #r, v = run(1, 1, 1, "false", "RK4", compile=True) 
    """Compile c++ file frequency_scan.cpp"""
    #run_frequency_scan(1, 1, 1, "false", 0, 1, 1, compile=True)

    N = 5000 #Timesteps
    T = 500 #Time
    n = 1 #particles
    f = 0.1 # Amplitude
    omega = 2.2 # Frequency
    N_array = np.array([4000, 8000, 16000, 32000]) # To use for relative error plot
    
    plot_compare_analytic = False # Compares Euler, RK4 and Euler for chosen parameters above and saviing figures
    phase_space_and_position = False # Phase space plots + 3D plots for two particles with and without interaction 
    compare_error_plot = False # Relative Error plots with convergence rate
    particle_escape = True #
    compare_RK4_analytic = False  

    """Uses N=5000, T=500, n=100"""
    wide_freq_scan = False
    narrow_freq_scan = False
    narrow_freq_scan_interaction = False 

    """Comparing to analytic solutions"""
    if plot_compare_analytic:
        compare_analytic(N, T, x_axis=0, y_axis=1, save="compare_analytic")
        compare_analytic(N, T, x_axis=0, y_axis=2, save="compare_analytic")

    if phase_space_and_position:
        """Phase space plot without interaction"""
        method = "RK4"
        r, v = run(N, T, 2, interaction="false", method=method)
        plot_phase_space(r, N, v, 0, save="phase_space_x_%s" %(method))
        plot_phase_space(r, N, v, 2, save="phase_space_z_%s" %(method))

        """3D plot"""
        plot_X(r, 0, 1)
        plt.savefig("../figures/2p_N%i_%s_xy.pdf" %(N, method), dpi=300, bbox_inches="tight")
        plt.show()
        plot_3D(r, N, save="3D_2_particles_%s" %(method))

        """Phase space plot with interaction"""
        r, v = run(N, T, 2, interaction="true", method=method)
        plot_phase_space(r, N, v, 0, save="phase_space_x_interaction_%s" %(method))
        plot_phase_space(r, N, v, 2, save="phase_space_z_interaction_%s" %(method))

        """3D plot"""
        plot_X(r, 0, 1)
        plt.savefig("../figures/2p_N%i_%s_xy_interaction.pdf" %(N, method), dpi=300, bbox_inches="tight")
        
        plt.show()

        plot_3D(r, N, save="3D_2_particles_%s_interaction" %(method))


    """Compare Error for different N"""
    if compare_error_plot:

        compare_error(N_array, T, "Euler", save="relative_error", norm=True)
        compare_error(N_array, T, "RK4", save="relative_error", norm=True)

    if particle_escape:
        r_t, v = run(N, T, n, interaction="true", method="RK4", time_dependency="true", f=f, omega=omega, compile=False)
        r, v = run(N, T, n, interaction="false", method="RK4", time_dependency="false", f=f, omega=omega, compile=False)
        t = np.linspace(0, T, N+1)
        y_axis = 2
        plot_Xt(r, t, y_axis, linestyle="-",ylim=True, first_index=3000, label="Constant $V_0$")
        plot_Xt(r_t, t, y_axis, linestyle="--",ylim=True, first_index=3000, label="Time dependent $V_0$ $f=$%.1f $\omega_V=$%.2f" %(f, omega))
        plt.savefig("../figures/%s.pdf" %("ressonance_p1_f%.2f_omega_%.2f_z_%i" %(f, omega, T)), dpi=300, bbox_inches="tight")
        plt.show()
        
        y_axis = 1
        plot_Xt(r, t, y_axis, linestyle="-",ylim=True, first_index=3000, label="Constant $V_0$")
        plot_Xt(r_t, t, y_axis, linestyle="--",ylim=True, first_index=3000, label="Time dependent $V_0$ $f=$%.1f $\omega_V=$%.2f" %(f, omega))
        plt.savefig("../figures/%s.pdf" %("ressonance_p1_f%.2f_omega_%.2f_y_%i" %(f, omega, T)), dpi=300, bbox_inches="tight")
        plt.show()
        
        y_axis = 0
        plot_Xt(r, t, y_axis, linestyle="-",ylim=True, first_index=3000, label="Constant $V_0$")
        plot_Xt(r_t, t, y_axis, linestyle="--",ylim=True, first_index=3000, label="Time dependent $V_0$ $f=$%.1f $\omega_V=$%.2f" %(f, omega))
        plt.savefig("../figures/%s.pdf" %("ressonance_p1_f%.2f_omega_%.2f_x_%i" %(f, omega, T)), dpi=300, bbox_inches="tight")
        plt.show()

        y_axis = "norm"
        plot_Xt(r, t, y_axis, linestyle="-",ylim=True, first_index=3000, label="Constant $V_0$")
        plot_Xt(r_t, t, y_axis, linestyle="-",linewidth=0.5, ylim=True, first_index=3000, label="Time dependent $V_0$ $f=$%.1f $\omega_V=$%.2f" %(f, omega))
        plt.savefig("../figures/%s.pdf" %("ressonance_p1_norm_f%.2f_omega_%.2f_x_%i" %(f, omega, T)), dpi=300, bbox_inches="tight")
        plt.show()

    if compare_RK4_analytic:
        compare_analytic(N, T, x_axis=0, y_axis=1, Euler=False, save="analytic_RK4")
        compare_analytic(N, T, x_axis=0, y_axis=2, Euler=False, save="analytic_RK4")

    N = 5000
    T = 500
    n = 100

    """Wide Frequency scan"""
    if wide_freq_scan:
        r = run_frequency_scan(N, T, n, interaction="false",
            omega_min=0.2, omega_max=2.5, omega_step=0.02, compile=True)
        plot_p_fraction_frequency(r, save="wide_freq_p_left_N_%s" %(N))
        plt.show()


    """Narrow Frequency scan without interaction"""
    if narrow_freq_scan:
        r = run_frequency_scan(N, T, n, interaction="false",
            omega_min=2.05, omega_max=2.35, omega_step=0.001, compile=True)
        plot_p_fraction_frequency(r, save="narrow_freq_p_left_N_%s_zoom" %(N)) 

        """
        if plotting comparison plot change parameters n_f and f in 
        PenningTrap::parcticles_left_for_omega_v in
        penningtrap.cpp #-mark plot_p_fraction_frequency above and unmark below
        """
        #plot_p_fraction_frequency(r, label="without interaction") 
        plt.show()



    """Narrow Frequency scan with interaction"""
    if narrow_freq_scan_interaction:
        r = run_frequency_scan(N, T, n, interaction="true",
            omega_min=2.05, omega_max=2.35, omega_step=0.001, compile=False)
        plot_p_fraction_frequency(r, save="narrow_freq_p_left_N_%s_interaction_zoom" %(N))

        """
        if plotting comparison plot change parameters n_f and f in 
        PenningTrap::parcticles_left_for_omega_v in
        penningtrap.cpp #-mark plot_p_fraction_frequency above and unmark below
        """
        #plot_p_fraction_frequency(r, label="with interation") 
        #plt.savefig("../figures/%s.pdf" %("narrow_freq_p_left_N_%s_zoom" %(N)), dpi=300, bbox_inches="tight")

        plt.show()

if __name__ == "__main__":
    main()
