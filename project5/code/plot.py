import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
from matplotlib.animation import FuncAnimation
from functools import partial
import matplotlib
plt.rcParams.update({"font.size": 12.5})
sns.set_style("darkgrid")


def run(h_in, dt_in, T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in, savename, slits, compile=False):
    """
    Compile and run c++ code for different parameters
    - h_in                  h constant and position step size      
    - dt_in                 time step size of simulation
    - T_in                  time to simulate
    - xc_in                 x center of initial wave packet
    - sigma_x_in            x width of initial wave packet
    - px_in                 wave packet x momentum
    - yc_in                 y center of initial wave packet
    - sigma_y_in            y width of initial wave packet
    - py_in                 wave packet y momentum
    - v0_in                 potential value inside walls
    - savename              savename of data files
    - slits                 number of slits 
    - compile=False         if True compile c++ file 
    """

    if compile == True:
        print("Compiling...\n")
        compile_cpp = os.system(
            "g++ -std=c++11 main.cpp slit_box.cpp -o main -larmadillo -O2"),
    start = time.time()
    run_cpp = os.system("./main %s %s %s %s %s %s %s %s %s %s %s %s" % (h_in, dt_in,
                        T_in, xc_in, sigma_x_in, px_in, yc_in, sigma_y_in, py_in, v0_in, slits, savename))
    print("time: ", time.time()-start)


def animation_plot(data, h, dt):
    sns.set_style("white")

    # Set up a 2D xy grid
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    t_points = np.arange(0, 1+dt, dt)

    # Some settings
    fontsize = 12
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data[:, 0, :]))
    print(np.shape(data))
    # Plot the first
    img = ax.imshow(data[:, 0, :], extent=[x_min, x_max,
                    y_min, y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel(r"$x$", fontsize=fontsize)
    plt.ylabel(r"$y$", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(r"Re$(u_i)$", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, partial(animation, img=img, data=data, t_min=t_min, dt=dt,
                         time_txt=time_txt), interval=30, frames=np.arange(0, len(data[0, :, 0]), 2),  repeat=True, blit=0)

    # Run the animation!
    plt.show()

    # # Save the animation
    anim.save('../figures/animation_real.gif',
              writer="writegif", bitrate=-1, fps=30)

    sns.set_style("darkgrid")
    plt.rcParams.update({"font.size": 12.5})

    # Function that takes care of updating the z data and other things for each frame


def animation(i, img, data, t_min, dt, time_txt):
    # Normalize the colour scale to the current frame?

    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data[:, 0, :]))
    img.set_norm(norm)

    # Update z data
    img.set_data(data[:, i, :])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img


def plot_p_diff(t, data, save):
    """
    Plot time development of p difference to 1
    args
    - t                 time array matching data 
    - data              data containing p_ij(t)
    - save              savename of plot 
    """
    p_sum = np.zeros(len(t))
    for i in range(len(t)):
        p_sum[i] = np.sum(data[:, i, :])

    plt.plot(t, p_sum-1)
    plt.ylabel(r"$\sum_{i,j} p_{ij}^{t} - 1$")
    plt.xlabel(r"$t$")
    plt.xticks(rotation=45)
    plt.savefig("../figures/%s.pdf" % (save), dpi=300, bbox_inches="tight")
    plt.show()


def plot_time(t, dt, data, V, save, barlabel):
    """
    Plot snapshot of data for any given number of t 
    args:
    - t                 list of t to plot
    - dt                dt used in simulation
    - data              [x, t, y] dataset 
    - V                 V potential used when simulating data
    - save              savename of plot
    - barlabel          label of colorbar
    """
    xy = np.linspace(0, 1, len(data[:, 0, 0]))
    x, y = np.meshgrid(xy, xy)

    for t_i in t:
        plt.figure()
        plt.title(r"$t=$%.3f" % (t_i))
        cnt = plt.contourf(
            x, y, data[:, int(t_i/dt), :], levels=250, cmap="viridis")

        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")

        for c in cnt.collections:
            c.set_edgecolor("face")
        plt.colorbar(label=r"%s" % (barlabel))
        plt.contourf(x, y, V.T, levels=np.linspace(1, 1e10, 2),
                     vmin=1e10-1, vmax=1e10,  corner_mask=True, cmap="gray")
        plt.savefig("../figures/%s_%.3f.pdf" %
                    (save, t_i), dpi=300, bbox_inches="tight")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")

    plt.show()


def plot_at_x(data, t, dt, x, save, title):
    """
    plot probability against y for given t and x
    args:
    - data              dataset containing p_ij 
    - t                 time to plot
    - dt                dt used in simulation
    - x                 x position to plot p
    - save              savename of plot
    - title             title of plot
    """
    x_idx = int(x*len(data[:, 0, 0]))
    p_sum = np.sum(data[:, int(t/dt), x_idx])
    data_normalized = data[:, int(t/dt), x_idx] / p_sum
    plt.title(title)
    plt.plot(np.linspace(0, 1, len(data[:, 0, 0])), data_normalized)
    plt.xlabel(r"$y$")
    plt.ylabel(r"$p(y|x=0.8;t=0.002)$")
    plt.savefig("../figures/%s_%.3f_x_%.1f.pdf" %
                (save, t, x), dpi=300, bbox_inches="tight")
    plt.show()


def main():

    runcpp = True  # RUN c++ file if True
    dt = 2.5*10**(-5)
    h = 0.005

    if runcpp:
        run(h_in=h,
            dt_in=dt,
            T_in=0.008,  # 0.002
            xc_in=0.25,
            sigma_x_in=0.05,
            px_in=200,
            yc_in=0.5,
            sigma_y_in=0.1,  # 0.05,  # 0.2,
            py_in=0,
            v0_in=1e10,
            slits=2,  # 3, #1, #0,
            savename="double_slit_7",
            compile=True)

    """Initializing datafiles to load"""
    DS = pa.cube()
    DS1 = pa.cube()
    DS2 = pa.cube()
    DS2_long = pa.cube()
    DS2_long_real = pa.cube()
    DS3 = pa.cube()
    DS2_real = pa.cube()
    DS2_imag = pa.cube()
    V_1 = pa.mat()
    V_2 = pa.mat()
    V_3 = pa.mat()

    """Loading data"""
    t = pa.mat()
    t.load("data/t8.dat")
    DS.load("data/no_slit_7.dat")
    DS1.load("data/single_slit_9.dat")
    DS2.load("data/double_slit_8.dat")
    DS2_long.load("data/double_slit_7.dat")
    DS2_long_real.load("data/double_slit_7_real.dat")

    DS3.load("data/tripple_slit_9.dat")

    # V for different number of slits
    V_1.load("data/V_1.dat")
    V_2.load("data/V_2.dat")
    V_3.load("data/V_3.dat")
    DS2_real.load("data/double_slit_8_real.dat")
    DS2_imag.load("data/double_slit_8_imag.dat")

    """Plotting normalized p difference from 1"""
    plot_p_diff(np.array(t), np.array(DS2_long), "double_slit_p_diff")
    plot_p_diff(np.array(t), np.array(DS), "no_slit_p_diff")
#
    """Animation plot"""
    animation_plot(np.array(DS2_long_real), h, dt)
    animation_plot(np.array(DS2_long), h, dt)

    """Plotting probability at chosen x and t"""
    plot_at_x(np.array(DS2), 0.002, dt, 0.8,
              save="prob_y_normalized", title="Double-slit")
    plot_at_x(np.array(DS1), 0.002, dt, 0.8,
              save="single_prob_y_normalized", title="Single-slit")
    plot_at_x(np.array(DS3), 0.002, dt, 0.8,
              save="tripple_prob_y_normalized", title="Tripple-slit")

    """Plotting slit experiments for chosen list of t"""
    plot_time([0, 0.001, 0.002], dt, np.array(DS1),
              np.array(V_1), "single_slit", "$p_{i_j}$")
    plot_time([0, 0.001, 0.002, 0.005, 0.006], dt, np.array(DS2_long),
              np.array(V_2), "no_slit", "$p_{i_j}$")
    plot_time([0, 0.001, 0.002], dt, np.array(DS2),
              np.array(V_2), "double_slit", "$p_{i_j}$")
    plot_time([0, 0.001, 0.002], dt, np.array(DS3),
              np.array(V_3), "tripple_slit", "$p_{i_j}$")
    plot_time([0, 0.001, 0.002], dt, np.array(DS2_real),
              np.array(V_2), "double_slit_real", "Re$(u_{ij})$")
    plot_time([0, 0.001, 0.002], dt, np.array(DS2_imag),
              np.array(V_2), "double_slit_imag", "Im$(u_{ij})$")


if __name__ == "__main__":
    main()
