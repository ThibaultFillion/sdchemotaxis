import ctypes
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import pathlib
import json

"""
This module allows to simulate the chemotactic migration of cells
degrading the chemoattractant with a very simple model. Self-driven chemotactic behavior has been 
extensively studied, both experimentally and numerically, by Tweedy et al., who 
showed how this mechanism allows cells to find their way in complex environments such as mazes [Tweedy2020]_.
Pleas read the README.rst file for more information. The simulate function allow to run n parallell simulations
with the same initial conditions, and the Output class store the trajectories and other relevant data, and allow
for data manipulation and visualization, as well as two functions, load_output and save_output, to load and save
Output objects.

References
----------

.. [Tweedy2020] Tweedy, L., Thomason, P. A., Paschke, P. I., Martin, K., Machesky, L. M., Zagnoni, M., & Insall, R. H. (2020). Seeing around corners: cells solve mazes and respond at a distance using attractant breakdown. Science, 369(6507), eaay9792. doi:10.1126/science.aay9792
"""

def _process_landscape_array(landscape_array):
    """
    Returns w, h, landscape from an array landscape.
    """

    h = len(landscape_array)
    
    if h==0:
        raise ValueError("landscape array cannot be empty.")
    
    for y in range(h):
        if len(landscape_array[y]) != len(landscape_array[0]):
            raise ValueError("all rows of the landscape must have the same size.")
    
    w = len(landscape_array[0])
    if w==0:
        raise ValueError("landscape rows must be larger than 0.")
        
    for y in range(h):
        if not isinstance(landscape_array[y], str):
            raise ValueError("landscape array must be an array of str only.")

    return w, h, str().join(landscape_array)

def _process_landscape_arguments(w, h, landscape):
    """
    Processes the three landscape arguments.
    """

    if isinstance(landscape, list) or isinstance(landscape, tuple) or isinstance(landscape, np.ndarray):
        w, h, landscape = _process_landscape_array(landscape)
    
    if w is None :
        raise ValueError("no value specified for w.")

    if h is None :
        raise ValueError("no value specified for h.")

    if landscape is None :
        raise ValueError("no value specified for landscape.")

    if w*h != len(landscape):
        raise ValueError("The landscape length does not match the given dimensions.")

    return w, h, landscape        

def count_cells(landscape):
    """
    Returns the number of cells in a given landscape.
    """
    return landscape.count("o")
        
def plot_landscape(
        w = None, 
        h = None, 
        landscape = None, 
        ax=None, 
        colors={
            "wall": "grey", 
            "free": "yellow", 
            "chstt": "orange", 
            "cell": "green"
            }
        ):
    """
    Plots the landscape (of dimensions w*h).
    ax(optionnal) is the Matplotlib axis to be used.
    See simulate for the descrition of the three arguments.
    """
    w, h, landscape = _process_landscape_arguments(w, h, landscape)
            
    if ax is None :
        ax = plt.gca()
    ax.set_title("landscape")
    ax.set_xlabel("x (nodes)")
    ax.set_ylabel("y (nodes)")
    numlandscape = []
    nodecode = {
        "#" : 0.5,
        "." : 1.5,
        "*" : 2.5,
        "o" : 3.5
        }
    for i in range(w*h):
        numlandscape.append(nodecode[landscape[i]])

    color_map = matplotlib.colors.ListedColormap([
        colors["wall"],
        colors["free"],
        colors["chstt"],
        colors["cell"]
        ])
    im = ax.imshow(np.array(numlandscape).reshape(h, w), cmap=color_map, vmin=0, vmax=4)
    color_bar = plt.colorbar(im)    
    color_bar.set_ticks([0.5, 1.5, 2.5, 3.5], labels=["wall", "free node", "chemostat", "cell"])

class Output:
    """
    Describes a set of independent trajectories
    with the same initial conditions. It is generated
    as the output of run_simulations.
    
    attributes:
    
    * w, h, landscape: width, height and landscape of the systeme
    * n_runs, n_nodes, n_cells: number of runs (independant realizations of the simulation),
        nodes and migrating cells.
    * rng_seed: seeds used for the different runs.
    """
    
    def __init__(
            self,
            # landscape infos
            w, 
            h,
            landscape,

            # flags
            has_state_traj,

            # dimensions
            n_nodes, 
            n_samples, 
            n_cells, 
            n_runs,

            # data
            cell_trajs, 
            state_trajs, 
            times,
            
            #rng seeds
            rng_seeds
            ):
        self.w         = w
        self.h         = h
        self.landscape = landscape

        self._has_state_traj = has_state_traj

        self.n_nodes    = n_nodes
        self._n_samples = n_samples
        self.n_cells    = n_cells
        self.n_runs     = n_runs
        
        self._cell_trajs  = cell_trajs 
        self._state_trajs = state_trajs
        self._times        = times
        
        self.rng_seeds = rng_seeds

    def n_samples(self, run_index=0):
        """
        Returns the number of trajectory samples for the given run.
        """
        return self._n_samples[run_index]

    def times(self, run_index=0):
        """
        Returns the trajectory sampling times for the given run.
        """
        return self._times[run_index]

    def get_L_bounds(self, run_index=0):
        """
        Return the lowest and highest ([min, max]) values of the chemoattractant state for a given run. 
        """
        return min(self._state_trajs[run_index]), max(self._state_trajs[run_index])
     
    def get_chem_state(self, sample_index, run_index=0, walls_as_nans=False):
        """
        Return the chemoattractant state (quantity in each node) for a given sample
        in a given run. If walls_as_nans is True, the state is set to NaN 
        """
        if not self._has_state_traj :
            raise NotImplementedError("chemoattractant state was not sampled")
 
        state = self._state_trajs[run_index].reshape((self.n_samples(run_index), self.n_nodes))[sample_index].copy()
        if walls_as_nans:
            for i in range(self.w*self.h):
                if self.landscape[i] == "#":
                    state[i] = math.nan
        return state

    def get_cell_trajectory_indices(self, cell_index=0, run_index=0):
        """
        Return the trajectory of a given cell for a given run.
        The trajectory fratures cell positions as node indices.
        """
        return self._cell_trajs[run_index].reshape((self.n_cells, self.n_samples(run_index)))[cell_index]

    def get_cell_trajectory(self, cell_index=0, run_index=0):
        """
        Return the trajectory of a given cell for a given run.
        The trajectory fratures cell positions as node coordinates, and is
        shaped as [[x0,y0], ...].
        """
        traj = self.get_cell_trajectory_indices(cell_index, run_index)
        return np.array([[int(traj[i]%self.w), int(traj[i]/self.w)] for i in range(self.n_samples(run_index))])

    def get_cell_position_index(self, sample_index, cell_index=0, run_index=0):
        """
        Return the position (node index) of a given cell for a given run
        for a given sample.
        """
        return self._cell_trajs[run_index].reshape((self.n_cells, self.n_samples(run_index)))[cell_index][sample_index]

    def get_cell_position(self, sample_index, cell_index=0, run_index=0):
        """
        Return the position ([x,y] node coordinates) of a given cell for a given run
        for a given sample.
        """
        pos = self.get_cell_position_index(sample_index, cell_index, run_index)
        return int(pos%self.w), int(pos/self.w)

    def plot_chem_state(self, sample_index, run_index=0, ax=None):
        """
        Plots the chemostat state for a given sample in a given run.
        """
        if not self._has_state_traj :
            raise NotImplementedError("chemoattractant state was not sampled")
            
        if ax is None :
            ax = plt.gca()

        ax.set_title("$t$="+str(self.times(run_index)[sample_index])+"s")
        ax.set_xlabel("x (nodes)")
        ax.set_ylabel("y (nodes)")
        lmin, lmax = self.get_L_bounds(run_index)
        im = ax.imshow(self.get_chem_state(sample_index, run_index, walls_as_nans=True).reshape((self.h, self.w)), cmap="magma", vmin=lmin, vmax=lmax)
        for cell_index in range(self.n_cells):
            cell_pos = self.get_cell_position(sample_index, cell_index, run_index)
            ax.plot(cell_pos[0], cell_pos[1], color="green", linestyle="", marker="o")
        plt.colorbar(im, label="[L] (molecule/node)")

    def plot_landscape(self, ax=None):
        """
        Plots the system landscape.
        """
        plot_landscape(self.w, self.h, self.landscape, ax)

    def plot_cell_path(self, cell_index=0, run_index=0, ax=None):
        """
        Plots path floowed by a given cell for a given run..
        """
        if ax is None :
            ax = plt.gca()
        
        colors = {
            "wall":  (0.9,0.9,0.9), 
            "free":  (1,1,0.8),
            "chstt": (1,0.9,0.8),
            "cell":  (0.8,1,0.8)
            }
        plot_landscape(self.w, self.h, self.landscape, ax, colors)
        traj = self.get_cell_trajectory(cell_index, run_index).transpose()
        ax.plot(traj[0], traj[1], label="cell "+str(cell_index))
        if self.n_cells>1:
            ax.legend()
    
    def build_cell_visit_map(self, cell_index=0, run_index=0):
        """
        Returns an array with values indicating, for each
        node, whether the given cell, in a given run, visited the node at least once.
        values for a node is 1 if visited, 0 otherwise.
        """
        vmap = np.zeros(self.n_nodes)
        traj = self.get_cell_trajectory_indices(cell_index, run_index)
        for i in range(self.n_samples(run_index)):
            vmap[traj[i]] = 1.0
        return vmap

    def plot_cell_visit_map(self, cell_index=0, run_index=0, ax=None):
        """
        Plots the cell visit map.
        """
        if ax is None :
            ax = plt.gca()
        
        ax.set_title("cell visit map")
        ax.set_xlabel("x (nodes)")
        ax.set_ylabel("y (nodes)")
        vmap = self.build_cell_visit_map(cell_index, run_index)        
        color_map = matplotlib.colors.ListedColormap([(68/255,1/255,84/255), (253/255,231/255,36/255)])
        im = ax.imshow(vmap.reshape((self.h, self.w)), cmap=color_map)
        color_bar = plt.colorbar(im)
        color_bar.set_ticks([0.25, 0.75], labels=["not visited", "visited"])

    def build_cell_visit_frequency_map(self, cell_index=0):
        """
        Returns an array with values indicating, for each
        node, the frequency at which the given cell visits the node across all runs.
        """
        pmap = np.zeros(self.n_nodes)
        for n in range(self.n_runs):
            vmap = self.build_cell_visit_map(cell_index, n)
            pmap += vmap
        return pmap/self.n_runs

    def plot_cell_visit_frequency_map(self, cell_index=0, ax=None):
        """
        Plots the cell visit frequency map.
        """
        if ax is None :
            ax = plt.gca()
        
        ax.set_title("cell visit frequency map")
        ax.set_xlabel("x (nodes)")
        ax.set_ylabel("y (nodes)")
        pmap = self.build_cell_visit_frequency_map(cell_index)        
        im = ax.imshow(pmap.reshape((self.h, self.w)))
        plt.colorbar(im, label="visit frequency")
        
def simulate(
    w=None, 
    h=None,
    landscape=None,
    sample_on_cell_jump = True,
    sample_final_state = True,
    sampling_interval = -1,
    sample_chem_state = True,
    end_on_t_max = False,
    t_max = 0,
    D = 50,
    ks = 1,
    v = 0.5/60, # cell velocity in µm/s
    L0 = 593.0/20.0**3, # chemoattractant density in molecules/µm3
    time_step = -1, # negative means auto time step
    rng_seeds = None,
    n_runs = 1
    ):
    """
    Run simulations of cell migration in a given landscape with some set of parameters.
    There are two ways to call this functions:
    with, the landscape as a str, and explicit width and height, i.e.

        simulate(2,2,"o..*") 
    
    or with the landscape as an array, in keyword argument, and no width and height (determined automatically), i.e.
    
        simulate(landscape=["o.",
                            ".*"]) 
    
    * w: width of the system (in nodes, default=None)
    * h: height of the system (in nodes, default=None)
    * landscape: description of the landscape at initial state (str, default=None).
    * sample_on_cell_jump: should the trajectories be sampled on cell jump? (bool, default=True)
    * sample_final_state: should the trajectories be sampled when the end condition is met? (bool, default=True)
    * sampling_interval: interval at which trajectories should be sampled. 
        If negative, no sampling on interval. (default=-1)
    * sample_chem_state: should the chemoattractant profile be sampled too? (bool, default=True)
    * end_on_t_max: should the simulation end whent t_max is reacched? (defaut=False)
    * t_max: time limit of the simulation (default=0).
    * D: diffusion coefficient of the chemoattractant (in µm2/s, default=50)
    * ks: cell sink rate (in s-1, default=1)
    * v: velocity of the cell (in µm/s, default=0.5/60).
    * L0: chemoattractant density at the initial state and at the chemostats (default=593.0/20.0**3).
    * time_step: time step. A negative value means is is determined automatically (default=-1)
    * rng_seeds: seeds for the pseudo-random number generator (list of n_runs values). if None, they are generated. (default=None)
    * n_runs: Number of runs to be performed (default=1).
    """
    w, h, landscape = _process_landscape_arguments(w, h, landscape)
    
    if rng_seeds is None:
        rng_seeds = [int(random.random()*1000000) for i in range(n_runs)]
    
    _rng_seeds = (ctypes.c_int * n_runs)()
    for i in range(n_runs):
        _rng_seeds[i] = rng_seeds[i]
    
    sl = ctypes.CDLL(str(pathlib.Path(__file__).parent/"sl.dll"))
    
    n_samples = (ctypes.c_int * n_runs)()
    
    sl.E__run(
        ctypes.c_int(w), 
        ctypes.c_int(h),
        ctypes.c_char_p(landscape.encode()),
        ctypes.c_int(sample_on_cell_jump),
        ctypes.c_int(sample_final_state),
        ctypes.c_double(sampling_interval),
        ctypes.c_int(sample_chem_state),
        ctypes.c_int(end_on_t_max),
        ctypes.c_double(t_max),
        ctypes.c_double(D),
        ctypes.c_double(ks),
        ctypes.c_double(v),
        ctypes.c_double(L0),
        ctypes.c_double(time_step),
        _rng_seeds,
        ctypes.c_int(n_runs),
        n_samples
        )

    n_samples = [int(n_samples[i]) for i in range(n_runs)]
    n_cells = count_cells(landscape)
    n_nodes = w*h

    trajs_state = [None for i in range(n_runs)]
    trajs_cells = [None for i in range(n_runs)]
    trajs_times = [None for i in range(n_runs)]
    
    for n in range(n_runs):
        trajs_state[n] = np.array([])
        trajs_cells[n] = (ctypes.c_int * int(n_cells * int(n_samples[n])))()
        trajs_times[n] = (ctypes.c_double * int(n_samples[n]))()
        sl.E__get_cells_trajectories(trajs_cells[n], ctypes.c_int(n))
        sl.E__get_trajectory_times(trajs_times[n], ctypes.c_int(n))
        
        if sample_chem_state:
            trajs_state[n] = (ctypes.c_double * int(n_nodes * int(n_samples[n])))()
            sl.E__get_chemoattractant_trajectory(trajs_state[n], ctypes.c_int(n))
            
        trajs_state[n] = np.array(trajs_state[n])
        trajs_cells[n] = np.array(trajs_cells[n])
        trajs_times[n] = np.array(trajs_times[n])
            
    return Output(
        w,
        h,
        landscape,
    
        sample_chem_state,
        
        n_nodes, 
        n_samples, 
        n_cells, 
        n_runs,
        
        trajs_cells, 
        trajs_state, 
        trajs_times,
        
        rng_seeds
        )

def save_output(output, path):
    """
    save an Output object.
    path is not the file name, but rather
    the base of the filename:
    upon saving, two files are created:
    
    * path_info.json, and
    * path_data.npz.
    """
    file = open(path + "_info.json", "w")
    info = {
        "w": output.w,
        "h": output.h,
        "landscape": output.landscape,
        "_has_state_traj": output._has_state_traj,
        "n_nodes": output.n_nodes,
        "_n_samples": list(output._n_samples),
        "n_cells": output.n_cells,
        "n_runs": output.n_runs,
        "rng_seeds": list(output.rng_seeds)
        }
    json.dump(info, file)
    file.close()
    
    data = {}
    for i in range(output.n_runs):
        data["c"+str(i)] = output._cell_trajs[i]
        data["s"+str(i)] = output._state_trajs[i]
        data["t"+str(i)] = output._times[i]        
    np.savez(path + "_data.npz", **data)

def load_output(path):
    """
    load an Output object.
    path is not the file name, but rather
    the base of the filename:
    loading path amounts to load path_info.json and
    path_data.npz.
    """
    file = open(path+"_info.json", "r")
    info = json.load(file)
    file.close()
    
    data = np.load(path + "_data.npz")
    _cell_trajs  = [None for i in range(info["n_runs"])]
    _state_trajs = [None for i in range(info["n_runs"])]
    _times       = [None for i in range(info["n_runs"])]
    for i in range(info["n_runs"]):
        _cell_trajs[i]  = data["c"+str(i)]
        _state_trajs[i] = data["s"+str(i)]
        _times[i]       = data["t"+str(i)]
        
    Output(
        info["w"],
        info["h"],
        info["landscape"],
        info["_has_state_traj"],
        info["n_nodes"],
        info["_n_samples"],
        info["n_cells"],
        info["n_runs"],
        _cell_trajs,
        _state_trajs,
        _times,
        info["rng_seeds"]
        )
