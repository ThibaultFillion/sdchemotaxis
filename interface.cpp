#include "simulation.hpp"
#include <chrono>
#include <iostream>
#include <omp.h>

Simulation * global_simulations;

extern "C" void E__run(
    //simulation input
    int w, 
    int h,
    const char * landscape,
    int sample_on_cell_jump,
    int sample_final_state,
    double sampling_interval,
    int sample_chemoattractant_state,
    int end_on_t_max,
    double t_max,
    double D,
    double ks,
    double v,
    double L0,
    double dt,
    int * rng_seeds,
    int n_runs,
    //output (number of samples for each run)
    int * n_samples
    )
    {
    global_simulations = new Simulation[n_runs];

    for(int n=0; n<n_runs; n++)
        {
        global_simulations[n].initialize(
            w,
            h,
            std::string(landscape),
            sample_on_cell_jump,
            sample_final_state,
            sampling_interval,
            sample_chemoattractant_state,
            end_on_t_max,
            t_max,
            D,
            ks,
            v,
            L0,
            dt,
            rng_seeds[n]
            );
        }
    
    #pragma omp parallel for
    for(int n=0; n<n_runs; n++)
        {
        std::cout<<"*";
        while(!global_simulations[n].get_complete())
            {
            global_simulations[n].iterate();
            }
        }

    for(int n=0; n<n_runs; n++)
        {
        n_samples[n] = global_simulations[n].get_n_samples();
        }
    }

extern "C" void E__get_cells_trajectories(int * out_cells_trajectories, int run_index)
    {
    const std::vector<std::vector<int>> & cells_trajectories = global_simulations[run_index].get_cells_trajectories();
    int n_samples = global_simulations[run_index].get_n_samples();
    int n_cells = global_simulations[run_index].get_n_cells();
 
    for(int i=0;i<n_cells;i++)
        {
        for(int j=0;j<n_samples;j++)
            {
            out_cells_trajectories[i * n_samples + j] = cells_trajectories[i][j];
            }
        }    
    }

extern "C" void E__get_trajectory_times(double * out_trajectory_times, int run_index)
    {
    const std::vector<double> & trajectory_times = global_simulations[run_index].get_trajectory_times();
    int n_samples = global_simulations[run_index].get_n_samples();
 
    for(int i=0;i<n_samples;i++)
        {
        out_trajectory_times[i] = trajectory_times[i];
        }
    }

extern "C" void E__get_chemoattractant_trajectory(double * out_chemoattractant_trajectory, int run_index)
    {
    const std::vector<std::vector<double>> & chemoattractant_trajectory = global_simulations[run_index].get_chemoattractant_trajectory();
    int n_samples = global_simulations[run_index].get_n_samples();
    int n_nodes = global_simulations[run_index].get_n_nodes();

    for(int i=0;i<n_samples;i++)
        {
        for(int j=0;j<n_nodes;j++)
            {
            out_chemoattractant_trajectory[i * n_nodes + j] = chemoattractant_trajectory[i][j];
            }
        }
    }

/*
Those source files contains the implementation of the simulation algorithm,
interfaced through the sdchemotaxis module.
Self-driven chemotactic behavior has been 
extensively studied, both experimentally and numerically, by Tweedy et al., who 
showed how this mechanism allows cells to find their way in complex environments such as mazes [Tweedy2020]_.
The implementation relies on the tau-leap algorithm [Gillespie2001]_ for reaction-diffusion, using 
Bernstein's approach for diffusion [Bernstein2005]_.
Please read the README.rst file for more information.

References:
===========
.. [Tweedy2020] Tweedy, L., Thomason, P. A., Paschke, P. I., Martin, K., Machesky, L. M., Zagnoni, M., & Insall, R. H. (2020). Seeing around corners: cells solve mazes and respond at a distance using attractant breakdown. Science, 369(6507), eaay9792. doi:10.1126/science.aay9792
.. [Gillespie2001] Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. The Journal of Chemical Physics, 115(4), 1716-1733. doi:10.1063/1.1378322
.. [Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical review. E, Statistical, nonlinear, and soft matter physics, 71(4), 041103. doi:10.1103/PhysRevE.71.041103
*/