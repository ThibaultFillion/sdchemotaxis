from sdchemotaxis import *

#simple example where one cell
#migrates through a fork with
#two chemostats (live-ends).
w = 4
h = 3
landscape = "#..*"\
            "o.##"\
            "#..*"

output = simulate(w, h, landscape)

# save the trajectory
# save_output(output, "_output")
# load the trajectory
# load_output("_output")

#plotting the state of the system at each sample
for i in range(output.n_samples()):
    output.plot_chem_state(i);
    plt.show()

#plotting the landscape of the system
output.plot_landscape()
plt.show()

#plotting the path of the cell
output.plot_cell_path()
plt.show()

# running n_runs realizations (runs) of the same stochastic simulation
# (the same as above), without sampling the chemoattractant state
n_runs = 10
output = simulate(w, h, landscape, sample_chem_state=False, n_runs=n_runs)

# save the trajectory
# save_output(output, "_output2")
# load the trajectory
# load_output("_output2")

# plotting the frequency at which the cell
# choses the left chemostat as a function
# of the number of runs
n = 0  # counts the number of time the cell chose the left chemostat
y = [] # frequency
x = [] # number of iterations
for i in range(n_runs):
    pos = output.get_cell_position(output.n_samples(i)-1, 0, i)
    if pos[1] == 0 : 
        n += 1
    
    y.append(n/(i+1))
    x.append(i+1)

plt.xlabel("fraction of left choice")
plt.ylabel("number of simulations")
plt.plot(x, y)
plt.show()

# plotting which nodes have been visited by the cell in the first run   
output.plot_cell_visit_map(run_index=0)
plt.show()

# plotting the frequency at which each node has been visited by the cell
# during the n_runs runs.
output.plot_cell_visit_frequency_map()
plt.show()