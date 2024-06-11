#include<random>
#include<vector>
#include<valarray>
#include<cmath>

struct COORDINATES {int x; int y;};
struct EDGE {int i; int j;};

struct CELL
    {
    int position;
    std::valarray<double> sensing;
    double last_jump_time;
    };

CELL create_cell(int position, double last_jump_time)
    {
    return {
        position,
        {0,0,0,0},
        last_jump_time
        };
    }

class Simulation
    {
    private :

    std::mt19937 rng; // pseudo-random number generator
    int w; // width of the migration environments
    int h; // height of the migration environments
    int n_cells; // number of cells
    int n_edges; // edges.size() as an int
    std::vector<bool> chstts;  // chemostats
    std::vector<double> chstt_factors;  // 1-chemostats as a double
    std::vector<bool> occupied;  // indicates whether a node is occupied by a cell or not
    std::vector<double> state; // chemoattractant distribution
    std::vector<double> nonwalls;  // 0 wall, 1 free node
    std::vector<CELL> cells; // migrating cells
    std::vector<std::vector<double>> state_trajectory;
    std::vector<std::vector<int>> cells_trajectories;
    std::vector<double> trajectory_times;
    std::vector<EDGE> edges;
    std::vector<double> diffusion_ij; //number of diffusion events for each edge, from i to j
    std::vector<double> diffusion_ji; //number of diffusion events for each edge, from j to i
    std::vector<double> sink; //number of sink reaction events in each node with a cell
    double t; // time
    double dt; // time step
    double ks; // sink rate (s-1)
    double kd; // diffusion rate (s-1)
    double cell_jump_interval; // interval at which cells jump (s-1)
    bool sample_on_cell_jump; // tells whether the state should be sampled on cell jumps
    bool sample_final_state; // tells whether the final state should be sampled
    double sampling_interval; // trajecry sampling interval. if <0, no sampling on interval
    bool sample_chemoattractant_state; //if true, the distribution of L is also sampled.
    bool cell_jumped; //tells if a cell jumps during the current iteration.
    double t_max; // time limit of the simulation
    bool end_on_t_max;
    bool complete; //tells if the simulation is complete

    std::uniform_real_distribution<double> unit_rand_distribution; //uniform random dist between 0 and 1

    double poisson(double lambda)
        {
        return static_cast<double>(std::poisson_distribution<int>(lambda)(rng));
        }

    double unit_rand()
        {
        return unit_rand_distribution(rng);
        }
    
    std::valarray<double> nullify_negative_values(std::valarray<double> v)
        {
        for(size_t i=0; i<v.size(); i++)
            {
            v[i] = std::max(v[i], 0.0);
            }
        return v;
        }

    COORDINATES get_coordinates(int node_index)
        {
        return {node_index%w, node_index/w};
        }

    bool is_within_bounds(int x, int y)
        {
        return (x>=0 && y>=0 && x<w && y<h);
        }

    double get_field_value(const std::vector<double> & field, int x, int y, double default_value=0)
        {
        // returns the value of the field at the given  
        // coordinates, or the default_value if out of bounds
        if(is_within_bounds(x, y))
            {
            return field[x+y*w];
            }
        else
            {
            return default_value;
            }            
        }

    std::valarray<double> get_neighborhood(const std::vector<double> & field, int i)
        {
        // returns the value of the field at the 
        // given coordinates, or 0 if out of bounds
        COORDINATES pos = get_coordinates(i);
        return {
            get_field_value(field, pos.x, pos.y+1), // up
            get_field_value(field, pos.x, pos.y-1), // down
            get_field_value(field, pos.x-1, pos.y), // left
            get_field_value(field, pos.x+1, pos.y), // right
            };             
        }

    void parse_landscape(const std::string & landscape, double L0)
        {
        chstts = std::vector<bool>(landscape.size(), false);
        occupied = std::vector<bool>(landscape.size(), false);
        chstt_factors = std::vector<double>(landscape.size(), 1.0);
        state = std::vector<double>(landscape.size(), 0);
        nonwalls = std::vector<double>(landscape.size(), 1);
        for(size_t i=0; i<landscape.size(); i++)
            {
            switch(landscape[i])
                {
                case '#' : 
                    nonwalls[i] = 0;
                    break;
                case '.' : 
                    state[i]  = L0;
                    break;
                case '*' : 
                    state[i]  = L0;
                    chstts[i] = true;
                    chstt_factors[i] = 0.0;
                    break;
                case 'o' : 
                    state[i]  = L0;
                    occupied[i] = true;
                    cells.push_back(create_cell(i, unit_rand()*cell_jump_interval));
                    break;
                }
            }
        n_cells = static_cast<int>(cells.size());
        }
    
    void build_edges()
        {
        //define edges (connexions between free nodes) without duplicates 
        for(int y=0; y<h; y++)
            {
            for(int x=0; x<w; x++)
                {
                int i = y*w+x;
                if(nonwalls[i])
                    {
                    int j = y*w+x+1;
                    if(x<w-1 && nonwalls[j])
                        {
                        edges.push_back({i, j});
                        }
                    j = (y+1)*w+x;
                    if(y<h-1 && nonwalls[j])
                        {
                        edges.push_back({i, j});
                        }
                    }
                }
            }
        n_edges = static_cast<int>(edges.size());
        }

    void init_event_arrays()
        {
        diffusion_ij.resize(n_edges);
        diffusion_ji.resize(n_edges);
        sink.resize(n_cells);
        }

    void tauleap_reaction_diffusion()
        {
        // reaction-diffusion handled with the tau-leap algorithm [Gillespie2001] 
        // with diffusion handled as descrived by [Bernstein2005]

        // compute the number of events
        for(int i=0; i<n_edges; i++)
            {
            diffusion_ij[i] = poisson(kd * state[edges[i].i] * dt);
            diffusion_ji[i] = poisson(kd * state[edges[i].j] * dt);
            }
        for(int i=0;i<n_cells; i++)
            {
            sink[i] = poisson(ks * state[cells[i].position] * dt);  
            }

        // apply events
        for(int i=0; i<n_edges; i++)
            {
            state[edges[i].i] += (diffusion_ji[i] - diffusion_ij[i])*chstt_factors[edges[i].i];
            state[edges[i].j] += (diffusion_ij[i] - diffusion_ji[i])*chstt_factors[edges[i].j];
            }
        for(int i=0;i<n_cells; i++)
            {
            state[cells[i].position] -= sink[i]*chstt_factors[cells[i].position];
            }
        }

    void cell_migration()
        {
        for(int i=0; i<n_cells; i++)
            {
            // sample neighboring concentrations
            cells[i].sensing += get_neighborhood(state, cells[i].position) * dt/cell_jump_interval;

            // jump
            if(t - cells[i].last_jump_time >= cell_jump_interval)
                {
                jump_cell(i);
                }
            }
        }

    void jump_cell(int cell_index)
        {
        CELL & cell = cells[cell_index];
        cell_jumped = true;

        // compute weights. negative values (artifacts) are treated as 0
        std::valarray<double> weights = pow(nullify_negative_values(cell.sensing), 6);
        
        // handle the special case where no chemoattractant
        // is sensed by eliminating wall nodes
        if(weights.max() == 0)
            {
            weights = get_neighborhood(nonwalls, cell.position);
            }

        // drawing the new direction
        int direction = std::discrete_distribution<int>(std::begin(weights), std::end(weights))(rng);
        
        // do the jump if the destination node is free
        int next_cell_position = compute_jump_pos(cell.position, direction); 
        if(!occupied[next_cell_position])
            {
            occupied[cell.position] = false;
            occupied[next_cell_position] = true;
            cell.position = next_cell_position;
            }
        
        // reset cell sensing and update last jump time
        cell.sensing = {0,0,0,0};
        cell.last_jump_time = t;
        }

    int compute_jump_pos(int position, int direction)
        {
        switch(direction)
            {
            case 0: return position+w;
            case 1: return position-w;
            case 2: return position-1;
            case 3: return position+1;
            }
        return -1;
        }

    void check_end_conditions()
        {
        // check t_max if necessary
        if(t>=t_max && end_on_t_max)
            {
            complete = true;
            return;
            }
        
        // check if cells have reached a chemostat
        bool all_cells_on_chemostats = true;
        for(int i=0; i<n_cells; i++)
            {
            all_cells_on_chemostats &= chstts[cells[i].position];
            }
        if(all_cells_on_chemostats)
            {
            complete = true;
            }
        }

    void sample_step()
        {
        // checks the sampling conditions and sample if necessary
        if(
            (cell_jumped && sample_on_cell_jump) ||
            (sampling_interval>=0 && t - trajectory_times[0]>=sampling_interval)||
            (sample_final_state && complete)
            )
            {
            sample();
            }
        }

    void sample()
        {
        // do the actual sampling
        trajectory_times.push_back(t);
        for(int i=0; i<n_cells; i++)
            {
            cells_trajectories[i].push_back(cells[i].position);
            }
        if(sample_chemoattractant_state)
            {
            state_trajectory.push_back(state);
            }
        }

    public :

    void initialize(
        // dimensions of the grid
        int w, 
        int h,

        // system landscape 
        std::string landscape,

        // sampling
        bool sample_on_cell_jump, // tells if sampling should be done after each cell jump
        bool sample_final_state,
        int sampling_interval, // if negative, no sampling on interval
        bool sample_chemoattractant_state, // if true the chemoattractant state is also sampled

        // end conditions
        bool end_on_t_max,
        double t_max,

        // kinetic parameters
        double D,  // chemoattractant diffusion coefficient (µm2/s)
        double ks, // cell sink rate (s-1)
        double v,  // cell speed (µm/s)
        double L0, // initial + chstt chemoattractant density (molecules/µm3)

        // time step
        double dt,

        // rng
        int rng_seed
        )
        {
        //kinetic parameters
        double node_edge = 20; // µm
        double node_volume = pow(node_edge, 3); // µm3
        this->cell_jump_interval = node_edge/v; // s/node_edge       
        this->ks = ks; // s-1
        this->kd = D/pow(h,2); // s-1, diffusion rate computed according to [Bernstein2005]

        // time and time step
        this->t = 0;
        if(dt<=0)
            {
            this->dt = 0.8/(4*kd+ks);
            }
        else
            {
            this->dt = dt;
            }        
                
        //system configuration and state
        this->w = w;
        this->h = h;
        parse_landscape(landscape, L0*node_volume); // L0 in molecules/node
        build_edges();
        init_event_arrays();

        //rng initialization
        this->rng = std::mt19937(rng_seed);
        this->unit_rand_distribution = std::uniform_real_distribution<double>(0,1);
        
        //sampling related info
        this->cells_trajectories.resize(n_cells);
        this->sample_on_cell_jump = sample_on_cell_jump;
        this->sample_final_state = sample_final_state; 
        this->sampling_interval = sampling_interval;
        this->sample_chemoattractant_state = sample_chemoattractant_state;
        this->end_on_t_max = end_on_t_max;
        this->t_max = t_max;
        this->complete = false;
        
        // sample initial state
        sample(); 
        }

    void iterate()
        {
        cell_jumped = false;
        tauleap_reaction_diffusion();
        cell_migration();
        check_end_conditions();
        t = t + dt; 
        sample_step();
        }
    
    const std::vector<double> & get_trajectory_times()
        {
        return trajectory_times;
        }

    const std::vector<std::vector<int>> & get_cells_trajectories()
        {
        return cells_trajectories;
        }

    const std::vector<std::vector<double>> & get_chemoattractant_trajectory()
        {
        return state_trajectory;
        }

    int get_n_cells()
        {
        return n_cells;
        }

    int get_n_samples()
        {
        return static_cast<int>(trajectory_times.size());
        }

    int get_n_nodes()
        {
        return w*h;
        }
    
    int get_complete()
        {
        return complete;
        }
    
    };

/*
References:
-----------
[Gillespie2001] Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. The Journal of Chemical Physics, 115(4), 1716-1733. doi:10.1063/1.1378322
[Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical review. E, Statistical, nonlinear, and soft matter physics, 71(4), 041103. doi:10.1103/PhysRevE.71.041103
*/