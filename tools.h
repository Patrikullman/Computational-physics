#pragma once

void velocity_verlet(int n_timesteps, int n_particles, double v[][3], double position[][3],
double position_evolution[][n_particles][3], double dt, double m, double force[][3], 
double cell_length);

void write_to_file(char *fname, double *times,
            double *pot_energy, double *kin_energy, 
            double *tot_energy, int n_points);