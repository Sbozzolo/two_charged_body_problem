#!/usr/bin/env python3
"""Solve the Newtonian two body problem with charges and radiation reaction.
Extract gravitational waves using the quadrupole approximation.

"""

# This package is a single-file module, so we adopt the unconventional choice
# of putting everything in the __init__

import numpy as np
from scipy.integrate import odeint


class two_charged_body_problem:
    """Solve the two body problem in the Newtonian approximation, including
    radiation reaction.

    We work in units in which G=c=M_tot=1.

    When output_in_cgs is set to True, a mass scale has to be provided with
    total_mass.

    """

    # Some useful constants to convert
    # (PEP 515, Underscores in Numeric Literals)
    C_CGS = 29_979_245_800  # Vacuum speed of light

    G_CGS = 6.673e-8  # Gravitational constant
    M_SOL_CGS = 1.98892e33  # Solar mass

    def __init__(self,
                 mass_ratio=1,
                 lambda1=0,
                 lambda2=0,
                 output_in_cgs=False,
                 total_mass=1,
                 enable_em_dipole_radiation_reaction=True,
                 enable_em_quadrupole_radiation_reaction=True,
                 enable_gw_quadrupole_radiation_reaction=True):

        self.output_in_cgs = output_in_cgs

        # Conversion to cgs
        self.geom_mass_to_M_sun = total_mass
        self.geom_energy_to_erg = total_mass * self.M_SOL_CGS * self.C_CGS**2
        self.geom_time_to_s = self.G_CGS * self.M_SOL_CGS / self.C_CGS**3 * total_mass
        self.geom_freq_to_Hz = 1 / self.geom_time_to_s
        self.geom_lum_to_erg_per_s = self.geom_energy_to_erg / self.geom_time_to_s

        # Now we can set M_tot to 1
        self.M_tot = 1

        # Mass ratio
        self.q = mass_ratio
        # Individial masses
        self.m1 = self.q * self.M_tot / (1 + self.q)
        self.m2 = self.M_tot - self.m1
        # Chirp mass
        self.M_c = ((self.m1**3 * self.m2**3) / self.M_tot)**(1 / 5)
        # Mass-to-charge ratios
        self.lambda1, self.lambda2 = lambda1, lambda2

        # Enable only some of the energy-loss mechanisms
        self.pemd = 1 if enable_em_dipole_radiation_reaction else 0
        self.pemq = 1 if enable_em_quadrupole_radiation_reaction else 0
        self.pgwq = 1 if enable_gw_quadrupole_radiation_reaction else 0

    def P_GW_quadrupole(self, R):
        return 32 * (1 - self.lambda1 * self.lambda2
                     )**3 * self.m1**2 * self.m2**2 * self.M_tot / (5 * R**5)

    def P_EM_dipole(self, R):
        return 2 / 3 * (self.lambda1 - self.lambda2)**2 * (
            1 -
            self.lambda1 * self.lambda2)**2 * self.m1**2 * self.m2**2 / R**4

    def P_EM_quadrupole(self, R):
        return self.m1**2 * self.m2**2 / self.M_tot**2 * (
            self.lambda1 / self.m1 +
            self.lambda2 / self.m2)**2 / 4 * self.P_GW_quadrupole(R)

    def P(self, R):
        return self.pemd * self.P_EM_dipole(
            R) + self.pemq * self.P_EM_quadrupole(
                R) + self.pgwq * self.P_GW_quadrupole(R)

    def GW_frequency(self, r):
        return 1 / np.pi * np.sqrt(
            (1 - self.lambda1 * self.lambda2) * self.M_tot / r**3)

    def GW_angular_velocity(self, r):
        return 2 * np.pi * self.GW_frequency(r)

    # The total energy of the system is
    # E = - (1 - lambda1 * lambda2) * m1 * m2 / (2 * R)

    # The equation of motion is \dot{E} = P_GW + P_EM
    # \dot{E} = (1 - lambda1 * lambda2) * m1 * m2 / (2 * R**2) * \dot{R}

    # So, \dot{R} = 2 * R**2 / ((1 - lambda1 * lambda2) * m1 * m2) * (P_GW + P_EM)

    # t is a dummy variable
    def evolution_equation(self, R, t):
        return -2 * R**2 / (self.m1 * self.m2 *
                            (1 - self.lambda1 * self.lambda2)) * self.P(R)

    def compute_gws(self):
        # To compute the gravitational wave phase we have to integrate the
        # gravitational wave frequency. The most computationally efficient way
        # to do this is to use a simple quadrature. Now, we create a fake point
        # before the evolution so that we can consider the points we already
        # have "midpoints"

        R0_fake = odeint(self.evolution_equation, self.R[0],
                         [self.t[0], self.t[0] + self.R[0] * self.dt_R])[-1][0]

        R_fake = np.append(R0_fake, self.R)
        # t is needed for total energies
        t_fake = np.append(self.t[0] + self.R[0] * self.dt_R, self.t)

        # The minus is needed to have everything positive.
        # Otherwise, we would have negative delta_R
        delta_R = -np.diff(R_fake)
        delta_t = np.diff(t_fake)

        self.GW_frequencies = self.GW_frequency(self.R)
        self.GW_angular_velocities = self.GW_angular_velocity(self.R)
        self.GW_phases = np.cumsum(self.GW_angular_velocities * delta_R /
                                   -self.evolution_equation(self.R, 0))

        # Maggiore (4.29)
        h0_rex = 4 * ((1 - self.lambda1 * self.lambda2) * self.M_c)**(
            5 / 3) * (np.pi * self.GW_frequencies)**(2 / 3)

        h_plus_rex = h0_rex * np.cos(self.GW_phases)
        h_cross_rex = h0_rex * np.sin(self.GW_phases)

        self.hp_rex = h_plus_rex
        self.hc_rex = h_cross_rex

        # Equation (4.22) in Maggiore
        self.GW_cycles = np.cumsum(self.GW_frequencies * delta_t)

        self.inst_P_EM_dipole = self.P_EM_dipole(self.R)
        self.inst_P_EM_quadrupole = self.P_EM_quadrupole(self.R)
        self.inst_P_GW_quadrupole = self.P_GW_quadrupole(self.R)
        self.E_EM_dipole = np.cumsum(self.inst_P_EM_dipole * delta_t)
        self.E_EM_quadrupole = np.cumsum(self.inst_P_EM_quadrupole * delta_t)
        self.E_GW_quadrupole = np.cumsum(self.inst_P_GW_quadrupole * delta_t)

        if (self.output_in_cgs):
            self.GW_frequencies *= self.geom_freq_to_Hz
            self.GW_angular_velocities *= self.geom_freq_to_Hz
            self.t *= self.geom_time_to_s
            self.inst_P_EM_dipole *= self.geom_lum_to_erg_per_s
            self.inst_P_EM_quadrupole *= self.geom_lum_to_erg_per_s
            self.inst_P_GW_quadrupole *= self.geom_lum_to_erg_per_s
            self.E_EM_dipole *= self.geom_energy_to_erg
            self.E_EM_quadrupole *= self.geom_energy_to_erg
            self.E_GW_quadrupole *= self.geom_energy_to_erg

        self.inst_P_EM = self.inst_P_EM_dipole + self.inst_P_EM_quadrupole
        self.inst_P_GW = self.inst_P_GW_quadrupole
        self.E_EM = self.E_EM_dipole + self.E_EM_quadrupole
        self.E_GW = self.E_GW_quadrupole

    def solve(
            self,
            initial_separation=10,  # Units of total_mass
            integration_step=0.005,  # Units of separation
            final_separation=4):  # Units of total_mass
        # Initial separation
        self.R_initial = initial_separation * self.M_tot
        # Integration step (in units of R)
        # The integration step varies with r, so that it is smaller
        # when the separation is smaller
        self.dt_R = integration_step
        # Final R is set to the ISCO
        self.R_final = final_separation * self.M_tot

        # Prepare array of results
        self.R = np.array([self.R_initial])
        self.t = np.array([0])

        while True:
            dt = self.R[-1] * self.dt_R
            R_new = odeint(self.evolution_equation, self.R[-1],
                           [self.t[-1], self.t[-1] + dt])[-1][0]
            self.R, self.t = np.append(self.R, R_new), np.append(
                self.t, self.t[-1] + dt)
            if R_new < self.R_final:
                break

        self.compute_gws()

    def time_at_frequency(self, freq):
        """
        Use cgs if output_cgs is True
        """
        return self.t[np.argmin(np.abs(self.GW_frequencies - freq))]

    def index_at_time(self, t):
        return np.argmin(np.abs(self.t - t))

    def quantity_in_t_range(self, quantity, t_min, t_max):
        """ Use cgs if output_cgs is True """
        i_min, i_max = self.index_at_time(t_min), self.index_at_time(t_max)
        return quantity[i_max] - quantity[i_min]

    def E_EM_dipole_t_range(self, t_min, t_max):
        return self.quantity_in_t_range(self.E_EM_dipole, t_min, t_max)

    def E_EM_quadrupole_t_range(self, t_min, t_max):
        return self.quantity_in_t_range(self.E_EM_quadrupole, t_min, t_max)

    def E_GW_quadrupole_t_range(self, t_min, t_max):
        return self.quantity_in_t_range(self.E_GW_quadrupole, t_min, t_max)

    def GW_cycles_t_range(self, t_min, t_max):
        return self.quantity_in_t_range(self.GW_cycles, t_min, t_max)
