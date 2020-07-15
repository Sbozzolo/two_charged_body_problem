#!/usr/bin/env python3

if __name__ == '__main__':

    from tcbp import two_charged_body_problem

    M_tot = 65

    f_min = 30
    f_mid = 55
    f_max = 120

    # lambdas = [0.01, 0.05, 0.1, 0.2, 0.3]
    lambdas = [0.1]

    sims = [(0, 0)]
    sims += [(l, l) for l in lambdas]
    sims += [(l, -l) for l in lambdas]
    sims += [(l, 0) for l in lambdas]

    print("# lambda1 lambda2 t_30_120 E_GW_30_120 E_EM_30_120 GW_cycles_30_120 t_30_55 E_GW_30_55 E_EM_30_55 GW_cycles_30_55")
    print("# [] [] [M] [M] [M] [#] [M] [M] [M] [#]")

    for s in sims:
        tcbp = two_charged_body_problem(total_mass=M_tot,
                                        lambda1=s[0],
                                        lambda2=s[1],
                                        mass_ratio=29/36)
        tcbp.solve(initial_separation=12.1)
        t_min = tcbp.time_at_frequency(f_min / tcbp.geom_freq_to_Hz)
        t_mid = tcbp.time_at_frequency(f_mid / tcbp.geom_freq_to_Hz)
        t_max = tcbp.time_at_frequency(f_max / tcbp.geom_freq_to_Hz)
        delta_t_30_120 = (t_max - t_min)
        delta_t_30_55 = (t_mid - t_min)

        E_EM_30_120 = tcbp.E_EM_dipole_t_range(t_min, t_max) + tcbp.E_EM_quadrupole_t_range(t_min, t_max)
        E_EM_30_55 = tcbp.E_EM_dipole_t_range(t_min, t_mid) + tcbp.E_EM_quadrupole_t_range(t_min, t_mid)
        E_GW_30_120 = tcbp.E_GW_quadrupole_t_range(t_min, t_max)
        E_GW_30_55 = tcbp.E_GW_quadrupole_t_range(t_min, t_mid)

        GW_cycles_30_120 = tcbp.GW_cycles_t_range(t_min, t_max)
        GW_cycles_30_55 = tcbp.GW_cycles_t_range(t_min, t_mid)
        print(f"{s[0]} {s[1]} {delta_t_30_120:.2f} {E_GW_30_120:.2e} {E_EM_30_120:.2e} {GW_cycles_30_120:.2f} {delta_t_30_55:.2f} {E_GW_30_55:.2e} {E_EM_30_55:.2e} {GW_cycles_30_55:.2f}")
