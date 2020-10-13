# equations for axial deformation of thermally expanding composite beam with fully bonded core and shell
import numpy as np


def get_shell_area(d2, d1):
    area = np.pi / 4 * (d2 ** 2 - d1 ** 2)
    return area


def get_core_area(d):
    area = np.pi /4 * d ** 2
    return area


def get_force(a_c, E_c, A_c, a_s, E_s, A_s, d_T):
    force = (a_s - a_c) * d_T / (1 / (A_c * E_c) + 1 / (A_s * E_s))
    return force


def get_mechanical_strain(force, area, E):
    strain = force / area * E
    return strain


def get_thermal_strain(alpha,delta_T):
    strain = alpha * delta_T
    return strain


def get_total_strain(thermal_strain, mechanical_strain):
    total_strain = thermal_strain + mechanical_strain
    return total_strain


def get_axial_shell_stress(a_c, E_c, A_c, a_s, E_s, A_s, d_T):
    stress = (a_c - a_s) * d_T / ((A_c / A_s) / E_s - 1 / E_c)
    return stress



def get_axial_deformation(core_metal, core_area, shell_metal, shell_area, delta_temp, bar_length):
    E_c = core_metal[0]
    a_c = core_metal[1]
    A_c = core_area
    E_s = shell_metal[0]
    a_s = shell_metal[1]
    A_s = shell_area
    d_T = delta_temp
    L = bar_length

    stress_s = -get_force(a_c, E_c, A_c, a_s, E_s, A_s, d_T)/ A_s
    m_strain_s = -get_force(a_c, E_c, A_c, a_s, E_s, A_s, d_T) / A_s / E_s
    t_strain_s = get_thermal_strain(a_s, d_T)
    total_strain_s = m_strain_s + t_strain_s

    stress_c = get_force(a_c, E_c, A_c, a_s, E_s, A_s, d_T)/ A_c
    m_strain_c = get_force(a_c, E_c, A_c, a_s, E_s, A_s, d_T) / A_c / E_c
    t_strain_c = get_thermal_strain(a_c, d_T)
    total_strain_c = m_strain_c + t_strain_c

    print("Shell")
    print(f"stress: {round(stress_s,5)}")
    print(f"mechanical strain: {round(m_strain_s,5)}")
    print(f"thermal strain: {round(t_strain_s,5)}")
    print(f"total strain: {round(total_strain_s,5)}")

    print("\nCore")
    print(f"stress: {stress_c}")
    print(f"mechanical strain: {round(m_strain_c,5)}")
    print(f"thermal strain: {round(t_strain_c,5)}")
    print(f"total strain: {round(total_strain_c,5)}")

    deformation = total_strain_c * L + L
    print(f"\nTotal Deformation Core: {round(deformation,5)}")

    deformation = total_strain_s * L + L
    print(f"Total Deformation Shell: {round(deformation,5)}")
    return


if __name__ == '__main__':

    # metals [Young's Modulus, thermal expansion coefficient]
    steel = [210 * 10 ** 9, 12 * 10 ** -6]
    aluminum = [70 * 10 ** 9, 24 * 10 ** -6]
    shell_metal = steel
    core_metal = aluminum

    # areas
    shell_area = 10000 * 10 ** -6  # m**2
    core_area = 5000 * 10 ** -6  # m**2

    # temperature
    delta_temp = 100  # degC

    # beam length
    bar_length = 1  # m

    get_axial_deformation(core_metal, core_area, shell_metal, shell_area, delta_temp, bar_length)
