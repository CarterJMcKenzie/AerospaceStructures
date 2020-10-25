from scipy import optimize
import numpy as np
import skfmm
from skimage import measure
import time
import matplotlib.pyplot as plt

############################################################

# burn-back equations

############################################################

def calculate_area_mach_relation(exit_mach, specific_heat, exit_area, throat_area):
    AMR = ((2 / (specific_heat + 1) * (1 + exit_mach ** 2 * (specific_heat - 1) / 2)) ** ((specific_heat + 1) /
                                             (2 * (specific_heat - 1)))) / exit_mach - exit_area / throat_area
    return AMR


def calculate_exit_mach(specific_heat, exit_area, throat_area):
    exit_mach = optimize.bisect(calculate_area_mach_relation, 1.01, 100, args=(specific_heat, exit_area, throat_area))
    exit_mach = float(exit_mach)
    return exit_mach


def calculate_nozzle_exit_temp(burn_temp, specific_heat, mach_exit_velocity):
    exit_temp = burn_temp / (1 + mach_exit_velocity ** 2 * (specific_heat - 1) / 2)
    return exit_temp


def calculate_char_velocity(specific_heat, gas_const, burn_temp):
    char_velocity = (specific_heat * gas_const * burn_temp) ** (1 / 2) /\
             specific_heat * ((2 / (specific_heat + 1)) ** ((specific_heat + 1) / (specific_heat - 1)))
    return char_velocity


def calculate_nozzle_exit_velocity(exit_mach, specific_heat, gas_const, temp_exit):
    exit_velocity = exit_mach * (specific_heat * gas_const * temp_exit) ** 0.5
    return exit_velocity


def calculate_chamber_pressure(area_burn, area_throat, burn_coef, density_propellant, c_star, burn_exp):
    chamber_pressure = (area_burn / area_throat * burn_coef * density_propellant * c_star) ** (1 / (1 - burn_exp))
    return chamber_pressure


def calculate_burn_rate(burn_coef, p_stag, burn_exp):
    burn_rate = burn_coef * p_stag ** burn_exp
    return burn_rate


def calculate_new_regression_dist(old_dist, burn_rate, dt):
    regression_dist = old_dist + burn_rate * dt
    return regression_dist


def calculate_nozzle_exit_pressure(chamber_pressure, specific_heat, mach_exit):
    exit_pressure = chamber_pressure / (1 + (specific_heat - 1) / 2 * mach_exit ** 2) ** (specific_heat / (specific_heat - 1))
    return exit_pressure


def calculate_mass_flow_rate(burn_area, burn_rate, propellant_density):
    mass_flow_rate = burn_area * burn_rate * propellant_density
    return mass_flow_rate


def calculate_vacuum_thrust(mass_flow, exit_vel, area_exit, pressure_exit):
    vacuum_thrust = mass_flow * exit_vel + area_exit * pressure_exit
    return vacuum_thrust


#######################################################################

# 2d meshing functions

#######################################################################

def make_coordinate_array_2d(motor_diameter, fidelity):
    x, y = np.meshgrid(np.linspace(-motor_diameter / 2, motor_diameter / 2, fidelity),
                                 np.linspace(-motor_diameter / 2, motor_diameter / 2, fidelity))
    return x, y


def add_perimeter_mask_2d(mesh, x, y, motor_diameter):
    mask = np.full_like(mesh, False, dtype=bool)
    mask[x ** 2 + y ** 2 > (motor_diameter/1.99) ** 2] = True
    masked_mesh = np.ma.MaskedArray(mesh, mask)
    return masked_mesh


def add_slice_mask_2d(mesh, x, y, N):
    mask = np.full_like(mesh, False, dtype=bool)
    # slope of slice border
    a = np.tan(np.pi / N)
    mask[y > a * x] = True
    mask[y < -a * x] = True
    masked_mesh = np.ma.MaskedArray(mesh, mask)
    return masked_mesh


def rotate_matrix_2d(x, y, RotRad=0):
    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],
                          [-np.sin(RotRad), np.cos(RotRad)]])
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))


def calculate_regression_2d(motor, map_ratio):
    return skfmm.distance(motor, map_ratio)


def calculate_port_area_2d(self, masked_motor):
    port_points = np.where(masked_motor == 0)
    port_area = len(port_points[0]) * self.motor_diam ** 2 / self.fidelity ** 2
    return port_area


def map_perimeter_2d(contour, fidelity, tolerance=3):
    """Returns the total length of all segments in a contour that aren't within 'tolerance' of the edge of a
    circle with diameter 'mapSize'"""
    offset = np.roll(contour.T, 1, axis=1)
    lengths = np.linalg.norm(contour.T - offset, axis=0)
    centerOffset = np.array([[fidelity / 2, fidelity / 2]])
    radius = np.linalg.norm(contour - centerOffset, axis=1)
    valid = radius < (fidelity / 2) - tolerance
    return np.sum(lengths[valid])


def calculate_burning_perimeter_2d(regression_map, regression_depth, map_ratio, fidelity):
    contours = measure.find_contours(regression_map, regression_depth, fully_connected='low')
    perimeter = 0
    for contour in contours:
        perimeter += (map_perimeter_2d(contour, fidelity)) * map_ratio
    return perimeter


def calculate_burning_area_2d(regression_map, regression_depth, motor_height, map_ratio, fidelity):
    perimeter = calculate_burning_perimeter_2d(regression_map, regression_depth, map_ratio, fidelity)
    return perimeter * motor_height


def plot_regression_2d(motor, motor_diameter, fidelity):
    """Uses the fast marching method to generate an image of how the grain regresses from the core map"""
    X, Y = make_coordinate_array_2d(motor_diameter, fidelity)
    Z = skfmm.distance(motor, motor_diameter/fidelity)
    mask = np.full_like(motor, False, dtype=bool)
    mask[X ** 2 + Y ** 2 > (motor_diameter / 2) ** 2] = True
    Z = np.ma.array(Z, mask=mask)
    fig, ax = plt.subplots()
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('2D FMM Map')
    fig.savefig('2D_FMM_Map.png')


def build_star_motor_2d(fidelity, outer_diameter, Ri, Rp, epsilon, N):
    theta = 2 * np.arctan(Rp * np.sin(np.pi * epsilon / N) * np.tan(np.pi * epsilon / N) /
                          (Rp * np.sin(np.pi * epsilon / N) - Ri * np.tan(np.pi * epsilon / N)))

    map_ratio = outer_diameter / fidelity  # in/delta
    x, y = make_coordinate_array_2d(outer_diameter, fidelity)

    motor = np.full_like(x, False, dtype=bool)

    for i in range(0, N):
        x, y = make_coordinate_array_2d(outer_diameter, fidelity)
        x, y = rotate_matrix_2d(x, y, 2 * np.pi / N * i)

        # top pointy part
        a = np.tan(theta / 2)
        h = a * Ri
        l1 = np.full_like(x, False, dtype=bool)
        l1[y < a * x - h] = True

        # bottom pointy part
        a = -np.tan(theta / 2)
        h = a * Ri
        l2 = np.full_like(x, False, dtype=bool)
        l2[y > a * x - h] = True

        # top flat part
        a = -np.tan(np.pi / 2 - np.pi / N)
        h = a * Rp / np.cos(np.pi / N)
        l3 = np.full_like(x, False, dtype=bool)
        l3[y > a * x - h] = True

        # bottom flat part
        a = np.tan(np.pi / 2 - np.pi / N)
        h = a * Rp / np.cos(np.pi / N)
        l4 = np.full_like(x, False, dtype=bool)
        l4[y < a * x - h] = True

        d1 = np.full_like(x, False, dtype=bool)
        d1[np.logical_and(l1, l2)] = True

        d2 = np.full_like(x, False, dtype=bool)
        d2[np.logical_or(l3, l4)] = True

        motor[np.logical_or(d1, d2)] = True

    motor = add_perimeter_mask_2d(motor, x, y, outer_diameter)
    regression_map = calculate_regression_2d(motor, map_ratio)
    plot_regression_2d(motor, outer_diameter, fidelity)
    return regression_map, map_ratio


def build_wagon_wheel_motor_2d(fidelity, outer_diameter, Ri, Rp, f, epsilon, N, theta):
    map_ratio = outer_diameter / fidelity  # in/delta
    x, y = make_coordinate_array_2d(outer_diameter, fidelity)

    motor = np.full_like(x, True, dtype=bool)
    motor[x ** 2 + y ** 2 < Rp ** 2] = False

    for i in range(0, N):
        x, y = make_coordinate_array_2d(outer_diameter, fidelity)
        x, y = rotate_matrix_2d(x, y, 2 * np.pi / N * i)

        # top horizontal
        h = Rp * np.sin(np.pi * epsilon / N) - f
        l1 = np.full_like(x, False, dtype=bool)
        l1[y < h] = True

        # bottom horizontal
        h = Rp * np.sin(np.pi * epsilon / N) - f
        l2 = np.full_like(x, False, dtype=bool)
        l2[y > - h] = True

        # negative slope
        a = -np.tan(theta/2)
        h = -Ri * a
        l3 = np.full_like(x, False, dtype=bool)
        l3[y > a * x + h] = True

        # positive slope
        a = np.tan(theta/2)
        h = -Ri * a
        l4 = np.full_like(x, False, dtype=bool)
        l4[y < a * x + h] = True

        d1 = np.full_like(x, False, dtype=bool)
        d1[np.logical_and(l1, l2)] = True

        d2 = np.full_like(x, False, dtype=bool)
        d2[np.logical_and(l3, l4)] = True

        dummy = np.full_like(x, False, dtype=bool)
        dummy[np.logical_and(d1, d2)] = True

        motor[np.logical_or(motor, dummy)] = True

    motor = add_perimeter_mask_2d(motor, x, y, outer_diameter)
    regression_map = calculate_regression_2d(motor, map_ratio)
    plot_regression_2d(motor, outer_diameter, fidelity)
    return regression_map, map_ratio


def build_finocyl_motor_2d(fidelity, outer_diameter, Rb, fin_width, fin_length, N):
    map_ratio = outer_diameter / fidelity  # in/delta
    x, y = make_coordinate_array_2d(outer_diameter, fidelity)

    motor = np.full_like(x, True, dtype=bool)
    motor[x ** 2 + y ** 2 < Rb ** 2] = False

    for i in range(0, N):
        x, y = make_coordinate_array_2d(outer_diameter, fidelity)
        x, y = rotate_matrix_2d(x, y, 2 * np.pi / N * i)

        # edge bore
        c1 = np.full_like(x, True, dtype=bool)
        c1[(x - fin_length) ** 2 + y ** 2 < (fin_width / 2) ** 2] = False

        motor[np.logical_or(~motor, ~c1)] = False

        # top horizontal
        h = fin_width/2
        l1 = np.full_like(x, True, dtype=bool)
        l1[y < h] = False

        # bottom horizontal
        h = fin_width/2
        l2 = np.full_like(x, True, dtype=bool)
        l2[y > - h] = False

        d1 = np.full_like(x, False, dtype=bool)
        d1[np.logical_or(l1, l2)] = True

        # back vertical
        h = 0
        l3 = np.full_like(x, True, dtype=bool)
        l3[x > h] = False

        # front vertical
        h = fin_length
        l4 = np.full_like(x, True, dtype=bool)
        l4[x < h] = False

        d2 = np.full_like(x, False, dtype=bool)
        d2[np.logical_or(l3, l4)] = True

        square = np.full_like(x, True, dtype=bool)
        square[np.logical_and(~d1, ~d2)] = False

        motor[np.logical_or(~motor, ~square)] = False

    motor = add_perimeter_mask_2d(motor, x, y, outer_diameter)
    regression_map = calculate_regression_2d(motor, map_ratio)
    plot_regression_2d(motor, outer_diameter, fidelity)
    return regression_map, map_ratio


def burn_motor_2d(regression_map, throat_area, a, density, c_star, gamma, n, exit_mach, exit_velocity, exit_area, dt,
               map_ratio, motor_height, core_diameter, motor_diameter, fidelity):
    elapsed_time = 0  # start the burn at zero seconds
    regression_depth = 0.0000001  # m
    max_regression = (motor_diameter - core_diameter) / 2
    thrust_initial = True
    count = 0

    thrust_list = np.array([0])
    elapsed_time_list = np.array([0])
    burning = True
    burning_area = calculate_burning_area_2d(regression_map, regression_depth, motor_height, map_ratio, fidelity)
    sim_start = time.time()
    while burning:
        count = count + 1

        # calculate values
        chamber_pressure = calculate_chamber_pressure(burning_area, throat_area, a, density, c_star, n)  # obsolete
        burn_rate = calculate_burn_rate(a, chamber_pressure, n)
        regression_depth = calculate_new_regression_dist(regression_depth, burn_rate, dt)
        burning_area = calculate_burning_area_2d(regression_map, regression_depth, motor_height, map_ratio, fidelity)
        m_dot = calculate_mass_flow_rate(burning_area, burn_rate, density)
        exit_pressure = calculate_nozzle_exit_pressure(chamber_pressure, gamma, exit_mach)
        thrust = calculate_vacuum_thrust(m_dot, exit_velocity, exit_area, exit_pressure)
        elapsed_time = elapsed_time + dt

        if count == 1:
            thrust_initial = thrust
        if thrust < thrust_initial * 0.1:
            burning = False
        if count >= 1000:
            burning = False
            print("2d sim took too many iterations")
        if time.time()-sim_start > 30:
            burning = False
            print("2d sim took too long")
        if regression_depth > max_regression:
            burning = False
            print(f"2d regression distance exceeded max of: {max_regression}")

        # record values
        elapsed_time_list = np.append(elapsed_time_list, elapsed_time)
        thrust_list = np.append(thrust_list, thrust)

    return elapsed_time_list, thrust_list


#######################################################################

# 3d meshing functions

#######################################################################


def make_coordinate_array_3d(motor_diameter, motor_height, fidelity, fidelity_z):
    x, y, z = np.meshgrid(np.linspace(-motor_diameter / 2, motor_diameter / 2, fidelity),
                          np.linspace(-motor_diameter / 2, motor_diameter / 2, fidelity),
                          np.linspace(-motor_height / 2, motor_height / 2, fidelity_z))
    return x, y, z


def add_uninhibited_ends_3d(mesh, fidelity_z, burn_top, burn_bottom):
    if burn_top:
        mesh[:, :, 0] = 0
    if burn_bottom:
        mesh[:, :, fidelity_z - 1] = 0
    return mesh


def add_perimeter_mask_3d(mesh, motor_diameter, fidelity_z, x, y):
    mask = np.full_like(mesh, False, dtype=bool)
    dummy = np.full_like(mesh[:, :, 0], False, dtype=bool)
    dummy[x[:, :, 0] ** 2 + y[:, :, 0] ** 2 > (motor_diameter / 1.9) ** 2] = True
    for i in range(0, fidelity_z):
        mask[:, :, i] = dummy
    masked_mesh = np.ma.MaskedArray(mesh, mask)
    return masked_mesh


def make_marching_cubes_mask_3d(mesh, motor_diameter, fidelity_z, x, y):
    mc_mask = np.full_like(mesh, True, dtype=bool)
    dummy = np.full_like(mesh[:, :, 0], True, dtype=bool)
    dummy[x[:, :, 0] ** 2 + y[:, :, 0] ** 2 > (motor_diameter / 2) ** 2] = False
    for i in range(0, fidelity_z):
        mc_mask[:, :, i] = dummy
    return mc_mask


def calculate_regression_3d(mesh):
    return skfmm.distance(mesh)


def calculate_burning_area_3d(regression_map, regression_depth, mc_mask, volume_ratio):
    regression_depth = regression_depth / volume_ratio
    verts, faces, norm, line = measure.marching_cubes(regression_map, level=regression_depth, mask=mc_mask, spacing=
                                                      (volume_ratio, volume_ratio, volume_ratio))  # (hide false values)
    burning_area = measure.mesh_surface_area(verts, faces)
    return burning_area


def build_bates_motor_3d(fidelity, height, outer_diameter, core_diameter):
    fidelity_z = int(fidelity * height / outer_diameter)  # discrete elements delta
    volume_ratio = outer_diameter / fidelity  # in/delta
    x, y, z = make_coordinate_array_3d(outer_diameter, height, fidelity, fidelity_z)

    motor = -1 * np.ones((fidelity, fidelity, fidelity_z))
    for i in range(0, fidelity_z):
        dummy = -1 * np.ones((fidelity, fidelity))
        dummy[x[:, :, 0] ** 2 + y[:, :, 0] ** 2 < (core_diameter / 2) ** 2] = 0
        motor[:, :, fidelity_z - i - 1] = dummy

    motor = add_uninhibited_ends_3d(motor, fidelity_z, False, False)
    motor = add_perimeter_mask_3d(motor, outer_diameter, fidelity_z, x, y)
    mc_mask = make_marching_cubes_mask_3d(motor, outer_diameter, fidelity_z, x, y)

    regression_map = calculate_regression_3d(motor)
    return regression_map, mc_mask, volume_ratio


def burn_motor_3d(propellant, regression_map, mc_mask, volume_ratio, throat_area, exit_area, dt):

    density = propellant.density
    a = propellant.a
    n = propellant.n
    T = propellant.T
    gamma = propellant.gamma
    c_star = propellant.c_star
    R = propellant.R

    burning_area = calculate_burning_area_3d(regression_map, -0.001, mc_mask, volume_ratio)
    exit_mach = calculate_exit_mach(gamma, exit_area, throat_area)
    exit_temp = calculate_nozzle_exit_temp(T, gamma, exit_mach)
    exit_velocity = calculate_nozzle_exit_velocity(exit_mach, gamma, R, exit_temp)

    thrust_list = np.array([0])
    elapsed_time_list = np.array([0])

    elapsed_time = 0  # start the burn at zero seconds
    regression_depth = 0.00001  # m
    thrust_initial = 10
    count = 0
    burning = True
    sim_start = time.time()
    while burning:
        count = count + 1
        chamber_pressure = calculate_chamber_pressure(burning_area, throat_area, a, density, c_star, n)  # obsolete
        burn_rate = calculate_burn_rate(a, chamber_pressure, n)
        regression_depth = calculate_new_regression_dist(regression_depth, burn_rate, dt)
        burning_area = calculate_burning_area_3d(regression_map, -regression_depth, mc_mask, volume_ratio)
        m_dot = calculate_mass_flow_rate(burning_area, burn_rate, density)
        exit_pressure = calculate_nozzle_exit_pressure(chamber_pressure, gamma, exit_mach)
        thrust = calculate_vacuum_thrust(m_dot, exit_velocity, exit_area, exit_pressure)
        elapsed_time = elapsed_time + dt

        if count == 1:
            thrust_initial = thrust

        elapsed_time_list = np.append(elapsed_time_list, elapsed_time)
        thrust_list = np.append(thrust_list, thrust)

        if thrust < thrust_initial * 0.1:
            burning = False
        if count >= 1000:
            burning = False
        if time.time()-sim_start > 60:
            burning = False
    return elapsed_time_list, thrust_list
