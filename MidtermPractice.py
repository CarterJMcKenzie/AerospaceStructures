import BeamBending as BB
import StressStrain as SS
import ThermalExpansion as TE

print("Problem 1")
diameter = 20
thickness = 0.1
pressure = 100
torque = 502640
E = 10000 * 10 ** 3
poissons = 0.3
tank = SS.CylindricalPressureVessel(pressure, diameter, thickness, torque, poissons, E)

print("\nProblem 2")
steel = [210 * 10 ** 9, 12 * 10 ** -6]
aluminum = [70 * 10 ** 9, 24 * 10 ** -6]
shell_metal = steel
core_metal = aluminum
shell_area = 10000 * 10 ** -6  # m**2
core_area = 5000 * 10 ** -6  # m**2
delta_temp = 100  # degC
bar_length = 1  # m
TE.get_axial_deformation(core_metal, core_area, shell_metal, shell_area, delta_temp, bar_length)

print("\nProblem 3.1")
M, angle = 70.7 * 10 ** 3, 60
A = [20, 20]
X = [5, 5]
Y = [1, 7]
B = [10, 2]
H = [2, 10]
BB.beam_with_rectangular_sections(M, angle, A, X, Y, B, H)

print("\nProblem 3.2")
M, angle = 70.7 * 10 ** 3, 60
A = [1150, 1000]
X = [5, 50]
Y = [115/2, 120]
B = [10, 100]
H = [115, 10]
BB.beam_with_rectangular_sections(M, angle, A, X, Y, B, H)

print("\nProblem 4.1")


