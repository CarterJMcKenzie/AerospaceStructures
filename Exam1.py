import math as ma
import matplotlib.pyplot as plt


def mohrs_circle_stress(sigma_x, sigma_y, tau_xy, theta, disp=False):
    oc = (sigma_x + sigma_y)/2
    radius = 0.5 * ((sigma_x-sigma_y)**2 + 4*tau_xy**2)**0.5
    sigma_1 = oc + radius
    sigma_2 = oc - radius
    tau_max = radius
    sigma_n = sigma_x*ma.cos(theta)**2 + sigma_y*ma.sin(theta)**2 + tau_xy*ma.sin(2*theta)
    tau = -((sigma_x - sigma_y)/2 * ma.sin(2*theta) - tau_xy*ma.cos(2*theta))
    beta = 180/ma.pi*ma.asin(tau_xy/radius)/2

    if disp:
        print(f"center:  {oc}")
        print(f"tau max: {tau_max}")
        print(f"sigma 1: {sigma_1}")
        print(f"sigma 2: {sigma_2}")
        print(f"\nsigma n: {sigma_n}")
        print(f"tau:     {tau}")
        print(f"\nbeta:    {beta}")
    return

def mohrs_circle_strain(epsilon_x, epsilon_y, gamma_xy, theta, disp=False):
    oc = (epsilon_x + epsilon_y) / 2
    radius = 0.5 * ((epsilon_x - epsilon_y) ** 2 + 4*gamma_xy ** 2) ** 0.5
    epsilon_1 = oc + radius
    epsilon_2 = oc - radius
    gamma_max = radius
    epsilon_n = epsilon_x * ma.cos(theta) ** 2 + epsilon_y * ma.sin(theta) ** 2 + gamma_xy * ma.sin(2 * theta)
    gamma = -((epsilon_x - epsilon_y) / 2 * ma.sin(2 * theta) - gamma_xy * ma.cos(2 * theta))


    if disp:
        print(f"center:  {oc}")
        print(f"gamma max: {gamma_max}")
        print(f"epsilon 1: {epsilon_1}")
        print(f"epsilon 2: {epsilon_2}")
        print(f"epsilon n\n: {epsilon_n}")
        print(f"gamma:     {gamma}")
    return

class CylindricalPressureVessel:

    def __init__(self, pressure, diameter, thickness):
        self.p = pressure
        self.d = diameter
        self.t = thickness

        self.hoop_stress = self.p*self.d/(self.t*2)
        self.axial_stress = self.p*self.d/(4*self.t)

    @staticmethod
    def torsion(torque, radius, moi):
        return torque*radius/moi

    @staticmethod
    def moi_cylinder(r, t):
        return 2*ma.pi*r**3*t


class ConicPressureVessel:

    def __init__(self, pressure, diameter1, diameter2, angle_degrees):
        self.p = pressure
        self.d1 = diameter1
        self.d2 = diameter2
        self.a = ma.radians(angle_degrees)

        self.yield_stress = 330 * 10**6

        self.p_max = self.yield_stress/2
        self.t = (self.p_max * self.d1)/(2*ma.cos(self.a)*(self.yield_stress))



if __name__ == "__main__":
    mohrs_circle_stress(30, -50, 30, 0, True)





# HW 2 Problem 1
# pressure = 100
# diameter = 20
# thickness = 0.1
# tank = CylindricalPressureVessel(pressure, diameter, thickness)
# sigma_x = tank.hoop_stress
# sigma_y = tank.axial_stress
# print(sigma_x)
# print(sigma_y)
# T = 502640
# J = tank.moi_cylinder(diameter/2, thickness)
# tau_xy = tank.torsion(T, diameter/2, J)
# print(tau_xy)
# gamma = 0 * ma.pi / 180
# ksi = 1000
# v = 0.3
# E = 10000 * ksi
# e_x = 1/E*(sigma_x - v*sigma_y)
# e_y = 1/E*(sigma_y - v*sigma_x)
# print(e_x)
# print(e_y)
# e_z = v/E*(sigma_y + sigma_x)
# G = E/(2*(1+v))
# gamma_xy = tau_xy / G
# print(gamma_xy)
# mohrs_circle_strain(e_x, e_y, gamma_xy, gamma_xy, True)


# in class aug 31st
# micro = 10**-6
# v = 0.3
# ep_x = 10 * micro
# ep_y = -v*ep_x
# theta = 30 * ma.pi/180
# mohrs_circle_strain(ep_x, ep_y, 0, theta, True)