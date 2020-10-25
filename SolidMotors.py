import numpy as np
import matplotlib.pyplot as plt
import configparser
import FMM


class Propellant:
    def __init__(self, propellant_name):
        config = configparser.ConfigParser()
        config.read('PropellantCatalog.ini')

        self.density = float(config[propellant_name]['density'])
        self.a = float(config[propellant_name]['a'])
        self.n = float(config[propellant_name]['n'])
        self.T = float(config[propellant_name]['T'])
        self.molar_mass = float(config[propellant_name]['molar_mass'])
        self.gamma = float(config[propellant_name]['gamma'])
        self.c_star = float(config[propellant_name]['c_star'])
        self.R =  float(config[propellant_name]['R'])


class Motor:
    def __init__(self, motor_name, dt, fidelity=100, burn_top=0, burn_bottom=0):

        config = configparser.ConfigParser()
        config.read('MotorCatalog.ini')

        geometry = config[motor_name]['geometry']
        propellant = config[motor_name]['propellant']
        propellant = Propellant(propellant)
        height = float(config[motor_name]['height'])
        outer_diameter = float(config[motor_name]['outer_diameter'])
        exit_area = float(config[motor_name]['exit_area'])
        throat_area = float(config[motor_name]['throat_area'])

        dt = dt
        fidelity = fidelity

        if geometry == "Bates":
            core_diameter = float(config[motor_name]['core_diameter'])
            regression_map, mc_mask, volume_ratio = FMM.build_bates_motor_3d(fidelity, height, outer_diameter, core_diameter)
            self.elapsed_time, self.thrust = FMM.burn_motor_3d(propellant, regression_map, mc_mask, volume_ratio, throat_area,
                                                    exit_area, dt)

        if geometry == "Star":
            Ri = float(config[motor_name]['bore_radius'])
            Rp = float(config[motor_name]['web_radius'])
            epsilon = float(config[motor_name]['epsilon'])
            N = int(config[motor_name]['number_of_points'])
            regression_map, map_ratio = FMM.build_star_motor_2d(fidelity, outer_diameter, Ri, Rp, epsilon, N)
            self.elapsed_time, self.thrust = FMM.burn_motor_2d(regression_map, propellant, throat_area, exit_area, dt,
                                                               map_ratio, height, Ri, outer_diameter, fidelity)

        if geometry == "WagonWheel":
            Ri = float(config[motor_name]['bore_radius'])
            Rp = float(config[motor_name]['web_radius'])
            f = float(config[motor_name]['fillet_radius'])
            epsilon = float(config[motor_name]['epsilon'])
            N = int(config[motor_name]['number_of_points'])
            theta = float(config[motor_name]["theta"])
            regression_map, map_ratio = FMM.build_wagon_wheel_motor_2d(fidelity, outer_diameter, Ri, Rp, f, epsilon, N, theta)
            self.elapsed_time, self.thrust = FMM.burn_motor_2d(regression_map, propellant, throat_area, exit_area, dt,
                                                               map_ratio, height, Ri, outer_diameter, fidelity)

        if geometry == "Finocyl":
            Ri = float(config[motor_name]['bore_radius'])
            fin_width = float(config[motor_name]['fin_width'])
            fin_length = float(config[motor_name]['fin_length'])
            N = int(config[motor_name]['number_of_fins'])
            regression_map, map_ratio = FMM.build_finocyl_motor_2d(fidelity, outer_diameter, Ri, fin_width, fin_length, N)
            self.elapsed_time, self.thrust = FMM.burn_motor_2d(regression_map, propellant, throat_area, exit_area, dt,
                                                               map_ratio, height, Ri, outer_diameter, fidelity)

        if geometry == "Custom":
            pass

    def get_time(self):
        time = self.elapsed_time

    def get_thrust(self):
        thrust = self.thrust


def plot_thrust(grain_list):
    fig, axs = plt.subplots(1, 1)
    count = 0
    for grain in grain_list:
        count = count + 1
        axs.plot(grain.elapsed_time, grain.thrust, label=f"Motor: {count}")
        axs.set_ylabel("Thrust [N]")
        axs.set_xlabel("Time [s]")
    plt.legend()
    fig.savefig("Thrust.png")
    plt.show()


if __name__ == '__main__':

    # experiment with values to change mesh fidelity and number of time steps
    fidelity = 60  # suggested ranges--- (3d FMM: 50-150), (2d FMM: 150-1000), (analytical: Not used)
    dt = 0.1  # suggested range 0.01 - 0.2

    # set motor
    # motor list: FinocylExample, WagonOneOCCAM, TubeOCCAM, StarOneOCCAM, StarTwoOCCAM, StarNickOCCAM
    motor1 = Motor("TubeOCCAM", dt, fidelity)

    # plot thrust curve for list of motors
    plot_thrust([motor1])

    thrust = motor1.get_thrust()
    time = motor1.get_time()
