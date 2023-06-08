# This is a vortex panel method in 2-D 
import json
import numpy as np
import math
import matplotlib.pyplot as plt

class plot_object:
    """This class takes any object and plots it based on x and y inputs"""
    def __init__(self, point_list, x_leading_edge, x_trailing_edge, step_size):
        self.point_list = point_list
        self.x_leading_edge = x_leading_edge
        self.x_trailing_edge = x_trailing_edge
        self.step_size = step_size

    def airfoil_plot(self, label):
        coord_list = np.array(self.point_list)
        plt.plot(coord_list[:,0], coord_list[:,1], label = label)
        plt.title("Airfoils")
        plt.xlim(self.x_leading_edge - 0.1, self.x_trailing_edge + 0.4)
        plt.ylim(-0.1, 0.14)
        plt.xlabel("x/c")
        plt.ylabel("y/c")
        plt.gca().set_aspect('equal')

    def vel_at_alpha_plot(self, alpha, angle_label):
        x_num_steps = int((self.x_trailing_edge-self.x_leading_edge)/self.step_size)
        x_values = np.linspace(self.x_leading_edge-self.x_trailing_edge/10, self.x_leading_edge, x_num_steps)
        vel_vec_list = []
        for x in x_values:
            vel_vec_list.append([x, x*math.tan(alpha)])
        vel_arrow_array = np.array(vel_vec_list)
        plt.plot(vel_arrow_array[:, 0], vel_arrow_array[:, 1], label = angle_label)

if __name__ == "__main__":
    with open("airfoils.json", 'r') as json_handle:
        input_vals = json.load(json_handle)
    
    # Extract the NACA airfoil paths
    NACA_2412_200_file = input_vals['airfoils'][0]
    NACA_0012_100_file = input_vals['airfoils'][1]
    NACA_2512_100_file = input_vals['airfoils'][2]
    NACA_4512_100_file = input_vals['airfoils'][3]
    NACA_6512_100_file = input_vals['airfoils'][4]
    NACA_2212_100_file = input_vals['airfoils'][5]
    NACA_2412_100_file = input_vals['airfoils'][6]
    NACA_2612_100_file = input_vals['airfoils'][7]
    NACA_2406_100_file = input_vals['airfoils'][8]
    NACA_2418_100_file = input_vals['airfoils'][9]


    alpha = (input_vals["alpha[deg]"])*math.pi/180

    # Read in the text file points
    with open(NACA_2412_200_file, 'r') as text_handle:
        NACA_2412_200_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_0012_100_file, 'r') as text_handle:
        NACA_0012_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_2512_100_file, 'r') as text_handle:
        NACA_2512_100_points = [list(map(float, line.strip().split())) for line in text_handle]
    
    with open(NACA_4512_100_file, 'r') as text_handle:
        NACA_4512_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_6512_100_file, 'r') as text_handle:
        NACA_6512_100_points = [list(map(float, line.strip().split())) for line in text_handle]
    
    with open(NACA_2212_100_file, 'r') as text_handle:
        NACA_2212_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_2412_100_file, 'r') as text_handle:
        NACA_2412_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_2612_100_file, 'r') as text_handle:
        NACA_2612_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_2406_100_file, 'r') as text_handle:
        NACA_2406_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    with open(NACA_2418_100_file, 'r') as text_handle:
        NACA_2418_100_points = [list(map(float, line.strip().split())) for line in text_handle]

    x_leading_edge = 0
    x_trailing_edge = 1
    step_size = 0.01
    angle_label = str(alpha*180/math.pi) + " [deg]"

    my_NACA_2412_200 = plot_object(NACA_2412_200_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_0012_100 = plot_object(NACA_0012_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2512_100 = plot_object(NACA_2512_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_4512_100 = plot_object(NACA_4512_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_6512_100 = plot_object(NACA_6512_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2212_100 = plot_object(NACA_2212_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2412_100 = plot_object(NACA_2412_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2612_100 = plot_object(NACA_2612_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2406_100 = plot_object(NACA_2406_100_points, x_leading_edge, x_trailing_edge, step_size)
    my_NACA_2418_100 = plot_object(NACA_2418_100_points, x_leading_edge, x_trailing_edge, step_size)    

        #############################################################################
        # COMMENT THE LINES BELOW IN OR OUT TO SHOW PLOTS OF THE DIFFERENT AIRFOILS #
        #############################################################################  

    my_NACA_2412_200.airfoil_plot("NACA_2412")
    my_NACA_2412_200.vel_at_alpha_plot(alpha, angle_label)
    my_NACA_0012_100.airfoil_plot("NACA_0012")
    my_NACA_2512_100.airfoil_plot("NACA_2512")
    my_NACA_4512_100.airfoil_plot("NACA_4512")
    my_NACA_6512_100.airfoil_plot("NACA_6512")
    # my_NACA_2212_100.airfoil_plot("NACA_2212")
    # my_NACA_2412_100.airfoil_plot("NACA_2412")
    # my_NACA_2612_100.airfoil_plot("NACA_2612")
    # my_NACA_2406_100.airfoil_plot("NACA_2406")
    # my_NACA_2418_100.airfoil_plot("NACA_2418")

    # for point in NACA_2412_points:
    #   x, y = point
    #   print(x, y)
    plt.legend(loc = "upper right")
    plt.show()


