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
    
    NACA_list = []
    for m in range(len(input_vals['airfoils'])): # Here you are filling 
        naca = input_vals['airfoils'][m] # naca is each type of airfoil in the 'airfoils' section of the json file
        NACA_list.append(naca)

    #############################################################################
    #  ADD TO THE END OF NACA_list BELOW IF YOU ADD ANY AIRFOILS IN THE JSON    #
    #############################################################################  

    NACA_vals = [] # this is a list that is filled with the txt file values for each airfoil. 
    alpha = (input_vals["alpha[deg]"])*math.pi/180

    # Read in the text file points
    for i in range(len(NACA_list)):
        with open(NACA_list[i], 'r') as text_handle:
            NACA_val = [list(map(float, line.strip().split())) for line in text_handle]
            NACA_vals.append(NACA_val) # this is how NACA_vals is being filled. It's being filled like this NACA_vals = [......,.......,......]
    
    x_leading_edge = 0
    x_trailing_edge = 1
    step_size = 0.01
    angle_label = str(alpha*180/math.pi) + " [deg]"

    NACA_objects = []

    for j in range(len(NACA_list)):
        NACA_ob = plot_object(NACA_vals[j], x_leading_edge, x_trailing_edge, step_size)
        NACA_objects.append(NACA_ob) 
        plot_object()

    #############################################################################
    # COMMENT THE LINES BELOW IN OR OUT TO SHOW PLOTS OF THE DIFFERENT AIRFOILS #
    #############################################################################  
    
    NACA_objects[0].vel_at_alpha_plot(alpha, angle_label)
    NACA_objects[0].airfoil_plot("NACA_2412")
    NACA_objects[1].airfoil_plot("NACA_0012")
    NACA_objects[2].airfoil_plot("NACA_2512")
    # NACA_objects[3].airfoil_plot("NACA_4512")
    # NACA_objects[4].airfoil_plot("NACA_6512")
    # NACA_objects[5].airfoil_plot("NACA_2212")
    # NACA_objects[6].airfoil_plot("NACA_2412")
    # NACA_objects[7].airfoil_plot("NACA_2612")
    # NACA_objects[8].airfoil_plot("NACA_2406")
    # NACA_objects[9].airfoil_plot("NACA_2418")

    # for point in NACA_2412_points:
    #   x, y = point
    #   print(x, y)
    plt.legend(loc = "upper right")
    plt.show()


