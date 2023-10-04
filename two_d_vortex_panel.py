import json
import numpy as np
from tabulate import tabulate

np.set_printoptions(precision=15)

class vortex_panels:
    """This class contains functions that calculates position of nodes, control points, li, xi, eta, phi, psi, P_matrix, A_matrix, gamma_vector, cartesian_velocity, C_p, C_L, C_mle, C_mc/4"""
    def __init__(self, json_file):
        self.json_file = json_file
        self.load_json()

    def load_json(self):
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.airfoils = input_vals['airfoils']
            self.alpha_deg = np.radians(input_vals["alpha[deg]"])
            self.vel_inf = input_vals["freestream_velocity"]

    def pull_nodes(self, i):
        foil = self.airfoils[i]
        with open(foil, 'r') as text_handle:
            points = [list(map(float, line.strip().split())) for line in text_handle]
        self.geometry = np.array(points)
    
    def pull_label(self, i):
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            foil = input_vals['airfoils'][i]  # index 10 is for the 10 point airfoil, use for testing.
        self.label = foil

    def pull_length(self):
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)
            alpha = np.radians(input_vals["alpha[deg]"]) ## INPUT_alpha ##
            length = len(alpha)
        self.length = length
    
    def pull_alpha(self, j):
        self.alpha = self.alpha_deg[j]
    
    def calc_control_points(self):
        """Given a list of nodes, this function calculates the control points by taking the average position of each pair adjacent nodes and returns them in a Nx2 list"""
        NACA_list = self.geometry
        control_points = []
        for i in range(0, len(NACA_list)-1):
            x1, y1 = NACA_list[i]
            x2, y2 = NACA_list[i+1]
            point = [(x1 + x2)/2, (y1 + y2)/2]
            control_points.append(point)
        control_points_array = np.array(control_points)
        self.control_points = control_points_array
   
    def calc_L(self):
        """Given a list of nodes, this function calculates each value of L_j and returns them in an 1xN list"""
        points_list = self.geometry
        x, y = points_list[:, 0], points_list[:, 1]
        diff_x = np.diff(x)
        diff_y = np.diff(y)
        L_vals = np.sqrt(diff_x**2 + diff_y**2)
        self.L = L_vals
   
    def calc_xi_eta_phi_psi(self, control_x, control_y, point_x_2, point_x_1, point_y_2, point_y_1, L):
        """This function calculates the xi, eta, phi, and psi values given a list of nodes and control points"""
        xi = (1/L)*((point_x_2-point_x_1)*(control_x-point_x_1)+(point_y_2-point_y_1)*(control_y-point_y_1)) # this is the dot product for solving for xi
        eta = (1/L)*(-(point_y_2-point_y_1)*(control_x-point_x_1)+(point_x_2-point_x_1)*(control_y-point_y_1)) # this is the dot product for solving for eta
        phi = np.arctan2(eta*L, eta**2+xi**2-xi*L)
        psi = 0.5*np.log((xi**2+eta**2)/((xi-L)**2+eta**2))
        rotate_xi_eta = np.array([[(L-xi)*phi+(eta*psi), (xi*phi)-(eta*psi)], [(eta*phi)-((L-xi)*psi)-L, (-eta*phi)-(xi*psi)+L]])
        return rotate_xi_eta
   
    def calc_p_matrix(self, mat_1, mat_2, L):
        """This function calculates the p_matrix given the two rotation matrices found in the calc_a_matrix function"""
        p_matrix = (1/(2*np.pi*(L**2)))*np.matmul(mat_1,mat_2)
        return p_matrix
        
    def calc_a_matrix(self):
        """This function finds the nxn a matrix given a list of nodes, control points, and correct functions that calculate xi,eta,phi,psi,li, and lj"""
        points_list = self.geometry
        n = int(len(points_list))
        x, y = points_list[:,0], points_list[:,1]
        x_control = self.control_points[:,0]
        y_control = self.control_points[:,1]
        a_vals = np.zeros((n, n))  # Initialize an empty array
        for i in range(0,n-1):
            for j in range(0,n-1):              
                l_j = self.L[j]                
                l_i = self.L[i]

                # define rotation matrix for x and y in the p_matrix calculation
                rotate_x_y = np.array([[x[j+1]-x[j], -(y[j+1]-y[j])],[y[j+1]-y[j], x[j+1]-x[j]]])

                # define rotation matrix for xi and eta
                rotate_xi_eta = self.calc_xi_eta_phi_psi(x_control[i], y_control[i], x[j+1], x[j], y[j+1], y[j], l_j)

                # Calculate the P matrix at i, j
                p_matrix = self.calc_p_matrix(rotate_x_y, rotate_xi_eta, l_j)

                a_vals[i,j] = a_vals[i,j] + ((x[i+1]-x[i])*p_matrix[1,0]-(y[i+1]-y[i])*p_matrix[0,0])/l_i
                a_vals[i,j+1] = a_vals[i,j+1] + ((x[i+1]-x[i])*p_matrix[1,1]-(y[i+1]-y[i])*p_matrix[0,1])/l_i    
        a_vals[n-1,0] = 1.0
        a_vals[n-1,n-1] = 1.0 
        self.A_matrix = a_vals
   
    def calc_B_matrix(self): # this matrix is the Nx1 matrix that's in equation 4.32 in the Aeronautics engineering handbook
        """This function finds the d_matrix given the a_matrix and vel_inf"""
        points_list = self.geometry
        n = int(len(points_list))
        x, y = points_list[:,0], points_list[:,1]
        B_matrix = np.zeros(n)
        for i in range(0,n-1):
            diff_x = x[i+1]-x[i]
            diff_y = y[i+1]-y[i]
            l_val = self.L[i]
            B_val = ((diff_y*np.cos(self.alpha))-(diff_x*np.sin(self.alpha)))/l_val
            B_matrix[i] = B_val
        self.B_matrix = B_matrix
   
    def calc_gammas(self):
        """This function finds the gamma values given matrix_a, matrix_d and vel_inf"""
        gammas = np.linalg.solve(self.A_matrix, self.vel_inf*self.B_matrix)
        self.gammas = gammas
   
    def calc_CL(self):
        """This function finds CL given gammas, vel_inf, a geometry, and l_i"""
        gammas = self.gammas
        vel_inf = self.vel_inf
        points_list = self.geometry
        l_i = self.L
        n = int(len(points_list))
        Coeff_L = 0.0
        for i in range(0, n-1):
            Co_L = l_i[i]*((gammas[i]+gammas[i+1])/(vel_inf))
            Coeff_L += Co_L
        self.Coeff_L = Coeff_L
   
    def calc_Cm_le(self):
        """This function finds the moment coefficient calculated at the leading edge"""
        points_list = self.geometry
        gammas = self.gammas
        l_i = self.L
        n = int(len(points_list))
        x, y = points_list[:,0], points_list[:,1]
        Cm_le = 0
        for i in range(0,n-1):
            cos_coeff = (2*x[i]*gammas[i]+x[i]*gammas[i+1]+x[i+1]*gammas[i]+2*x[i+1]*gammas[i+1])
            sin_coeff = (2*y[i]*gammas[i]+y[i]*gammas[i+1]+y[i+1]*gammas[i]+2*y[i+1]*gammas[i+1])
            Cmle = l_i[i]*(cos_coeff*np.cos(self.alpha)+sin_coeff*np.sin(self.alpha))
            Cm_le = Cm_le + Cmle
        Cm_le = (-1/3)*Cm_le/self.vel_inf
        self.Coeff_mle = Cm_le
   
    def calc_Cm_4(self):
        """This function calculates the coefficient of lift at the quarter chord"""
        Cm_4 = self.Coeff_mle+ 0.25*self.Coeff_L*np.cos(self.alpha)
        self.Cm_4 = Cm_4
    
    def program(self, i):
        """This is a run function that uses all of the functions above."""
        self.pull_nodes(i)
        self.pull_length()
        self.load_json()
        points = []
        for j in range(self.length):
            self.pull_alpha(j)
            self.calc_L()
            self.calc_control_points()
            self.calc_a_matrix()
            self.calc_B_matrix()
            self.calc_gammas()
            self.calc_CL()
            self.calc_Cm_le()
            self.calc_Cm_4()
            data = [np.degrees(self.alpha), self.Coeff_L, self.Coeff_mle, self.Cm_4]
            points.append(data)
        self.pull_label(i)
        Label_modified = self.label.replace(".txt", " ")
        print("\n                   ", Label_modified)
        print(tabulate(points, headers=["alpha(deg)", 'CL', "Cm_le", 'Cm_c/4']), "\n")
        


   