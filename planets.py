from constants import *


class Planets :
    def __init__(self, name, mass, x, y, z, vx, vy, vz, x_list, y_list, z_list, Ec, Ep):
        self.name = name
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.x_list = x_list
        self.y_list = y_list
        self.z_list = z_list
        self.Ec = Ec
        self.Ep = Ep
    
    def __eq__(self, other):
        return self.name == other.name

planets = []
