from abc import ABC, abstractmethod
import numpy as np

class FailureCriterion(ABC):
    @abstractmethod
    def mohr(self):
        pass 

    #@abstractmethod
    #def wellbore_wall(self):
    #    pass


class MohrCoulomb(FailureCriterion):
    def __init__(self, mu, cohesion):
        self.mu = mu
        self.cohesion = cohesion
    
    def equation(self, normal_stress):
        return normal_stress * self.mu + self.cohesion

    def mohr(self, smin, smax, so=0):
        normal = np.array([smin, smax])
        shear = self.equation(normal)
        return normal, shear
        




