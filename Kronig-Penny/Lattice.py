import numpy as np

class Rectangle:
    def __init__(self,ax,ay,theta=np.pi/2,space="2d"):
        labels = ["G", "X", "M", "XP"]
        x, y = 1.0/ax,1.0/ay
        a, b = np.cos(theta), np.sin(theta)
        dx, yp = y*a, y*x
        points = {
            "G" : np.array([0, 0]),
            "X" : np.array([x, 0])*np.pi,
            "M" : np.array([x+dx, yp])*np.pi,
            "XP": np.array([dx, yp])*np.pi
        }
        self.points = points
        self.labels = labels
        
    def trayectory(self,targets,n=100):
        T = len(targets)
        trip = np.empty(( n*(T-1), 2 ))
        for i in xrange(1,T):
            A,B = targets[i-1],targets[i]
            x0,y0 = self.points[A]
            x1,y1 = self.points[B]
            trip[n*(i-1):n*i,0] = np.linspace(x0,x1,n,endpoint=True)
            trip[n*(i-1):n*i,1] = np.linspace(y0,y1,n,endpoint=True)
        return trip      
