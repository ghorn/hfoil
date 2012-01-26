from numpy import *
import matplotlib.pyplot as pylab

class Naca4():
    def __init__(self,m,p,t):
        # m: max camber in hundredths of chord
        # p: position of max camber in tenths of chord
        # t: max thickness in hundredths of chord
        self.m = 0.01*m
        self.p = 0.1*p
        self.t = 0.01*t
    
    def yc(self,xc):
        # xc: x/chord
        if xc < 0:
            raise ValueError("xc < 0")
        elif  xc <= self.p:
            return self.m/(self.p*self.p)*(2*self.p*xc - xc*xc)
        elif xc <= 1:
            p = self.p
            return self.m/((1-p)*(1-p))*((1-2*p) + 2*p*xc - xc*xc)
        else:
            raise ValueError("xc > 0")

    def dyc(self,xc):
        # xc: x/chord
        if xc < 0:
            raise ValueError("xc < 0")
        elif  xc <= self.p:
            return self.m/(self.p*self.p)*(2*self.p - 2*xc)
        elif xc <= 1:
            p = self.p
            return self.m/((1-p)*(1-p))*(2*p - 2*xc)
        else:
            raise ValueError("xc > 0")

    def yt(self,xc):
        t = self.t
        return 5*t*(0.2969*sqrt(xc) - 0.1260*(xc) - 0.3537*(xc)**2 + 0.2843*(xc)**3 - 0.1015*(xc)**4)

    def coords(self,xc):
        yt = self.yt(xc)

        if self.m == 0:
            yt = self.yt(xc)
            return ((xc,yt), (xc,-yt))
        
        yc = self.yc(xc)
        theta = arctan(self.dyc(xc))

        yu = yc + yt*cos(theta)
        yl = yc - yt*cos(theta)

        xu = xc - yt*sin(theta)
        xl = xc + yt*sin(theta)

        return ((xu,yu), (xl,yl))
        
if __name__=="__main__":
    foil = Naca4(2,4,12)

    xcs = linspace(0,1,500)
    coords = [foil.coords(xc) for xc in xcs]
    top    = [c[0] for c in coords]
    bottom = [c[1] for c in coords]
    pylab.plot([t[0] for t in top], [t[1] for t in top])
    pylab.plot([b[0] for b in bottom], [b[1] for b in bottom])
    pylab.axis('equal')
    pylab.show()
