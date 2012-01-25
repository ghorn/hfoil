from numpy import *
import sys
import matplotlib.pyplot as pylab
import scipy.integrate

class Vortex():
    def __init__(self, position):
        self._position = array([float(p) for p in position])

    def vel(self, testPos):
        return self.velPerGamma(testPos)*self._gamma

    def velPerGamma(self, testPosition):
        dx = testPosition[0] - self._position[0]
        dy = testPosition[1] - self._position[1]
        r2 = dx*dx + dy*dy
        return array([-dy/r2, dx/r2]) / (2*pi)

class Source():
    def __init__(self, position):
        self._position = array([float(p) for p in position])

    def vel(self, testPos):
        return self.velPerGamma(testPos)*self._gamma

    def velPerGamma(self, testPosition):
        dx = testPosition[0] - self._position[0]
        dy = testPosition[1] - self._position[1]
        r2 = dx*dx + dy*dy
        return array([dx/r2, dy/r2]) / (2*pi)

class Flow():
    def __init__(self, uinf, panels):
        self.uinf = array([float(x) for x in uinf])

        # don't put vortex at trailing edge
        maxX = 0.0
        for k,p in enumerate(panels):
            if p._position0[0] > maxX:
                maxX = p._position0[0]
                maxK = k
        panels = panels[:maxK] + panels[(maxK+1):]

        self.vortices = [Vortex(panel._position0) for panel in panels]
        self.sources = [Source((0.6*maxX,0.0))]

    def primitives(self):
        return self.vortices + self.sources

    def vel(self, pos):
        return sum(vstack([p.vel(pos) for p in self.primitives()]), axis=0) + self.uinf

    def dUdGamma(self, panels):
        xMat = mat([[vs.velPerGamma(panel.centerPos())[0] for vs in self.primitives()] for panel in panels])
        yMat = mat([[vs.velPerGamma(panel.centerPos())[1] for vs in self.primitives()] for panel in panels])
        return vstack((xMat, yMat))

    def solveGammas(self, panels):
        sys.stdout.write("Solving flow... ")
        sys.stdout.flush()
        np = len(panels)
        uinf = vstack((mat(ones((np,1)))*self.uinf[0], mat(ones((np,1)))*self.uinf[1]))

        A = self.dUdGamma(panels)

        nxs = [p.length()*p.normal()[0] for p in panels]
        nys = [p.length()*p.normal()[1] for p in panels]
        N = mat(hstack((diag(nxs),diag(nys))))

        NA = dot(N,A)
#        print ""
#        print linalg.matrix_rank(A)
#        print linalg.matrix_rank(dot(N,A))

        Nuinf = dot(N,uinf)

#        piNA = dot( linalg.inv(dot(NA.T,NA)), NA.T)
#        piNA = linalg.pinv(NA)
#        gammas = dot(piNA, -Nuinf)
        gammas = linalg.lstsq(NA, -Nuinf)[0]
#        gammas = linalg.solve(NA, -Nuinf)

        for vs,gamma in zip(self.primitives(),gammas):
            vs._gamma = gamma[0,0]
        print "done"

        # print error:
        print "mean error: " + str(abs((dot(NA, gammas) + dot(N, uinf))).mean())
        print "max  error: " + str(abs((dot(NA, gammas) + dot(N, uinf))).max())

    def force(self, panels, rho):
        sys.stdout.write("Calculating panel forces... ")
        sys.stdout.flush()
        force = array([0.0,0.0])
        Pt = 0.5*rho*dot(self.uinf, self.uinf)
        for p in panels:
            vel = self.vel(p.centerPos())
            pressure = Pt - 0.5*rho*dot(vel,vel)
            force += pressure*(-p.normal())*p.length()
        print "done"
        return force

    def plotCps(self,panels):
        uinf2 = dot(self.uinf, self.uinf)

        Cps = []
        xs = []
        ys = []

        for p in panels:
            u = self.vel(p.centerPos())
            u2 = dot(u,u)
            Cps.append( 1 - u2/uinf2 )
            xs.append(p.centerPos()[0])
        
        pylab.plot(xs[:len(xs)/2],-array(Cps[:len(xs)/2]),'r')
        pylab.plot(xs[len(xs)/2:],-array(Cps[len(xs)/2:]),'g')
        print "done"

    def plotForces(self, panels, rho):
        sys.stdout.write("Drawing forces... ")
        sys.stdout.flush()

        xs = []
        ys = []
        us = []
        vs = []

        uinf2 = dot(self.uinf, self.uinf)
        Pt = 0.5*rho*uinf2

        for p in panels:
            vel = self.vel(p.centerPos())
            u2 = dot(vel, vel)
            pressure = Pt - 0.5*rho*u2
            Cp = 1 - u2/uinf2
#            print pressure
            force = -p.normal()*Cp#pressure
            
            xs.append(p.centerPos()[0])
            ys.append(p.centerPos()[1])
            us.append(force[0])
            vs.append(force[1])
            
        pylab.quiver(xs,ys,us,vs)
        print "done"

    def plot(self, thicknessYFun, xRange=linspace(-0.3,1.3,30), yRange=linspace(-0.4,0.4,26)):
        sys.stdout.write("Drawing flow... ")
        sys.stdout.flush()
        xs = []
        ys = []
        us = []
        vs = []
        for x in xRange:
            for y in yRange:
                if abs(y) > thicknessYFun(x)+0.001:
                    xs.append(x)
                    ys.append(y)
                    uv = self.vel([x,y])
                    us.append(uv[0])
                    vs.append(uv[1])
        pylab.quiver(xs,ys,us,vs)

#        xs = array(xs)
#        ys = array(ys)

#        us = array( sum([vort.vel([xs,ys[:,newaxis]])[0] for vort in self.primitives()] ) + self.uinf[0] )
#        vs = array( sum([vort.vel([xs,ys[:,newaxis]])[0] for vort in self.primitives()] ) + self.uinf[0] )
#        streamplot(xs,ys,us,vs)
        print "done"

    def plotPrimitives(self):
        pylab.plot([v._position[0] for v in self.primitives()], [v._position[1] for v in self.primitives()], 'rx')
