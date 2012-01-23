from numpy import *
import matplotlib.pyplot as pylab
import sys

#from scipy.interpolate import interp1d
#from streamplot import streamplot

def thicknessY(tau, xc):
    # tau = thickness/chord
    # xc = x/c
    return 5*tau*(0.2969*sqrt(xc) - 0.1260*(xc) - 0.3537*(xc)**2 + 0.2843*(xc)**3 - 0.1015*(xc)**4)

def sourceU(strength, x, y):
    r2 = x*x + y*y
    return strength/(2*pi) * array([x/r2, y/r2])

class Panel():
    def __init__(self, position0, position1):
        self._position0 = array(position0)
        self._position1 = array(position1)

    def centerPos(self):
        return 0.5*(self._position0 + self._position1)

    def normal(self):
        dxy = self._position1 - self._position0
        normal = array([-dxy[1], dxy[0]])
        return normal / sqrt(dot(normal,normal))

    def length(self):
        delta = self._position0 - self._position1
        return dot(delta, delta)

class Vortex():
    def __init__(self, position):
        self._position = array(position)

    def vel(self, testPos):
        return self.velPerGamma(testPos)*self._gamma

    def velPerGamma(self, testPosition):
        dx = testPosition[0] - self._position[0]
        dy = testPosition[1] - self._position[1]
        r2 = dx*dx + dy*dy
        return array([-dy/r2, dx/r2]) / (2*pi)

class Source():
    def __init__(self, position):
        self._position = array(position)

    def vel(self, testPos):
        return self.velPerGamma(testPos)*self._gamma

    def velPerGamma(self, testPosition):
        dx = testPosition[0] - self._position[0]
        dy = testPosition[1] - self._position[1]
        r2 = dx*dx + dy*dy
        return array([dx/r2, dy/r2]) / (2*pi)


class Flow():
    def __init__(self, uinf, panels):
        self.uinf = array(uinf)

        # don't put vortex at trailing edge
        maxX = 0
        for k,p in enumerate(panels):
            if p._position0[0] > maxX:
                maxX = p._position0[0]
                maxK = k
        panels = panels[:maxK] + panels[(maxK+1):]

        self.vortices = [Vortex(panel._position0) for panel in panels]
        self.sources  = [Source((0.6*maxX,0.0))]
#        self.vortices.append(Vortex(array([0.4,0.1])))

    def primitives(self):
        return self.vortices + self.sources

    def vel(self, pos):
        return sum([v.vel(pos) for v in self.primitives()]) + self.uinf

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
        sys.stdout.write("Calculating forces... ")
        sys.stdout.flush()
        force = array([0.0,0.0])
        Pt = 0.5*rho*dot(self.uinf, self.uinf)
        for p in panels:
            vel = self.vel(p.centerPos())
            pressure = Pt - 0.5*rho*dot(vel,vel)
            force += pressure*(-p.normal())*p.length()
        print "done"
        return force


    def plot(self, thicknessYFun):
        sys.stdout.write("Drawing flow... ")
        sys.stdout.flush()
        xs = []
        ys = []
        us = []
        vs = []
        for x in linspace(-0.3,1.3,30):
            for y in linspace(-0.4,0.4, 26):
                if abs(y) > thicknessYFun(x):
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
        pylab.plot([v._position[0] for v in self.primitives()], [v._position[1] for v in self.primitives()], 'rx')
#        for v in self.primitives():
#            print v._gamma
        print "done"


class Airfoil():

    def __init__(self, tau=0.12, chord=1.0, nPanels=6):
        self.tau = tau
        self.chord = chord

        self.nPanels = nPanels
        # make number of panels even by adding one if odd
        if self.nPanels % 2 == 1:
            self.nPanels += 1
            print "setting number of panels to " + str(self.nPanels) + " to make them even"
    
        self.nChordPoints = self.nPanels/2 + 1
        #print "nChordPoints: " + str(self.nChordPoints)
        xcs = linspace(0,1,self.nChordPoints)
    
    #    xcsLeading  = xcs[0]
    #    xcsTrailing = xcs[-1]
        xcsTop      = xcs[1:-1]
        xcsBottom   = xcsTop[::-1]
    
        positions = [[0,0]] + \
                    [[chord*xc,  self.halfThickness(chord*xc)] for xc in xcsTop] + \
                    [[chord, 0]] + \
                    [[chord*xc, -self.halfThickness(chord*xc)] for xc in xcsBottom]

        self.panels = []
        for k in range(len(positions)-1):
            self.panels.append(Panel(positions[k], positions[k+1]))
        self.panels.append(Panel(positions[len(positions)-1], positions[0]))

    def halfThickness(self, x):
        if x < 0 or x > self.chord:
            return 0.0

        xc = x/self.chord
        return self.chord*thicknessY(self.tau, xc)

    def plot(self):
        xs = [p._position0[0] for p in self.panels + [self.panels[0]]]
        ys = [p._position0[1] for p in self.panels + [self.panels[0]]]
        pylab.plot(xs,ys,'-')

        xcs = [p.centerPos()[0] for p in self.panels + [self.panels[0]]]
        ycs = [p.centerPos()[1] for p in self.panels + [self.panels[0]]]
        ucs = [p.normal()[0] for p in self.panels + [self.panels[0]]]
        vcs = [p.normal()[1] for p in self.panels + [self.panels[0]]]
#        pylab.quiver(xcs,ycs,ucs,vcs,scale=15)

normUinf = 20
alpha = 5*pi/180
rho = 1.2
uinf = array( [cos(alpha)*normUinf, sin(alpha)*normUinf] )

foil = Airfoil(tau=0.12, chord=1.0, nPanels=200)
flow = Flow(uinf, foil.panels)

flow.solveGammas(foil.panels)


# calculate forces by "integrating" pressure
forces = flow.force(foil.panels, rho)
forcesWind = dot(matrix([[cos(alpha), sin(alpha)],[-sin(alpha), cos(alpha)]]), forces)
print ""
print "forces calculated from pressure:"
print "body frame force:              " + str(forces)
print "wind frame force (drag, lift): " + str(forcesWind)
print "CL: " + str(forcesWind[0,1]/(0.5*rho*normUinf*normUinf))
print "CD: " + str(forcesWind[0,0]/(0.5*rho*normUinf*normUinf))

# calculate forces from "vorticity" pressure
print ""
print "forces calculated from vorticity:"
vorticity = sum([v._gamma for v in flow.vortices])
print "vorticity: " +str(vorticity)
vortForce = rho*cross(array( [cos(alpha)*normUinf, sin(alpha)*normUinf, 0] ), array([ 0, 0, vorticity]))
vortForceWind = dot(matrix([[cos(alpha), sin(alpha),0],[-sin(alpha), cos(alpha),0],[0,0,1]]), vortForce)
print "body frame force:                " + str(vortForce)
print "wind frame force (drag, lift,_): " + str(vortForceWind)
print "CL: " + str(vortForceWind[0,1]/(0.5*rho*normUinf*normUinf))
print "CD: " + str(vortForceWind[0,0]/(0.5*rho*normUinf*normUinf))
foil.plot()
flow.plot(lambda x: foil.halfThickness(x))

pylab.axis('equal')
pylab.show()


#xcs = linspace(0,1,100)
#T = thicknessY(0.12,xcs)
#
#pylab.plot(xcs,T)
#pylab.plot(xcs,-T)
#pylab.axis('equal')
#pylab.show()
