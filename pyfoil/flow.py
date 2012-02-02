from numpy import *
import sys
import matplotlib.pyplot as pylab
import scipy.integrate

class HSPanel():
    def __init__(self, pos0, pos1):
        self.pos0 = pos0
        self.pos1 = pos1

        delta = self.pos1 - self.pos0

        self.length = sqrt(dot(delta, delta))
        self.normal = array([-delta[1], delta[0]])/self.length
        self.tangent = delta/self.length
        self.midpoint = 0.5*(self.pos1 + self.pos0)

        self.cosTheta = delta[0]/self.length
        self.sinTheta = delta[1]/self.length

    def beta(self, pos):
        delta1 = pos - self.pos1
        delta0 = pos - self.pos0

        v1 = arctan2(delta1[1], delta1[0])
        v0 = arctan2(delta0[1], delta0[0])
        return v1 - v0

    def lnrr(self, pos):
        delta1 = pos - self.pos1
        delta0 = pos - self.pos0

        r1 = sqrt(dot(delta1, delta1))
        r0 = sqrt(dot(delta0, delta0))
        return -log(r1/r0)

class HSPanelFlow():
    def __init__(self, uinf, panels):
        self.uinf = array([float(x) for x in uinf])
        self.hsPanels = [HSPanel(p.pos0(), p.pos1()) for p in panels]
        self.setUV()

    def test(self):
        print "testing..."
        self.uinf = 0
        self.gamma = 0.0
        self.qs = array([0.0 for k in range(len(self.hsPanels))])
#        self.qs[2] = 50.0
        self.qs[1] = 50.0
#        self.qs[49] = 50.0
#        self.qs[99] = 50.0
#        self.gamma = 100.0

        self.plotField(xRange=linspace(-0.6, -0.1, 15), yRange=linspace(-0.3, 0.3, 15))

#        np = len(self.hsPanels)
#        qsg = zeros((np+1,1))
#        print self.hsPanels[0].vel(array([1,3]))

        pylab.axis('equal')
        pylab.show()
        exit(0)

    def vel(self, pos):
        betas = array([p.beta(pos) for p in self.hsPanels])
        lnrrs = array([p.lnrr(pos) for p in self.hsPanels])
        coss = array([p.cosTheta for p in self.hsPanels])
        sins = array([p.sinTheta for p in self.hsPanels])

        vxStar = (self.qs*lnrrs + self.gamma*betas)
        vyStar = (self.qs*betas - self.gamma*lnrrs)
        vx = (vxStar*coss - vyStar*sins).sum()/(2*pi)
        vy = (vxStar*sins + vyStar*coss).sum()/(2*pi)

        return array([vx, vy]) + self.uinf

    def setUV(self):
        sys.stdout.write("Setting up Hess Smith sensitivities... ")
        sys.stdout.flush()

        # (Uij, Vij) is flow on panel i from source at panel j
        np = len(self.hsPanels)
        U = zeros((np, np+1))
        V = zeros((np, np+1))

        for i,p_i in enumerate(self.hsPanels):
            for j,p_j in enumerate(self.hsPanels):
                if i == j:
                    lnrr = 0.0
                    beta = pi
                else:
                    lnrr = p_j.lnrr(p_i.midpoint)
                    beta = p_j.beta(p_i.midpoint)
                # source
                U[i,j] =  p_j.cosTheta*lnrr - p_j.sinTheta*beta
                V[i,j] =  p_j.sinTheta*lnrr + p_j.cosTheta*beta
                # vortex
                U[i,-1] += p_j.cosTheta*beta - p_j.sinTheta*(-lnrr)
                V[i,-1] += p_j.sinTheta*beta + p_j.cosTheta*(-lnrr)
        self.UV = vstack((U,V))/(2*pi)
        print "done"

    def solve(self):
        sys.stdout.write("Solving flow... ")
        sys.stdout.flush()

        np = len(self.hsPanels)

        ## flow tangency conditions
        nxs = [p.normal[0] for p in self.hsPanels]
        nys = [p.normal[1] for p in self.hsPanels]
        N = mat(hstack((diag(nxs),diag(nys))))

        ## kutta conditions
        T = zeros((1,2*np))
        
        T[0,0]    = self.hsPanels[0].tangent[0]
        T[0,np]   = self.hsPanels[0].tangent[1]
        T[0,np-1] = self.hsPanels[-1].tangent[0]
        T[0,-1]   = self.hsPanels[-1].tangent[1]

        NT = vstack((N,T))
        NTUV = dot(NT,self.UV)

        print ""
        print "rank UV:  "+str(linalg.matrix_rank(self.UV))
        print "rank NUV: "+str(linalg.matrix_rank(dot(N,self.UV)))
        print "rank NTUV:  "+str(linalg.matrix_rank(NTUV))

        uinf = vstack((mat(ones((np,1)))*self.uinf[0], mat(ones((np,1)))*self.uinf[1]))

        Nuinf = dot(N,uinf)
        NTuinf = vstack((Nuinf, array(-dot(self.uinf, self.hsPanels[0].tangent+self.hsPanels[-1].tangent))))

        qsGamma = linalg.solve(NTUV, -NTuinf)
#        qsGamma = linalg.lstsq(dot(N,self.UV), -Nuinf)[0]


        self.qs = array([qsGamma[k,0] for k in range(np)])
        self.gamma = qsGamma[-1,0]
        panelVels = dot(self.UV, qsGamma) + uinf
        for k,p in enumerate(self.hsPanels):
            p.solutionVel = array([panelVels[k,0], panelVels[k+np,0]])
#        for k,p in enumerate(self.hsPanels):
#            print ""
#            print p.solutionVel
#            print self.vel(p.midpoint + p.normal*1e-3)

        print "done"
        print "mean error: " + str(abs((dot(NTUV, qsGamma) + NTuinf)).mean())
        print "max  error: " + str(abs((dot(NTUV, qsGamma) + NTuinf)).max())


    def force(self, rho):
        sys.stdout.write("Calculating panel forces... ")
        sys.stdout.flush()
        force = array([0.0,0.0])
        Pt = 0.5*rho*dot(self.uinf, self.uinf)
        for p in self.hsPanels:
            pressure = Pt - 0.5*rho*dot(p.solutionVel,p.solutionVel)
            force += pressure*(-p.normal)*p.length
        print "done"
        return force

    def plotCps(self):
        uinf2 = dot(self.uinf, self.uinf)

        Cps = []
        xs = []
        ys = []

        for p in self.hsPanels:
            u2 = dot(p.solutionVel,p.solutionVel)
            Cps.append( 1 - u2/uinf2 )
            xs.append(p.midpoint[0])
        
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
            u2 = dot(p.solutionVel, p.solutionVel)
            pressure = Pt - 0.5*rho*u2
            Cp = 1 - u2/uinf2
#            print pressure
            force = -p.normal*Cp#pressure
            
            xs.append(p.midpoint[0])
            ys.append(p.midpoint[1])
            us.append(force[0])
            vs.append(force[1])
            
        pylab.quiver(xs,ys,us,vs)
        print "done"

    def plotField(self, thicknessYFun=lambda x: 0.0, xRange=linspace(-0.3,1.3,30), yRange=linspace(-0.4,0.4,26)):
        sys.stdout.write("Drawing flow field... ")
        sys.stdout.flush()
        xs = []
        ys = []
        us = []
        vs = []
        for x in xRange:
            for y in yRange:
                if abs(y) > thicknessYFun(x)+0.004:
                    xs.append(x)
                    ys.append(y)
                    uv = self.vel(array([x,y]))
                    us.append(uv[0])
                    vs.append(uv[1])
        pylab.quiver(xs,ys,us,vs)
        print "done"

    def plotStreamlines(self):
        sys.stdout.write("Drawing streamlines... ")
        sys.stdout.flush()
        for y0 in linspace(-0.25, 0.25, 20):
            y = scipy.integrate.odeint(lambda x,t: self.vel(x), array([-0.3,y0]), linspace(0,0.17,100),rtol=1e-6,atol=1e-6)
            pylab.plot(y[:,0],y[:,1],'g')
        print "done"

    def plotSurfaceVelocities(self):
        sys.stdout.write("Drawing surface velocities... ")
        sys.stdout.flush()
        
        xs = []
        ys = []
        us = []
        vs = []
        for p in self.hsPanels:
            xs.append(p.midpoint[0])
            ys.append(p.midpoint[1])
            us.append(p.vel[0])
            vs.append(p.vel[1])
        pylab.quiver(xs,ys,us,vs)
        print "done"

    def plotNormals(self):
        xcs = [p.midpoint[0] for p in self.hsPanels]
        ycs = [p.midpoint[1] for p in self.hsPanels]
        ucs = [p.normal[0] for p in self.hsPanels]
        vcs = [p.normal[1] for p in self.hsPanels]
        pylab.quiver(xcs,ycs,ucs,vcs,scale=15)

    def plotTangents(self):
        xcs = [p.midpoint[0] for p in self.hsPanels]
        ycs = [p.midpoint[1] for p in self.hsPanels]
        ucs = [p.tangent[0] for p in self.hsPanels]
        vcs = [p.tangent[1] for p in self.hsPanels]
        pylab.quiver(xcs,ycs,ucs,vcs,scale=15)
