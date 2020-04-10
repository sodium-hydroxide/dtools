# Set of functions for dynamical systems written by Noa3h J. Blair
# Imports packages used
import numpy as np
import sympy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

### Tools for solving the driven duffing equation
# Numerically solves the driven duffing equation over one period
# q'(t) = p(t)
# p'(t) = q(t) - (q(t))^3 + g*cos(2 pi t)
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
def duffing_solution(qo,po,g):
    dt = (1.0)/50.0
    q = qo
    p = po
    for j in range(1,50):
        # uses runga kutta 4 solve equation
        t = j*dt
        k1q = dt*(p)
        k1p = dt*(q - ((q)**3) + g*np.cos(2.0*np.pi*t))
        k2q = dt*(p + (k1p/2.0))
        k2p = dt*((q + k1q) - ((q + (k1q/2.0))**3) + g*np.cos(2.0*np.pi*(t + dt/2.0)))
        k3q = dt*(p + (k2p/2.0))
        k3p = dt*((q + k2q/2.0) - ((q + k2q/2.0)**3) + g*np.cos(2.0*np.pi*(t + dt/2.0)))
        k4q = dt*(p + (k3p))
        k4p = dt*((q + k3q) - ((q + k3q)**3) + g*np.cos(2.0*np.pi*(t + dt)))
        q += (1.0/6.0)*(k1q + 2*k2q + 2*k3q + k4q)
        p += (1.0/6.0)*(k1p + 2*k2p + 2*k3p + k4p)
    return q,p

# Generates the first n values of the poincare map for the duffing equation
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
    # n:= number of periods
def duffing_poincare(qo,po,g,n):
    qvals = np.zeros(n)
    pvals = np.zeros(n)
    qvals[0] = qo
    pvals[0] = po
    for k in range(1,n):
        q,p = duffing_solution(qvals[k-1],pvals[k-1],g)
        qvals[k] = q
        pvals[k] = p
    return qvals, pvals

# Plots the poincare section
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
    # n:= number of periods
def duffing_poincare_plot(qo,po,g,n):
    qvals, pvals = duffing_poincare(qo,po,g,n)
    # plots the separatrix for undamped solution
    delta = 0.025
    Q, P = np.meshgrid(np.arange(-1.5,1.5,delta),np.arange(-1.5,1.5,delta))
    H = 0.5*(P**2) - 0.5*(Q**2) + 0.25*(Q**4)
    plt.contour(Q,P,H,0)
    plt.scatter(qvals,pvals)
    plt.show()

### Tools for driven linear oscillator
# Function for the closed form poincare map for a driven linear oscillator
# q''+w^2 p=g cos(2pi t)
    # q:= initial position
    # p:= initial momentum
    # w:= natural frequency of oscillator 
        # NEVER SET THIS EQUAL TO 2 PI
        # THE CODE WILL BREAK
    # b:= damping parameter
    # g:= driving strength
    # All these are floats
#def driven_linear_osc(q,p,w,b,g):
#    const = 1.0/((w**2 - 4.0*(np.pi**2))**2 + 4.0*(b**2))
#    cterm = np.cos( np.sqrt())
#    qnew = g*const + (q + g*const)*np.cos(w)+(p/w)*np.sin(w)
#    pnew = p*np.cos(w) + w*(const - q)*np.sin(w)
#    return qnew, pnew

# Plots the poincare map for the first n periods
#def poincare_lin_driven(q,p,w,g,n):
#    qvals = np.zeros(n)
#    pvals = np.zeros(n)
#    qvals[0] = q
#    pvals[0] = p
#    for k in range(1,n):
#        qnew, pnew = driven_linear_osc(qvals[k-1],pvals[k-1],w,g)
#        qvals[k] = qnew
#        pvals[k] = pnew
#    #return qvals, pvals
#    plt.scatter(qvals,pvals)
#    plt.show()

### Set of tools for damped driven pendulum
# Numerically solves oscillator over one period
# q'(t) = p(t)
# p'(t) = -2*b p(t) - w^2 sin(q) + g*cos(2 pi t)
    # qo:= the initial position
    # po:= the initial momentum
    # b:= damping coefficent
    # w:= natural frequency of undamped unforced oscillator
    # g:= the strength of driving force
def ddpend(qo,po,b,w,g):
    dt = (2*np.pi)/5000.0
    q = qo
    p = po
    for j in range(1,5000):
        # uses runga kutta 4 solve equation
        t = j*dt
        k1q = p
        k1p = -2.0*b*p - (w**2)*np.sin(q) + g*np.cos(2.0*np.pi*t)
        k2q = p + (k1q/2.0)
        k2p = -2.0*b*(k1p / 2.0) - (w**2)*np.sin((k1q / 2.0)) + g*np.cos(2.0*np.pi*(t + dt/2.0))
        k3q = p + (k2q/2.0)
        k3p = -2.0*b*(k2p / 2.0) - (w**2)*np.sin((k2q / 2.0)) + g*np.cos(2.0*np.pi*(t + dt/2.0))
        k4q = p + (k3q)
        k4p = -2.0*b*(k2p) - (w**2)*np.sin((k2q)) + g*np.cos(2.0*np.pi*(t + dt/2.0))
        q += (dt/6.0)*(k1q + 2*k2q + 2*k3q + k4q)
        p += (dt/6.0)*(k1p + 2*k2p + 2*k3p + k4p)
    return q,p

# Generate poincare map for damped driven pendulum over n periods
    # qo:= the initial position
    # po:= the initial momentum
    # b:= damping coefficent
    # w:= natural frequency of undamped unforced oscillator
    # g:= the strength of driving force
    # n:= number of periods
def poincare_ddpend(qo,po,b,w,g,n):
    qvals = np.zeros(n)
    pvals = np.zeros(n)
    qvals[0] = qo
    pvals[0] = po
    for k in range(1,n):
        q,p = ddpend(qvals[k-1],pvals[k-1],b,w,g)
        qvals[k] = q
        pvals[k] = p
    return qvals, pvals

# Plots poincare map for damped driven pendulum over n periods
    # qo:= the initial position
    # po:= the initial momentum
    # b:= damping coefficent
    # w:= natural frequency of undamped unforced oscillator
    # g:= the strength of driving force
    # n:= number of periods
def ddpend_poincare_plot(qo,po,b,w,g,n):
    qvals, pvals = poincare_ddpend(qo,po,b,w,g,n)
    plt.scatter(qvals,pvals)
    plt.show()

### Henon Map
def henon_map(a,b,n):
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = 1.2 
    y[0] = 0.75
    for j in range(1,n):
        x[j] = 1.0 - a*(x[j-1]**2) + y[j-1]
        y[j] = b*x[j - 1]
    return x,y

def henon_plot(a,b,n):
    x,y = henon_map(a,b,n)
    plt.plot(x,y)
    plt.title('Phase Plot of Henon Map')
    plt.show()

### Tools for Logistic map
# Interates on the logistic map x_n+1 = r x_n(1-x_n)
    # x:= current population
    # r:= growth parameter
    # n:= number of iterations to do
def logistic_iteration(x,r,n):
    pop = np.zeros(n)
    pop[0] = x
    for j in range(1,n):
        pop[j] = pop[0] + r*pop[j-1]*(1.0 - pop[j - 1])
    return pop

# Plots the time series and histogram for the logistic map
    # x:= current population
    # r:= growth parameter
    # n:= number of iterations to do
def logistic_plot(x,r,n):
    pop = logistic_iteration(x,r,n)
    time = np.linspace(0,n,num =n)
    plt.subplot(2,1,1)
    plt.plot(time,pop)
    plt.title('Time Series for Logistic Map')
    plt.subplot(2,1,2)
    num_bins = round(n/15.0)
    n, bins, patches = plt.hist(pop,num_bins)
    plt.title('Histogram for Logistic Map')
    plt.show()

### Tools for Langevin equation 
# Provides realization for 
# x_n+1 = -a*x_n + g eta_n
# eta_n is a random number given by a gaussian distribution
    # n:= number of steps to solve for
    # a:= bias towards previous values
    # g:= strength of random driving force
def langevin_sequence(n,a,g):
    x = np.zeros(n)
    x[0] = np.random.randn()
    for j in range(1,n):
        x[j] = -1.0*a*x[j - 1] + g* np.random.rand()
    return x

def langevin_plots(n,a,g):
    x = langevin_sequence(n,a,g)
    Xreal = np.real(np.fft.fft(x))
    Ximag = np.imag(np.fft.fft(x))
    # plot values of x
    plt.subplot(2,2,1)
    plt.plot(x)
    plt.title("Values of X")
    # histogram of x
    plt.subplot(2,2,2)
    num_bins = round(n/15.0)
    n, bins, patches = plt.hist(x,num_bins)
    plt.title("Distribution of X values")
    # real fourier transform
    plt.subplot(2,2,3)
    plt.plot(Xreal)
    plt.title("Real FFT(X)")
    # imag fourier transform
    plt.subplot(2,2,4)
    plt.plot(Ximag)
    plt.title("Imag FFT(X)")
    plt.show()

### Tools for stochastic population model
# provides realization of stochastic population model
# x_n+1 = (r + g*eta_n)*x_n(1-x_n) (eta is random variable given by gaussian distribution)
    # n:= number of timesteps
    # r:= mean growth rate
    # g:= strength of stochastic component
def stochastic_pop(n,r,g):
    x = np.zeros(n)


### Set of tools for a periodically impulsed particle
# Provides the next position and momentum for equation
# q_n+1 = p_n + q_n + P(q_n)
# p_n+1 = p_n + P(q_n)
# P(q) = aq+bq^2
    # qo:= initial coordinate
    # po:= initial momentum
    # a:= coefficient on linear term
    # b:= coefficeint on quadratic term
def periodic_impulse(qo,po,a,b):
    # J is the impulse
    J = a*qo + b*(qo ** 2)
    q = qo + po + J
    p = po + J
    return q,p

# Iterates the periodically impulsed system
    # q:= initial coordinate
    # p:= initial momentum
    # a:= coefficient on linear term
    # b:= coefficeint on quadratic term
    # n:= number of steps
def per_imp_interation(qo,po,a,b,n):
    qvec = np.zeros(n)
    pvec = np.zeros(n)
    qvec[0] = qo
    pvec[0] = po
    for j in range(1,n):
        q,p = periodic_impulse(qvec[j - 1],pvec[j-1],a,b)
        qvec[j] = q
        pvec[j] = p
    return qvec,pvec

# Plots the periodically impulsed system
    # q:= initial coordinate
    # p:= initial momentum
    # a:= coefficient on linear term
    # b:= coefficeint on quadratic term
    # n:= number of steps
def per_imp_plot(qo,po,a,b,n):
    qvec,pvec = per_imp_interation(qo,po,a,b,n)
    plt.plot(qvec,pvec)
    plt.show()




















##################### Ignore this
### Tools for driven van der pol oscillator
# Numerically solves oscillator over one period
# q'(t) = p(t)
# p'(t) = -p(t)(1-q(t)^2)-q(t) g*cos(2 pi t)
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
def vanderpol(qo,po,g):
    dt = (2*np.pi)/5000.0
    q = qo
    p = po
    for j in range(1,5000):
        # uses runga kutta 4 solve equation
        t = j*dt
        k1q = p
        k1p = p*(1.0-(q**2))-q- g*np.cos(2.0*np.pi*t)
        k2q = p + (k1q/2.0)
        k2p = (k1p / 2.0)*(1.0 - (k1q**2)) - k1q + g*np.cos(2.0*np.pi*(t + dt/2.0))
        k3q = p + (k2q/2.0)
        k3p = (k2p / 2.0)*(1.0 - (k2q**2)) - k2q + g*np.cos(2.0*np.pi*(t + dt/2.0))
        k4q = p + (k3q)
        k4p = k3p*(1.0-(k3q**2))-k3q- g*np.cos(2.0*np.pi*t)
        q += (dt/6.0)*(k1q + 2*k2q + 2*k3q + k4q)
        p += (dt/6.0)*(k1p + 2*k2p + 2*k3p + k4p)
    return q,p

# Generates poincare map for first n values
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
    # n:= number of periods to solve for (integer)
def poincare_vanderpol(qo,po,g,n):
    qvals = np.zeros(n)
    pvals = np.zeros(n)
    qvals[0] = qo
    pvals[0] = po
    for k in range(1,n):
        q,p = vanderpol(qvals[k-1],pvals[k-1],g)
        qvals[k] = q
        pvals[k] = p
    return qvals, pvals

# Plots poincare map for first n values
    # qo:= the initial position
    # po:= the initial momentum
    # g:= the strength of driving force
    # n:= number of periods to solve for (integer)
def vanderpol_poincare_plot(qo,po,g,n):
    qvals, pvals = poincare_vanderpol(qo,po,g,n)
    plt.scatter(qvals,pvals)
    plt.show()
