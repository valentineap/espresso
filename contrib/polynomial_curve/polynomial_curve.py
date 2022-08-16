from ..espresso_problem import EspressoProblem

class PolynomialCurve(EspressoProblem):
    def __init__(self, example_number=0):
        super().__init__(self,example_number)

    def suggested_model(self):


    def data(self):

    
    @abstractmethod
    def forward(self, model, with_jacobian=False):
        if self.example_number == 0:
            curveFittingFwd(model,)

    def jacobian(self, model):
        raise NotImplementedError

    def plot_model(self, model):
        raise NotImplementedError
    
    def plot_data(self, data):
        raise NotImplementedError
    
    def __getattr__(self, key):
        if key in self.params:
            return self.params[key]
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{key}'"
            )

# The following arrays define the 'sampling points' used in various examples
# In principle we could generate these on-the-fly but using hard-coded values
# guarantees repeatibility across different machines, architectures, prngs.

# The first four are uniformly-distributed points between 0 and 1
xp1 = np.array([0.02814768, 0.0885654 , 0.11748064, 0.12294414, 0.16288759,
                0.16877703, 0.1843951 , 0.39127594, 0.46038946, 0.4914394 ,
                0.54084596, 0.5641636 , 0.58900649, 0.58976937, 0.6547135 ,
                0.85518813, 0.86779304, 0.91481368, 0.951955  , 0.98545064])
xp2 = np.array([0.02061079, 0.02454563, 0.19277473, 0.2784715 , 0.29414602,
                0.32944733, 0.42601344, 0.42966824, 0.45568706, 0.4727606 ,
                0.47873683, 0.50393437, 0.59354021, 0.60901867, 0.68636357,
                0.70977542, 0.73149133, 0.78673413, 0.82707107, 0.84599549,
                0.88076936, 0.9171245 , 0.95712395, 0.97514086, 0.99942137])
xp3 = np.array([0.02129797, 0.0680587 , 0.07325879, 0.07540233, 0.0815562 ,
                0.11653215, 0.13302591, 0.17036322, 0.17948113, 0.18668185,
                0.22055489, 0.22751703, 0.26470851, 0.31516965, 0.34038959,
                0.34051367, 0.34721832, 0.35517515, 0.36644436, 0.42221368,
                0.4780765 , 0.50384201, 0.54011153, 0.56945004, 0.65217677,
                0.65762461, 0.66908502, 0.68851014, 0.71684459, 0.717921  ,
                0.72093096, 0.73367274, 0.73389493, 0.75033087, 0.75890497,
                0.76225345, 0.76552936, 0.83833901, 0.86436217, 0.88042872,
                0.90603222, 0.91094849, 0.9336106 , 0.9400528 , 0.96037445,
                0.965113  , 0.96609565, 0.97968897, 0.98044997, 0.99266562])
xp4 = np.array([0.17833636, 0.19050344, 0.26430464, 0.51520092, 0.93174146])
# This set is made up of 25 samples uniformly-distributed in [0,0.3] and a further
# 5 samples uniformly-distributed in [0.9,1]
xp5 = np.array([0.01394879, 0.02194952, 0.02925089, 0.03395566, 0.05389559,
                0.0754086 , 0.08002495, 0.08598053, 0.09475561, 0.12191561,
                0.12315494, 0.14364847, 0.18455025, 0.19247459, 0.19436732,
                0.20205877, 0.20425212, 0.21203449, 0.21238893, 0.2316263 ,
                0.2381924 , 0.24158175, 0.260489  , 0.26776801, 0.27662425,
                0.90844754, 0.92687324, 0.95818959, 0.98446575, 0.99313721])


def curveFittingFwd(model, xpts, basis='polynomial',domainLength = 1.):
    """
    A fairly generic forward model for a range of 1D curve-fitting problems.
    model        - array of model parameters
    xpts         - array of sampling points
    basis        - one of "polynomial", "fourier" or "discrete"
    domainLength - float; all valid values for xpts should be in the range [0,domainLength]
    """
    singlePoint = False
    if domainLength<=0: raise ValueError("Argument 'domainLength' must be positive")
    # Convert list inputs to arrays
    if type(model) is type([]): model = np.array(model)
    if type(xpts) is type([]): xpts = np.array(xpts)
    # Check model parameters
    try:
        nModelParameters = model.shape[0]
    except AttributeError:
        raise ValueError("Argument 'model' should be a 1-D array")
    if len(model.shape)>1: 
        raise ValueError("Argument 'model' should be a 1-D array")
    # Check x-values
    try:
        npts = xpts.shape[0]
        if len(xpts.shape)!=1: raise ValueError("Argument 'xpts' should be a 1-D array")
    except AttributeError:
        singlePoint = True
        npts = 1
    if basis == 'polynomial':
        # y = m[0] + m[1]*x + m[2]*x**2 +...
        y = model[0]*np.ones([npts])
        for iord in range(1,nModelParameters):
            y += model[iord]*xpts**iord
    elif basis == 'fourier':
        if nModelParameters%2==0: 
            raise ValueError("Fourier basis requires odd number of model parameters")
        if not np.all(0<= xpts) and np.all(xpts<= domainLength): 
            raise ValueError("For Fourier basis, all sample points must be in interval (0,domainLength)")
        # y = m[0]/2 + m[1] sin (pi x/L) + m[2] cos (pi x/L) + m[3] sin(pi x/L) + ...
        y = np.ones([npts])*0.5*model[0]
        n = 1
        for i in range(1,nModelParameters,2):
            y += model[i]*np.sin(2*n*np.pi*xpts/domainLength) + model[i+1]*np.cos(2*n*np.pi*xpts/domainLength)
            n+=1
    elif basis == 'discrete':
        if not np.all(0<= xpts) and np.all(xpts<= domainLength): 
            raise ValueError("For discrete basis, all sample points must be in interval (0,domainLength)")
        bounds = np.linspace(0,domainLength,nModelParameters+1)
        y = np.zeros([npts])
        for ipt in range(0,npts):
            y[ipt] = model[max(0,np.searchsorted(bounds,xpts[ipt])-1)]
    else:
        raise ValueError("Unrecognised  value for 'basis'; please specify one of: 'polynomial', 'fourier' or 'discrete'")
    if singlePoint:
        return y[0]
    else:
        return y

def curveFittingInv(xpts,ypts,nModelParameters, basis='polynomial',domainLength=1.,regularisation=None,priorModel=None,returnPosteriorCovariance=False):
    singlePoint=False
    if domainLength<0:raise ValueError("Argument 'domainLength' must be positive")
    if type(xpts) is type([]):xpts=np.array(xpts)
    if type(ypts) is type([]):ypts=np.array(ypts)
    try:
        npts = xpts.shape[0]
        if len(xpts.shape)!=1: raise ValueError("Argument 'xpts' should be a 1-D array")
    except AttributeError:
        singlePoint = True
        npts = 1
    try:
        if ypts.shape[0] != npts: raise ValueError("Argument 'ypts' should have same dimension as 'xpts'")
        if len(ypts.shape)!=1: raise ValueError("Argument 'ypts' should be a 1-D array")
    except AttributeError:
        if not singlePoint: raise ValueError("Argument 'ypts' should have same dimension as 'xpts'")

    if basis == 'polynomial':
        G = np.zeros([npts,nModelParameters])
        for iord in range(0,nModelParameters):
            G[:,iord] = xpts**iord
    elif basis == 'fourier':
        if nModelParameters%2==0: raise ValueError("Fourier basis requires an odd number of model parameters")
        G = np.zeros([npts,nModelParameters])
        G[:,0] = 0.5
        n = 1
        for i in range(1,nModelParameters,2):
            G[:,i] = np.sin(2*n*np.pi*xpts/domainLength)
            G[:,i+1] = np.cos(2*n*np.pi*xpts/domainLength)
            n+=1
    elif basis == 'discrete':
        G = np.zeros([npts, nModelParameters])
        bounds = np.linspace(0,domainLength,nModelParameters+1)
        y = np.zeros([npts])
        for ipt in range(0,npts):
            G[ipt,max(0,np.searchsorted(bounds,xpts[ipt])-1)] = 1.
    GTG = G.T.dot(G)
    if regularisation is not None:
        if regularisation<0: raise ValueError("Argument 'regularisation' should be positive or None")
        GTG+=regularisation*np.eye(GTG.shape[0])
        if priorModel is None:
            mp = np.zeros(nModelParameters)
        else:
            if type(priorModel) is type([]): priorModel = np.array(priorModel)
            try:
                if priorModel.shape[0]!=nModelParameters: raise ValueError ("Argument 'priorModel' should match requested number of model parameters")
                if len(priorModel.shape)!=1: raise ValueError ("Argument 'priorModel' should be a 1-D array")
            except AttributeError:
                if nModelParameters>1:raise ValueError("Argument 'priorModel' should match requested number of model parameters")
            mp = priorModel
    else:
        mp = np.zeros(nModelParameters)
    if returnPosteriorCovariance:
        return mp+np.linalg.inv(GTG).dot(G.T.dot(ypts-G.dot(mp))), np.linalg.inv(GTG)
    else:
        return mp+np.linalg.inv(GTG).dot(G.T.dot(ypts-G.dot(mp)))

def generateExampleDatasets():
    np.random.seed(42)
    # Example 1: Straight line
    model = np.array([0.5,2])
    xpts = np.random.uniform(0,1,size=20)
    xpts.sort()
    ypts = curveFittingFwd(model,xpts,'polynomial')+np.random.normal(0,0.1,size=20)
    fp = open('curve_fitting_1.dat','w')
    fp.write('# x       y        sigma\n')
    for i in range(0,20):
        fp.write("%.3f    %.3f    %.3f\n"%(xpts[i],ypts[i],0.1))
    fp.close()

    # Example 2: Cubic
    model = np.array([1.,-0.2,0.5,0.3])
    xpts = np.random.uniform(0,1,size=25)
    xpts.sort()
    ypts = curveFittingFwd(model,xpts,'polynomial')+np.random.normal(0,0.1,size=25)
    fp = open('curve_fitting_2.dat','w')
    fp.write('# x       y        sigma\n')
    for i in range(0,25):
        fp.write("%.3f    %.3f    %.3f\n"%(xpts[i],ypts[i],0.1))
    fp.close()


    # Example 3: Sinusoid
    model = np.array([1,0,0.3,-0.2,0,0,0.7,0.,0.,0.3,0.,0.,0.,0.,-.2])
    xpts = np.random.uniform(0,1,size=50)
    xpts.sort()
    ypts = curveFittingFwd(model,xpts,'fourier')+np.random.normal(0,0.1,size=50)
    fp = open('curve_fitting_3.dat','w')
    fp.write('# x       y        sigma\n')
    for i in range(0,50):
        fp.write("%.3f    %.3f    %.3f\n"%(xpts[i],ypts[i],0.1))
    fp.close()

    # Example 4: Small dataset
    model = np.array([0.,2,-.3,-.6])
    xpts = np.random.uniform(0,1,size=5)
    xpts.sort()
    ypts = curveFittingFwd(model,xpts,'polynomial')+np.random.normal(0,0.05,size=5)
    fp = open('curve_fitting_4.dat','w')
    fp.write('# x       y        sigma\n')
    for i in range(0,5):
        fp.write("%.3f    %.3f    %.3f\n"%(xpts[i],ypts[i],0.05))
    fp.close()

    # Example 5: Incomplete dataset
    model = np.array([0.3,0.3,0.,-0.2,0.5,-.8,0.1,0.125])
    xpts = np.zeros([30])
    xpts[0:25] = np.random.uniform(0,0.3,size=25)
    xpts[25:] = np.random.uniform(0.9,1.,size=5)
    xpts.sort()
    ypts = curveFittingFwd(model,xpts,'polynomial')+np.random.normal(0,0.1,size=30)
    fp = open('curve_fitting_5.dat','w')
    fp.write('# x       y        sigma\n')
    for i in range(0,30):
        fp.write("%.3f    %.3f    %.3f\n"%(xpts[i],ypts[i],0.1))
    fp.close()


    return xpts,ypts
