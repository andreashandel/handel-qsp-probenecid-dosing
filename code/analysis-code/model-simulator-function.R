#this specifies a function that runs the underlying ODE model when called
#initial conditions and parameter values, as well as drug treatment information and time
#is passed in to the function
#the function returns a time-series with the solution of the ODE
# it needs: deSolve
library(deSolve)

# to test run this
#simulate_model(Ad = 0, Ac = 0, At = 0, U = 1e7, I = 0, V = 1, F = 0, A = 1, S = 0, b = 1e-8, cI = 1, k = 1e-6, p = 1e3, kF = 1e-3, cV = 10, gF = 1, hV = 1e4, Fmax = 10, cF = 5, hF = 10, gA = 1, gS = 1, cS = 2, Emax_V = 1, C50_V = 10, Emax_F = 0.5, C50_F = 10, ka = 2, Vc = 4, Vt = 20, Q = 5, Vmax = 1.6, Km = 60, fmax = 0.5, f50 = 0.12, Ad0 = 100, txstart = 1, txinterval = 0.5, txend = 4, tstart = 0, tfinal = 10, dt = 0.01, solvertype = "lsoda", tols = 1e-9)

# don't provide any defaults to ensure they are all passed in
simulate_model <- function(Ad, Ac, At, U, I, V, F, A, S, 
                           b, cI, k, p, kF, cV, gF, hV, Fmax, cF, hF, gA, gS, cS, 
                           Emax_V, C50_V, Emax_F, C50_F, ka, Vc, Vt, Q, Vmax, Km, fmax, f50, 
                           Ad0, txstart, txinterval, txend, tstart, tfinal, dt, solvertype, tols)
 {
  #inner function that specifies the ode model
  odemodel <- function(t, y, parms) {
    with(
      as.list(c(y, parms)), #lets us access variables and parameters stored in y and parms by name
      {
        # drug PK
        # dosing occurs at discrete times via the callback function in the main script
        dAd = -ka * Ad #depot
        dAc = ka * Ad - Q / Vc * Ac + Q / Vt * At - Vmax * Ac / (Km * Vc + Ac) #central compartment
        dAt = Q / Vc * Ac - Q / Vt * At #target site

        # drug PD
        Ct = At / Vt
        fu = fmax * Ct / (f50 + Ct)
        Cu = fu * Ct #drug concentration in tissue compartment
        fV = Emax_V * Cu / (C50_V + Cu) #effect of drug on virus
        fF = Emax_F * Cu / (C50_F + Cu) #effect of drug on innate
        #fV = 0.9  #effect of drug on virus
        #fF = 0.9  #effect of drug on innate
       
        
        # avoid virus to go to very low levels and then rebound
        # is biologically unreasonable since very low virus means clearance
        # this is the 'nanofox' problem when doing ODEs with low quantities
        # theoretically, <1 virus particle means extinction
        # however, we run the modeil in units of PFU/ml of virus in lung
        # unclear how that exactly translates to number of infectious virus particles
        # thus setting some low but reasonable threshold
        # see my 2007 PCB paper for more discussions
        # if (V < 0) {
        #   V <- 0
        #   dV <- 0
        # } else {
        #   dV = (1 - fV) * p * I / (1 + kF * F) - cV * V #virus
        # }
        # if (V<0) {
        #   browser()
        # }

        # define system of ODEs
        dU = -b * U * V #uninfected cells
        dI = b * U * V - cI * I - k * I * A #infected cells
        dV = (1 - fV) * p * I / (1 + kF * F) - cV * V #virus
        dF = (1 - fF) * gF * (V / (V + hV)) * (Fmax - F) - cF * F #innate response
        dA = V * F / (V * F + hF) + gA * A #adaptive response
        dS = gS * F - cS * S #symptoms

        list(c(dAd, dAc, dAt, dU, dI, dV, dF, dA, dS))
      }
    ) #close with statement
  } #end function specifying the ODEs

  #function that specifies addition of drug at the indicated time
  #drug doses, actual amount in mg used here
  #division by 50 to scale from dose per kg to mouse weight, which is about 20g
  #so 1000g/20g = 50
  adddrug <- function(t, y, parms) {
    y['Ad'] = y['Ad'] + parms['Ad0']/50
    return(y)
  }

  # main part of the function that sets up the system
  # then runs the ODE model specified above

  #combine initial conditions into a vector
  Y0 = c(Ad = Ad, Ac = Ac, At = At, U = U, I = I, V = V, F = F, A = A, S = S)
  timevec = seq(tstart, tfinal, by = dt) #vector of times for which solution is returned (not that internal timestep of the integrator is different)

  #combining parameters into a parameter vector
  odepars = c(b = b, cI = cI, k = k, p = p, kF = kF, cV = cV, gF = gF, hV = hV, Fmax = Fmax, cF = cF, hF = hF, gA = gA, gS = gS, cS = cS, 
              Emax_V = Emax_V, C50_V = C50_V, Emax_F = Emax_F, C50_F = C50_F, ka = ka, Vc = Vc, Vt = Vt, Q = Q, Vmax = Vmax, Km = Km, fmax = fmax, f50 = f50, Ad0 = Ad0)
 
  drugtimes = seq(txstart, txend, by = txinterval) #times at which drug is administered (in days)

  #this line runs the simulation, i.e. integrates the differential equations describing the infection process
  #the result is saved in the odeoutput matrix, with the 1st column the time, all other column the model variables
  #in the order they are passed into Y0 (which needs to agree with the order in virusode)
  odeoutput = deSolve::ode(
    y = Y0,
    times = timevec,
    func = odemodel,
    events = list(func = adddrug, time = drugtimes),
    parms = odepars,
    method = solvertype,
    atol = tols,
    rtol = tols
  )
  return(odeoutput)
}
