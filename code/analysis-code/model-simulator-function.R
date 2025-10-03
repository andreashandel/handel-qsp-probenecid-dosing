#this specifies a function that runs the underlying ODE model when called
#initial conditions and parameter values, as well as drug treatment information and time
#is passed in to the function
#the function returns a time-series with the solution of the ODE
# required R packages for the function are loaded in the main script
# it needs: deSolve

simulate_model <- function(
  Ad = 0,
  Ac = 0,
  At = 0,
  U = 1e7,
  I = 0,
  V = 1,
  F = 0,
  A = 1,
  S = 0,
  b = 1e-8,
  cI = 1,
  k = 1e-6,
  p = 1e3,
  kF = 1e-3,
  cV = 10,
  gF = 1,
  hV = 1e4,
  Fmax = 10,
  cF = 5,
  hF = 10,
  gA = 1,
  gS = 1,
  cS = 2,
  Emax_V = 1,
  C50_V = 10,
  Emax_F = 0.5,
  C50_F = 10,
  ka = 2,
  Vc = 4,
  Vt = 20,
  Q = 5,
  Vmax = 1.6,
  Km = 60,
  fmax = 0.5,
  f50 = 0.12,
  txstart = 1,
  txinterval = 0.5,
  txend = 4,
  Ad0 = 100,
  tstart = 0,
  tfinal = 10,
  dt = 0.01
) {
  #function that specifies the ode model
  odemodel <- function(t, y, parms) {
    with(
      as.list(c(y, parms)), #lets us access variables and parameters stored in y and parms by name
      {
        # drug PK
        # dosing occurs at discrete times via the callback function in the main script
        dAd = -ka * Ad #depot
        dAc = ka * Ad - Q / Vc * Ac + Q / Vt * At - Vmax * Ac / (Km * Vc + Ac) #central compartment
        dAt = Q / Vc * Ac - Q / Vt * At #taret site

        # drug PD
        Ct = At / Vt
        fu = fmax * Ct / (f50 + Ct)
        Cu = fu * Ct #drug concentration in tissue compartment
        f_V = Emax_V * Cu / (C50_V + Cu) #effect of drug on virus
        f_F = Emax_F * Cu / (C50_F + Cu) #effect of drug on innate

        # define system of ODEs
        dU = -b * U * V #uninfected cells
        dI = b * U * V - cI * I - k * I * A #infected cells
        dV = (1 - f_V) * p * I / (1 + kF * F) - cV * V #virus
        dF = (1 - f_F) * gF * (V / (V + hV)) * (Fmax - F) - cF * F #innate response
        dA = V * F / (V * F + hF) + gA * A #adaptive response
        dS = gS * F - cS * S #symptoms

        list(c(dAd, dAc, dAt, dU, dI, dV, dF, dA, dS))
      }
    ) #close with statement
  } #end function specifying the ODEs

  #function that specifies addition of drug at the indicated time
  adddrug <- function(t, y, parms) {
    y['Ad'] = y['Ad'] + parms['Ad0']
    return(y)
  }

  # main part of the function that sets up the system
  # then runs the ODE model specified above

  #combine initial conditions into a vector
  Y0 = c(Ad = Ad, Ac = Ac, At = At, U = U, I = I, V = V, F = F, A = A, S = S)
  timevec = seq(tstart, tfinal, by = dt) #vector of times for which solution is returned (not that internal timestep of the integrator is different)

  #combining parameters into a parameter vector
  pars = c(
    b = b,
    cI = cI,
    k = k,
    p = p,
    kF = kF,
    cV = cV,
    gF = gF,
    hV = hV,
    Fmax = Fmax,
    cF = cF,
    hF = hF,
    gA = gA,
    gS = gS,
    cS = cS,
    Emax_V = Emax_V,
    C50_V = C50_V,
    Emax_F = Emax_F,
    C50_F = C50_F,
    ka = ka,
    Vc = Vc,
    Vt = Vt,
    Q = Q,
    Vmax = Vmax,
    Km = Km,
    fmax = fmax,
    f50 = f50,
    Ad0 = Ad0
  )

  drugtimes = seq(txstart, txend, by = txinterval) #times at which drug is administered (in days)

  #this line runs the simulation, i.e. integrates the differential equations describing the infection process
  #the result is saved in the odeoutput matrix, with the 1st column the time, all other column the model variables
  #in the order they are passed into Y0 (which needs to agree with the order in virusode)
  odeoutput = deSolve::ode(
    y = Y0,
    times = timevec,
    func = odemodel,
    events = list(func = adddrug, time = drugtimes),
    parms = pars,
    atol = 1e-8,
    rtol = 1e-8
  )

  return(odeoutput)
}
