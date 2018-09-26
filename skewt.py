# makes a SkewT plot.
# writtenn by K Barfus 9/2018


def calc_dtdp_dry(T,p):
  # not applying mixed-phase model !
  # input is
  # T: temperature [K]
  R0 = 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  cp0 = specific_heat_dry_air(T)  
  dtdp = (T*R0)/(p*cp0)
  return dtdp
  
def calc_dtdp_wet(T, p, rF):
  # not applying mixed-phase model !
  # input is
  # T: temperature [K]
  # rF is liquid mixing ratio <- here 0.0 because of an irreversible process
  R0 = 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  R1 = 461.5   # gas constant for water vapour [J * kg**-1 * K**-1]
  pF1 = saturation_vapour_pressure(T) # hPa
  p0 = p - pF1                        # hPa
  rF1 = calc_rF1(pF1,p0)  # saturation mixing ratio in g/g
  lF1 = latent_heat_gas_to_liquid(T) #J/kg
  LLF1 = pF1 * lF1                   # hPa * (J/kg)
  cp0 = specific_heat_dry_air(T)      # J/(kg*K)
  cp1 = specific_heat_water_vapour(T)  # J/(kg*K)
  cp2 = specific_heat_liquid_water(T)  # J/(kg*K)
  Cp = cp0 + cp1 * rF1 + cp2 * rF  # J/(kg*K)
  v = (rF1 * lF1)/pF1 * (1.0 + (R1/R0)*rF1) * (LLF1/(R1*T**2.0))
  dtdp = ((rF1*R1*T/pF1) * (1.0 + (rF1*lF1/(R0*T))))/(Cp + v)
  return dtdp

  

def calc_rF1(pF1,p0):  # Frueh and Wirth, Eq. 4
  # input variables:
  # pF1 is saturation vapour pressure [hPa]
  # p0 is partial pressure of dry air [hPa]  
  R0 = 287.058 # gas constant for dry air [J * kg**-1 * K**-1]
  R1 = 461.5   # gas constant for water vapour [J * kg**-1 * K**-1]

  res = (R0 * pF1) / (R1 * p0)
  return res
  
def saturation_vapour_pressure(T,ice=False):
  import math  
  # calculates the saturation vapour pressure in hPa using the Clausius-Claperon equation
  # incoming variables are
  # T, temperature in [K]
  # keyword ice, indicates if even in case of temperatures lower than 273.15 K es is calculated with
  # respect to liquid water (then ice must not been set)
  # output is in hPa
  # written by K.Barfus 12/2009

  e0 = 0.611 # [kPa]
  T0 = 273.15 # [K]
  Rv = 461.0 # [J K**-1 kg**-1] gas constant for water vapour

  if(ice == True):
    if(T > 273.15):  # water
      L = 2.5 * 10.0**6.0 # J kg**-1
    else:
      L = 2.83 * 10.0**6.0  # J kg**-1
  else:
    L = 2.5 * 10.0**6.0 # J kg**-1

  es = e0 * math.exp((L/Rv)*(1.0/T0-1.0/T))
  es = es * 10.0
  return es
  
def latent_heat_gas_to_liquid(T):
  # latent heat of condensation due to Rogers and Yau in J/kg
  # valid for 248.15 K < T < 313.15 K
  # input parameters:
  T  # temperature in [K]
  T_temp = T - 273.15
  latent_heat = 2500.8 - 2.36 * T_temp + 0.0016 * T_temp**2.0 - 0.00006 * T_temp**3.0
  res = latent_heat * 1000.0
  return res
  
  # alternative approach
  # calculates the latent heat of condensation (gas -> liquid) due to
  # Fleagle, R.G. and J.A. Businger, (1980)
  # An Introduction to Atmospheric Physics.  2d ed.  Academic Press, 432 pp.
  # input
  # T in K
  # output in J kg^-1 K^-1
  #t_temp = T - 273.15
  #Lv = (25.00 - 0.02274 * t_temp) * 10.0^5.0
  
def specific_heat_dry_air(T):
  # source is unknown
  # input:: T  ![K]
  # T should be: -40°C < T < 40^C
  # output is in [J kg^-1 C^-1]

  t_temp = T - 273.15
  C_pd = 1005.60 + 0.017211 * t_temp + 0.000392 * t_temp**2.0
  return C_pd

def specific_heat_water_vapour(T):
  # due to
  # Reid, R.C., J.M. Prausnitz, and B.E. Poling (1987)
  # The Properties of Gases and Liquids.  4th ed.  McGraw-Hill, 741 pp.
  # input: T temperature [K]
  # output is in J kg^-1 K^-1
  t_temp = T - 273.15
  c_pv = 1858.0 + 3.820 * 10.0**(-1.0) * t_temp + 4.220 * 10.0**(-4.0) * t_temp**2.0 - \
   1.996 * 10.0**(-7.0) * T**3.0
  return c_pv 
  
def specific_heat_liquid_water(T):
  # input: T  ! temperature [K]
  # output is in J kg^-1 K^-1
  t_temp = T - 273.15
  c_pw =  4217.4 - 3.720283 * t_temp +0.1412855 * t_temp**2.0 - 2.654387 * 10.0**(-3.0) * t_temp**3.0 \
       + 2.093236 * 10.0**(-5.0) * t_temp**(4.0)
  return c_pw 

def vapor_pressure(pressure, mixing):
    r"""Calculate water vapor (partial) pressure.

    Given total `pressure` and water vapor `mixing` ratio, calculates the
    partial pressure of water vapor.

    Parameters
    ----------
    pressure : total atmospheric pressure
    mixing : dimensionless mass mixing ratio

    Returns
    -------
    The ambient water vapor (partial) pressure in the same units as `pressure`.

    Notes
    -----
    This function is a straightforward implementation of the equation given in many places,
    such as [Hobbs1977]_ pg.71:

    .. math:: e = p \frac{r}{r + \epsilon}

    See Also
    --------
    saturation_vapor_pressure, dewpoint

    """
    epsilon = 0.622
    return pressure * mixing / (epsilon + mixing)

def dewpoint(e):
    r"""Calculate the ambient dewpoint given the vapor pressure.

    Parameters
    ----------
    e : Water vapor partial pressure

    Returns Dew point temperature

    Notes
    -----
    This function inverts the [Bolton1980]_ formula for saturation vapor
    pressure to instead calculate the temperature. This yield the following
    formula for dewpoint in degrees Celsius:

    .. math:: T = \frac{243.5 log(e / 6.112)}{17.67 - log(e / 6.112)}

    """
    import numpy as np
    sat_pressure_0c = 6.112
    val = np.log(e / sat_pressure_0c)
    return 0.  + 243.5 * val / (17.67 - val)
    

def plot_moist_adiabats(ax, t0=None, p=None, **kwargs):
  r"""Plot moist adiabats.

  Adds saturated pseudo-adiabats (lines of constant equivalent potential
  temperature) to the plot. The default style of these lines is dashed
  blue lines with an alpha value of 0.5. These can be overridden using
  keyword arguments.

  Parameters
  ----------
  t0 : array_like, optional
     Starting temperature values in Kelvin. If none are given, they will be
     generated using the current temperature range at the bottom of
     the plot.
  p : array_like, optional
    Pressure values to be included in the moist adiabats. If not
    specified, they will be linearly distributed across the current
    plotted pressure range.
    kwargs

  """
  import numpy as np
  import math
  from matplotlib.collections import LineCollection
  # Determine set of starting temps if necessary
  if t0 is None:
    tmin, tmax = ax.get_xlim()
    t0 = np.concatenate((np.arange(tmin, 0, 15), np.arange(0,20,15), 
          np.arange(20, tmax, 5))) #* units.degC

  # Get pressure levels based on ylims if necessary
  if p is None:
    pmax, pmin = ax.get_ylim()
    nn = math.floor((pmax-pmin) /2)
    p = np.linspace(pmax, pmin, nn) #* units.mbar

        # Assemble into data for plotting
  #t = moist_lapse(p, t0[:, np.newaxis])#.to(units.degC)
  #linedata = [np.vstack((ti, p)).T for ti in t]
  npr = len(p)
  nt = len(t0)
  t_plot = np.zeros(npr)
  rF = 0.0
  segs = []
  #t_res = np.zeros(npr)     
  for it in range(0, nt):
    seg_temp = []  
    for ip in range(0, npr):
      if(ip == 0):
        t_res = t0[it]
        t_plot, p_plot = TP_to_plot(ax, t_res, p[ip], pmax)
      else:
        dt =  calc_dtdp_wet(t_res+273.15, p[ip], rF)
        t_res = t_res + (p[ip] - p[ip-1]) * dt
        t_plot, p_plot = TP_to_plot(ax, t_res, p[ip], pmax)
      seg_temp.append((t_plot, p_plot))
    #plt.plot(t_plot,p, color='b', linestyle='dashed', alpha=0.5)   
    segs.append(seg_temp) 
  line_segments = LineCollection(segs, colors='b', linestyle='dashed', alpha=0.5)
  ax.add_collection(line_segments)
               
 
def plot_dry_adiabats(ax, t0=None, p=None, **kwargs):
  r"""Plot dry adiabats.

  Adds dry adiabats (lines of constant potential temperature) to the
  plot. The default style of these lines is dashed red lines with an alpha
  value of 0.5. These can be overridden using keyword arguments.


  """
  import numpy as np
  from matplotlib.collections import LineCollection
  # Determine set of starting temps if necessary
  if t0 is None:
    xmin, xmax = ax.get_xlim()
    t0 = np.arange(xmin, (xmax + 20) + 1, 15) #* units.degC

  # Get pressure levels based on ylims if necessary
  if p is None:
    pmax, pmin = ax.get_ylim()  
    p = np.linspace(pmax, pmin, 100) #* units.mbar

  # Assemble into data for plotting
  npr = len(p)
  nt = len(t0)
  #t_plot = np.zeros(npr)
  #t_res = np.zeros(npr)     
  segs = []
  for it in range(0, nt):
    seg_temp = []  
    for ip in range(0, npr):
      if(ip == 0):
        t_res = t0[it]
        t_plot, p_plot = TP_to_plot(ax, t_res, p[ip], pmax)
      else:
        dt =  calc_dtdp_dry(t_res + 273.15, p[ip])
        #print(dt)
        t_res = t_res + (p[ip] - p[ip-1]) * dt
        t_plot, p_plot = TP_to_plot(ax, t_res, p[ip], pmax)
      seg_temp.append((t_plot, p_plot))
    segs.append(seg_temp) 
    #plt.plot(t_plot,p, color='r', linestyle='dashed', alpha=0.5)
  line_segments = LineCollection(segs, colors='r', linestyle='dashed', alpha=0.5)
  ax.add_collection(line_segments)
               

def plot_mixing_lines(ax, w=None, p=None, **kwargs):
  r"""Plot lines of constant mixing ratio.

  Adds lines of constant mixing ratio (isohumes) to the
  plot. The default style of these lines is dashed green lines with an
  alpha value of 0.8. These can be overridden using keyword arguments.

  Parameters
  ----------
  w : array_like, optional
        Unitless mixing ratio values to plot. If none are given, default
        values are used.
  p : array_like, optional
        Pressure values to be included in the isohumes. If not
        specified, they will be linearly distributed across the current
        plotted pressure range up to 600 mb.
  kwargs
      Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`

  """
  # Default mixing level values if necessary
  import numpy as np
  import matplotlib.pyplot as plt
  from matplotlib.collections import LineCollection
  if w is None:
    w = np.array([0.0004, 0.001, 0.002, 0.004, 0.007, 0.01,
                          0.016, 0.024, 0.032]).reshape(-1, 1)

  # Set pressure range if necessary
  if p is None:
    pmax, pmin = ax.get_ylim()   
    p = np.linspace(pmax, 600)

  # Assemble data for plotting
  segs = []
  for iw in range(0, len(w)):
    #td = np.zeros(len(p))
    seg_temp = []
    for ip in range(0, len(p)):
      td = np.asscalar(dewpoint(vapor_pressure(p[ip],w[iw])))
      td_plot, p_plot = TP_to_plot(ax, td, p[ip], pmax)
      seg_temp.append((td_plot,p_plot))
    #ax.semilogy(td, p, color='g', linestyle='dashed', alpha=0.8)  
    segs.append(seg_temp)
  line_segments = LineCollection(segs, colors='g', linestyle='dashed', alpha=0.8)
  ax.add_collection(line_segments)
  p_annot = 1000.0
  for iw in range(0, len(w)):
    temp = np.asscalar(w[iw]) * 1000.0
    if(temp < 1):
      annot = '{:3.1f}'.format(temp)
    else:
      annot = '{:d}'.format(int(temp))
    if(iw == len(w)-1):
      annot = annot + 'g/kg'
    td = np.asscalar(dewpoint(vapor_pressure(p_annot,w[iw])))
    td_plot, p_plot = TP_to_plot(ax, td, p_annot, pmax)
    plt.text(td_plot, p_plot, annot, color='g', horizontalalignment='center', alpha=0.8)  
    
def TP_to_plot(ax, T, P, P0):
  xy = ax.transData.transform((T,P))
  xy0 = ax.transData.transform((T,P0))
  dy = xy[1] - xy0[1]
  x_plot = xy[0] + dy
  y_plot = xy[1]
  inv = ax.transData.inverted()
  TP_plot = inv.transform((x_plot,y_plot))
  T_plot = TP_plot[0]
  P_plot = TP_plot[1]
  return T_plot, P_plot  

def plot_to_TP(ax, Tplot, Pplot, P0):
  xy = ax.transData.transform((Tplot,Pplot))
  xy0 = ax.transData.transform((Tplot,P0))
  dy = xy[1] - xy0[1]
  x_plot = xy[0] - dy
  y_plot = xy[1]
  inv = ax.transData.inverted()
  TP_plot = inv.transform((x_plot,y_plot))
  T = TP_plot[0]
  P = TP_plot[1]
  return T, P     
    
    
def skewt(pressure, temperature, dewpoint, t_ascent):
  import matplotlib.pyplot as plt
  from matplotlib.collections import LineCollection
  import numpy as np
  import math

  # Create a new figure. The dimensions here give a good aspect ratio
  pmax = 1050.0
  pmin = 100.0
  tmin = -40.0
  tmax = 50.0
  fig = plt.figure(figsize=(6.5875, 6.2125))
  ax = fig.add_subplot(111)
  ax.set_xlim([tmin, tmax])
  ax.set_ylim([pmax, pmin])
  dp = 100
  pmin_trunc = math.ceil(pmin / dp) * dp
  pmax_trunc = math.floor(pmax / dp) * dp
  n_pticks = int((pmax_trunc-pmin_trunc)/dp +1)
  pticks = np.linspace(pmin_trunc, pmax_trunc, n_pticks)
  pticks_label = [str(math.floor(pp)) for pp in pticks]
  plt.yscale('log')
  plt.yticks(pticks, pticks_label) 

  # plot grid
  xmin, xmax = ax.get_xlim()
  T_temp, P_temp = plot_to_TP(ax, tmin, pmin, pmax)
  dx = 5.0
  xmin_trunc = math.ceil(T_temp / dx) * dx
  xmax_trunc = math.floor(xmax / dx) * dx
  n_temperature0 = int((xmax_trunc - xmin_trunc)/dx +1)
  temperature0 = np.linspace(xmin_trunc, xmax_trunc, n_temperature0)
  #p = [pmin, pmax]
  segs = []
  for i_temperature in range(0, n_temperature0):
    t_plot0, p_plot0 = TP_to_plot(ax, temperature0[i_temperature], pmax, pmax)
    t_plot1, p_plot1 = TP_to_plot(ax, temperature0[i_temperature], pmin, pmax)
    segs.append([(t_plot0, p_plot0), (t_plot1, p_plot1)])    
  line_segments = LineCollection(segs, colors='k', linewidth=0.5, alpha=0.5)
  ax.add_collection(line_segments)
  
  plot_moist_adiabats(ax)
  plot_dry_adiabats(ax)
  plot_mixing_lines(ax)
    
  # plot temperature
  t_plot = []
  dew_plot = []
  t_asc_plot = []
  for i_pressure in range(0, len(t_ascent)):
    t_plot0, p_plot0 = TP_to_plot(ax, temperature[i_pressure], pressure[i_pressure], pmax)
    t_plot.append(t_plot0)
    t_plot0, p_plot0 = TP_to_plot(ax, dewpoint[i_pressure], pressure[i_pressure], pmax)  
    dew_plot.append(t_plot0)
    t_plot0, p_plot0 = TP_to_plot(ax, t_ascent[i_pressure], pressure[i_pressure], pmax)
    t_asc_plot.append(t_plot0)
  plt.plot(t_plot, pressure[0:len(t_ascent)], color='r')
  plt.plot(dew_plot, pressure[0:len(t_ascent)], color='b')
  plt.plot(t_asc_plot, pressure[0:len(t_ascent)], color='k')  
  
  plt.show()
