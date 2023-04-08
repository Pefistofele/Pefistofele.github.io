# Calculate outflow mass and energy
import astropy.units as u
import numpy as np

# Convert Jy/beam*km/s to K*km/s

bmaj, bmin = 2.20*u.arcsec, 1.76*u.arcsec
fwmh2sigma = 1. / (8 * np.log(2))**0.5

# Get solid angle of the beam
omega_B = 2*np.pi/(8*np.log(2))*bmaj*bmin

# Convert Jy/beam to K
freq = 230.538*u.GHz # rest freq of CO 2-1
# bright_T= (3.8*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))

# Calculate CO outflow mass of red and blue lobe
tau0 = 2
# tau0/(1-np.exp(-tau0))
Tex = 20
J = 2

# Get red/blue lobe maxium column density
red_tbdv_max = ((35.902*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))*u.km/u.s).value

N_red_max = (4.33e13*Tex/J**2*np.exp(2.77*J*(J+1)/Tex)*tau0/(1-np.exp(-tau0))*red_tbdv_max)*u.cm**(-2)

blue_tbdv_max = ((10.414*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))*u.km/u.s).value

N_blue_max = (4.33e13*Tex/J**2*np.exp(2.77*J*(J+1)/Tex)*tau0/(1-np.exp(-tau0))*blue_tbdv_max)*u.cm**(-2)

# Get red/blue lobe mean column density

red_tbdv_mean = ((11*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))*u.km/u.s).value

N_red_mean = (4.33e13*Tex/J**2*np.exp(2.77*J*(J+1)/Tex)*tau0/(1-np.exp(-tau0))*red_tbdv_mean)*u.cm**(-2)

blue_tbdv_mean = ((6*u.Jy/u.beam).to(u.K, u.brightness_temperature(freq, beam_area=omega_B))*u.km/u.s).value

N_blue_mean = (4.33e13*Tex/J**2*np.exp(2.77*J*(J+1)/Tex)*tau0/(1-np.exp(-tau0))*blue_tbdv_mean)*u.cm**(-2)

# Area of red/blue lobe

area_red = 1.2137e4*(0.2*1.4e3/206265)**2*u.pc**2
area_blue = 1.079e3*(0.2*1.4e3/206265)**2*u.pc**2

Mgas_red = (2.25e-16*area_red.value*N_red_mean.value)*u.Msun
Mgas_blue = (2.25e-16*area_blue.value*N_blue_mean.value)*u.Msun

# Calculate dynamic time of red/blue lobe outflow
red_size = 44.2*1.4e3*u.au
blue_size = 9.6*1.4e3*u.au

red_vmax = 35*u.km/u.s # vsys = 5km/s, vmax = 40 km/s
blue_vmax = 23*u.km/u.s # vsys = 5km/s, vmax = -18km/s

tdyn_red = (red_size/red_vmax).to(u.yr)
tdyn_blue = (blue_size/blue_vmax).to(u.yr)

# Calculate mass-loss rate
Mdot_red = Mgas_red/tdyn_red
Mdot_blue = Mgas_blue/tdyn_blue

# Calculate momentum of red/blue lobe
vrange_red = 28*u.km/u.s # 12~40 km/s
vrange_blue = 18*u.km/u.s # -18~0 km/s

p_red = Mgas_red*vrange_red
p_blue = Mgas_blue*vrange_blue

# Calculate momentum rate of red/blue lobe

pdot_red = p_red/tdyn_red
pdot_blue = p_blue/tdyn_blue

# Calculate energy of the red/blue outflow

E_red = (1/2*Mgas_red*vrange_red**2).to(u.erg)
E_blue = (1/2*Mgas_blue*vrange_blue**2).to(u.erg)

# Calculate mechanical luminosity
L_red = (E_red/tdyn_red).to(u.Lsun)
L_blue = (E_blue/tdyn_blue).to(u.Lsun)