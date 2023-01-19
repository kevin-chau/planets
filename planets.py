#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = "2022-03-19"     # simulating a solar system starting from this date
sim_duration = 100 * 365            # (int) simulation duration in days
m_sun = 1.98847e30                # Mass of the Sun in Kg
# mass_earth = 5.9722e24               # Mass of Earth in Kg
# mass_moon = 7.3477e22                # Mass of the Moon in Kg
G = 6.67408e-11                        # gravitational constant
meters_per_AU = 1.496e+11
seconds_per_day = 86400
# convert gravitational parameter from SI to AU/day^2
MU = G*m_sun*(seconds_per_day**2) / (meters_per_AU**3) # roughly 2.959e-4 AU/day^2


class Object:                     # defines an astronomical object, ie the sun, planets
    def __init__(self, name, rad, color, position, velocity):
        self.name = name
        self.r_position = np.array(position, dtype=np.cfloat) # radial position vector r
        self.velocity = np.array(velocity, dtype=np.cfloat) # velocity vector v
        self.xs = []
        self.ys = []
        self.plot = axes.scatter(self.r_position[0], self.r_position[1], color=color, s=rad**2, edgecolors=None, zorder=10)
        self.line, = axes.plot([], [], color=color, linewidth=1.4)

class SolarSystem:
    def __init__(self, thesun):
        self.thesun = thesun
        self.planets = []
        self.time = None
        self.timestamp = axes.text(.03, .94, 'Date: ', color='w', transform=axes.transAxes, fontsize='x-large')
    def add_planet(self, planet):
        self.planets.append(planet)
    def evolve(self):           # evolve the trajectories
        dt = 1.0
        self.time += dt
        plots = []
        lines = []
        for planet in self.planets:
            planet.r_position += planet.velocity * dt
            acceleration = -1*MU * planet.r_position / np.sum(planet.r_position**2)**(3./2)  # in units of AU/day^2
            planet.velocity += acceleration * dt
            planet.xs.append(planet.r_position[0])
            planet.ys.append(planet.r_position[1])
            planet.plot.set_offsets(planet.r_position[:2])
            planet.line.set_xdata(planet.xs)
            planet.line.set_ydata(planet.ys)
            plots.append(planet.plot)
            lines.append(planet.line)
        self.timestamp.set_text('Date: ' + self.time.iso)
        return plots + lines + [self.timestamp]

plt.style.use('dark_background')
fig = plt.figure(figsize=[6, 6])
axes = plt.axes([0., 0., 1., 1.], xlim=(-1.8, 1.8), ylim=(-1.8, 1.8))
axes.set_aspect('equal')
axes.axis('off')
solar_system = SolarSystem(Object("Sun", 28, 'yellow', [0, 0, 0], [0, 0, 0]))
solar_system.time = Time(sim_start_date)
solar_system.time.format='jd'
colors = ['orange', 'green', 'blue', 'red']
sizes = [0.38, 0.95, 1., 0.53]
names = ['Mercury', 'Venus', 'Earth', 'Mars']
texty = [.47, .73, 1, 1.5]
for i, nasaid in enumerate([1, 2, 3, 4]):  # The 1st, 2nd, 3rd, 4th planet in solar system
    obj = Horizons(id=nasaid, location="@sun", epochs=solar_system.time, id_type='id').vectors()
    solar_system.add_planet(Object(nasaid, 20 * sizes[i], colors[i],
                         [np.double(obj[xi]) for xi in ['x', 'y', 'z']],
                         [np.double(obj[vxi]) for vxi in ['vx', 'vy', 'vz']]))
    axes.text(0, - (texty[i] + 0.1), names[i], color=colors[i], zorder=1000, ha='center', fontsize='large')
def animate(i):
    return solar_system.evolve()
ani = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=False, interval=1,)
plt.show()
