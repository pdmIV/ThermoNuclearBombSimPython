import os
import requests
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

# Set up the map
map = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)
map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color='#d9d9d9')
map.drawparallels(np.arange(-90, 90, 30))
map.drawmeridians(np.arange(-180, 180, 60))
map.drawmapboundary(fill_color='#ffffff')
map.drawmapboundary(linewidth=0.25)

# Define the parameters of the simulation
def get_parameters():
    megatonnage = float(input("Enter the megatonnage/kilotonage of the bomb: "))
    detonation_height = float(input("Enter the detonation height (in meters): "))
    wind_direction = float(input("Enter the wind direction (in degrees): "))
    temperature = float(input("Enter the ambient temperature (in degrees Celsius): "))
    return megatonnage, detonation_height, wind_direction, temperature

# Simulate the explosion
def simulate_explosion(megatonnage, detonation_height, wind_direction, temperature):
    # Calculate the blast radius based on the megatonnage of the bomb
    blast_radius = (megatonnage * 1000) ** 0.5 * 1000  # km^2

    # Calculate the pressure wave speed based on the detonation height
    pressure_wave_speed = 340.29 * (detonation_height / 1000) ** 0.5  # m/s

    # Calculate the shock wave speed based on the pressure wave speed and wind direction
    shock_wave_speed = pressure_wave_speed * np.cos(wind_direction * np.pi / 180)  # m/s

    # Calculate the blast wave speed based on the shock wave speed and blast radius
    blast_wave_speed = shock_wave_speed * blast_radius / (4 * np.pi * blast_radius)  # m/s

    # Calculate the time to reach the maximum overpressure based on the blast wave speed and distance to the center of the explosion
    time_to_max_overpressure = (blast_wave_speed * (4 * np.pi * blast_radius) / 340.29) ** 0.5  # s

    # Calculate the time to reach the peak impulse based on the time to maximum overpressure and the impulse duration
    time_to_peak_impulse = time_to_max_overpressure * 0.01  # s

    # Calculate the impulse based on the time to peak impulse and the blast wave speed
    impulse = blast_wave_speed * time_to_peak_impulse  # N*s

    # Calculate the peak overpressure based on the impulse and the time to peak overpressure
    peak_overpressure = impulse / (4 * np.pi * blast_radius * time_to_max_overpressure)  # Pa

    # Calculate the peak impulse based on the impulse and the time to peak impulse
    peak_impulse = impulse * time_to_peak_impulse  # N*s

    # Calculate the peak dynamic pressure based on the peak overpressure and the peak impulse
    peak_dynamic_pressure = peak_overpressure * peak_impulse  # Pa

    # Calculate the peak radiation pressure based on the peak dynamic pressure and the speed of light
    peak_radiation_pressure = peak_dynamic_pressure * (340.29 * (detonation_height / 1000) ** 0.5) * (1 + (0.5 * (wind_direction - 90) / 180) ** 2)  # Pa

    # Calculate the peak thermal radiation based on the peak radiation pressure and the Stefan-Boltzmann constant
    stefan_boltzmann_constant = 5.67e-8  # W/(m^2 K^4)
    peak_thermal_radiation = peak_radiation_pressure / (stefan_boltzmann_constant * (temperature + 273.15) ** 4)  # W/m^2

    return peak_overpressure, peak_impulse, peak_dynamic_pressure, peak_radiation_pressure, peak_thermal_radiation

# Get the parameters from the user
megatonnage, detonation_height, wind_direction, temperature = get_parameters()

# Simulate the explosion
peak_overpressure, peak_impulse, peak_dynamic_pressure, peak_radiation_pressure, peak_thermal_radiation = simulate_explosion(megatonnage, detonation_height, wind_direction, temperature)

# Print the results
print("Peak Overpressure:", peak_overpressure, "Pa")
print("Peak Impulse:", peak_impulse, "N*s")
print("Peak Dynamic Pressure:", peak_dynamic_pressure, "Pa")
print("Peak Radiation Pressure:", peak_radiation_pressure, "Pa")
print("Peak Thermal Radiation:", peak_thermal_radiation, "W/m^2")

# Designate the target city (e.g., New York City) with its latitude and longitude coordinates
target_city_lat = 40.7128
target_city_lon = -74.0060

# Calculate the blast radius based on the parameters of the simulation
blast_radius_km = (megatonnage * 1000) ** 0.5  # km

# Convert the blast radius from kilometers to degrees of longitude (assuming the Earth is a perfect sphere)
degrees_per_km_longitude = 360 / (40075.016686 * np.cos(np.radians(target_city_lat)))  # Approximate circumference of the Earth at the given latitude
blast_radius_deg_longitude = blast_radius_km * degrees_per_km_longitude

# Plot the target city on the map as a point
x_target, y_target = map(target_city_lon, target_city_lat)
map.plot(x_target, y_target, 'ro', markersize=8, label='Target City')

# Create a circle showing the radius of the blast
map.drawgreatcircle(target_city_lon, target_city_lat, target_city_lon + blast_radius_deg_longitude, target_city_lat, color='r', linewidth=2, label='Blast Radius')

# Add a legend to the map
plt.legend(loc='upper left')

# Finally, show the map
plt.show()
