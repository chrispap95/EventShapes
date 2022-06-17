import uproot_methods
import numpy as np
import matplotlib.pyplot as plt


def particleGeneratorUniformHemisphere(N=100):
    # N particles uniformly distributed on a hemisphere with r = 100
    particles_theta = 2 * np.pi * np.random.rand(N)
    particles_u = np.random.rand(N)
    particles_x = 100 * np.sqrt(1 - particles_u**2) * np.cos(particles_theta)
    particles_y = 100 * np.sqrt(1 - particles_u**2) * np.sin(particles_theta)
    particles_z = 100 * particles_u
    particles_E = np.sqrt(
        particles_x**2 + particles_y**2 + particles_z**2 + 0.13957**2
    )
    particles = uproot_methods.TLorentzVectorArray.from_cartesian(
        particles_x, particles_y, particles_z, particles_E
    )
    return particles


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

particles = particleGeneratorUniformHemisphere(500)
ax.scatter(particles.x, particles.y, particles.z)
plt.show()
