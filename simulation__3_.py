import numpy as np
from sevnpy.sevn import SEVNmanager, Star
from fireworks.nbodylib.nunits import Nbody_units as NU
from fireworks.nbodylib.potentials import MultiPotential, NFW, MyamotoNagai, TruncatedPlaw
from fireworks.particles import Particles
import pandas as pd
import warnings

warnings.filterwarnings("ignore")


class MilkyWayPotential(MultiPotential):
    Mscale = 1e10
    rscale = 1000
    NU = NU(M=Mscale, L=rscale)

    def __init__(self):
        Mvir = 0.8e12 / self.Mscale
        c = 15.3
        rs = 16
        halo = NFW(Mvir=self.NU.m_to_Nbody(Mvir), a=self.NU.pos_to_Nbody(rs), c=c)
        a = 3
        b = 0.280
        Mdisc = 6.8
        disc = MyamotoNagai(Mass=self.NU.m_to_Nbody(Mdisc), a=self.NU.pos_to_Nbody(a), b=self.NU.pos_to_Nbody(b))
        ra = 1
        rt = 1.9
        gamma = 1.8
        Mspher = 0.5
        bulge = TruncatedPlaw(Mass=self.NU.m_to_Nbody(Mspher), rt=self.NU.pos_to_Nbody(rt), a=self.NU.pos_to_Nbody(ra), gamma=gamma)
        super().__init__([halo, disc, bulge])

    def get_acceleration(self, particles):
        acc, _, _ = self.acceleration(particles)
        return acc


def generate_stellar_population(N_stars, metallicity, NU, potential):
    Rd = 2.6
    zd = 0.3
    R = np.random.exponential(scale=Rd, size=N_stars)
    phi = np.random.uniform(0, 2 * np.pi, size=N_stars)
    z = np.random.exponential(scale=zd, size=N_stars) * np.random.choice([-1, 1], N_stars)
    positions = np.array([R * np.cos(phi), R * np.sin(phi), z]).T

    velocities = np.zeros_like(positions)
    for i in range(N_stars):
        pos = positions[i]
        radius = np.sqrt(pos[0]**2 + pos[1]**2)
        acc = potential.get_acceleration(Particles(position=np.array([pos]), velocity=np.zeros((1, 3)), mass=np.ones(1)))[0]
        a_R = np.linalg.norm(acc[:2])
        v_circ = np.sqrt(radius * a_R)
        velocities[i, 0] = -v_circ * np.sin(phi[i])
        velocities[i, 1] = v_circ * np.cos(phi[i])
 
    def kroupa_imf(m_min=10, m_max=150, alpha1=1.3, alpha2=2.3):
        m_break = 0.5
        def imf(m):
            return np.where(m < m_break, (m / m_break)**-alpha1, (m / m_break)**-alpha2)
        masses = np.linspace(m_min, m_max, 1000)
        cdf = np.cumsum(imf(masses))
        cdf /= cdf[-1]
        return np.interp(np.random.rand(N_stars), cdf, masses)

    masses = kroupa_imf()
    stars = []
    for i in range(N_stars):
        star = Star(Mzams=masses[i], Z=metallicity)
        star.pos = np.array(positions[i])
        star.vel = np.array(velocities[i])
        stars.append(star)
    particles = Particles(position=positions, velocity=velocities, mass=np.array(masses))
    return stars, particles


def integrator_symplectic_leapfrog(particles, acceleration_estimator, dt):
    acc = acceleration_estimator(particles)
    particles.vel += 0.5 * dt * acc
    particles.pos += dt * particles.vel
    acc = acceleration_estimator(particles)
    particles.vel += 0.5 * dt * acc
    return particles


def integrate_orbits(stars, potential, dt_nbody, T_nbody):
    t = 0
    all_positions = []
    while t < T_nbody:
        positions = np.array([star.pos for star in stars])
        velocities = np.array([star.vel for star in stars])
        masses = np.array([star.getp("Mass", mode="last").iloc[0] for star in stars])
        particles = Particles(position=positions, velocity=velocities, mass=masses)
        particles = integrator_symplectic_leapfrog(particles, potential.get_acceleration, dt_nbody)
        for i, star in enumerate(stars):
            star.pos = particles.pos[i]
            star.vel = particles.vel[i]
        t += dt_nbody
        all_positions.append(particles.pos.copy())
    return stars, np.array(all_positions)


def evolve_stars(stars, t_max, NU, potential):
    pos = []
    tl = []
    mass_list = []
    dt_nb = 0.001
    part = Particles(
        position=np.array([star.pos for star in stars]),
        velocity=np.array([star.vel for star in stars]),
        mass=np.array([NU.m_to_Nbody(star.getp("Mass", mode="last").iloc[0]) for star in stars])
    )
    t = 0
    while t < t_max:
        t += dt_nb
        part = integrator_symplectic_leapfrog(part, lambda p: potential.get_acceleration(p), dt_nb)
        pos.append(part.pos.copy())
        for i, star in enumerate(stars):
            dt_physics = NU.t_to_physical(dt_nb)
            star.evolve_for(dt_physics)
            updated_mass = star.getp("Mass", mode="last").iloc[0]
            part.mass[i] = NU.m_to_Nbody(updated_mass)
            remnant = star.get_remnant().iloc[0]
            remnant_type = int(remnant["RemnantType"])
            if remnant_type in [4, 5, 6]:
                vkick = star.get_SN_kick()["Vkick"]
                if not np.any(np.isnan(vkick)):
                    part.vel[i] += vkick
        tl.append(t)
        mass_list.append(part.mass.copy())
    return np.array(tl), np.array(pos), np.array(mass_list)


def analyze_results(stars, pos):
    neutron_star_positions = []
    black_hole_positions = []
    for star, p in zip(stars, pos[-1]):
        remnant = star.get_remnant().iloc[0]
        remnant_type = int(remnant["RemnantType"])
        if remnant_type == 5:
            neutron_star_positions.append(p)
        elif remnant_type == 6:
            black_hole_positions.append(p)
    neutron_star_positions = np.array(neutron_star_positions)
    black_hole_positions = np.array(black_hole_positions)
    return neutron_star_positions, black_hole_positions


def save_results_to_csv(results, all_positions, initial_positions, filename):
    data = []
    for result in results:
        for position in all_positions:
            config, metallicity, neutron_positions, black_hole_positions = position
            initial_pos = initial_positions[config['sn_kicks']][metallicity]
            if result['config'] == config and result['metallicity'] == metallicity:
                for i, pos in enumerate(neutron_positions):
                    data.append({
                        'Type': 'Neutron Star',
                        'Initial X': initial_pos[i][0],
                        'Initial Y': initial_pos[i][1],
                        'Initial Z': initial_pos[i][2],
                        'Final X': pos[0],
                        'Final Y': pos[1],
                        'Final Z': pos[2],
                        'Metallicity': metallicity,
                        'SN Kicks': config['sn_kicks'],
                        'Kick Std': config['sn_kick_velocity_stdev']
                    })
                for i, pos in enumerate(black_hole_positions):
                    data.append({
                        'Type': 'Black Hole',
                        'Initial X': initial_pos[i][0],
                        'Initial Y': initial_pos[i][1],
                        'Initial Z': initial_pos[i][2],
                        'Final X': pos[0],
                        'Final Y': pos[1],
                        'Final Z': pos[2],
                        'Metallicity': metallicity,
                        'SN Kicks': config['sn_kicks'],
                        'Kick Std': config['sn_kick_velocity_stdev']
                    })
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)


def main():
    configurations = [
        {"sn_kicks": "hobbs", "sn_kick_velocity_stdev": 265.0},
        {"sn_kicks": "hobbs_pure", "sn_kick_velocity_stdev": 265.0},
    ]
    metallicities = [0.02, 0.002]
    results = []
    N_stars = 100
    all_positions = []
    initial_positions = {config['sn_kicks']: {metallicity: [] for metallicity in metallicities} for config in configurations}

    for config in configurations:
        for metallicity in metallicities:
            SEVNmanager.init(config)
            potential = MilkyWayPotential()
            nu_instance = MilkyWayPotential.NU

            stars, particles = generate_stellar_population(N_stars=N_stars, metallicity=metallicity, NU=nu_instance, potential=potential)
            local_initial_positions = [star.pos.copy() for star in stars]
            initial_positions[config['sn_kicks']][metallicity].extend(local_initial_positions)

            t_max = 1.1 * nu_instance.t_to_Nbody(np.min([star.getp("Worldtime", mode="last").iloc[0] for star in stars]))
            tl, pos, masses = evolve_stars(stars, t_max=t_max, NU=nu_instance, potential=potential)
            T_nbody = t_max + 5.0
            stars, final_positions = integrate_orbits(stars, potential, dt_nbody=0.001, T_nbody=T_nbody)
            neutron_positions, black_hole_positions = analyze_results(stars, final_positions)
            
            local_positions = (config, metallicity, neutron_positions, black_hole_positions)
            all_positions.append(local_positions)

            results.append({"config": config, "metallicity": metallicity})

            SEVNmanager.close()

    save_results_to_csv(results, all_positions, initial_positions, 'stellar_simulation_results.csv')


if __name__ == "__main__":
    main()
