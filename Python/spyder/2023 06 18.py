import meep as mp

sim = mp.Simulation()

eps = 13    # dielectric constant of waveguide
w = 1.2     # width of waveguide
r = 0.36    # radius of holes

# The cell dimensions
sy = 12     # size of cell in y direction (perpendicular to wvg.)
dpml = 1    # PML thickness (y direction only!)
cell = mp.Vector3(1, sy)

b = mp.Block(size=mp.Vector3(1e20, w, 1e20), material=mp.Medium(epsilon=eps))
c = mp.Cylinder(radius=r)

resolution=20

pml_layers = mp.PML(dpml, direction=mp.Y)

fcen = 0.25  # pulse center frequency
df = 1.5     # pulse freq. width: large df = short impulse

s = mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Hz,
              center=mp.Vector3(0.1234,0))

sym = mp.Mirror(direction=mp.Y, phase=-1)

kx = 0.4
sim.k_point = mp.Vector3(kx)

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(mp.Harminv(mp.Hz, mp.Vector3(0.1234), fcen, df)),
        until_after_sources=300)

sim.run(mp.at_every(1/fcen/20, mp.output_hfield_z), until=1/fcen)

k_interp = 19

sim.run_k_points(300, mp.interpolate(k_interp, [mp.Vector3(0), mp.Vector3(0.5)]))

