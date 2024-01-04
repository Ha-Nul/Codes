import fdtd
import numpy as np
import matplotlib.pyplot as plt

fdtd.set_backend("numpy")

grid = fdtd.Grid(

    shape = (25e-6, 15e-6, 1),
    grid_spacing = 1.55e-07,
    courant_number= 0.70

)

grid[7.5e-6 : 8.0e-6 , 11.8e-6:13.0e-6, 0] = fdtd.LineSource(
    period = 1550e-9 / (3e8),
    name = "source"
)

grid.run(
    total_time = 100
)

grid.visualize(z=0,animate=False, show=False)