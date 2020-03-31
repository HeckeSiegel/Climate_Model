# Climate Model
> This model grew organically over the course of 4 months as part of an Atmospheric Physics class.
It is a 1D model of the earth's atmosphere which calculates temperatures of the different atmospheric layers from the surface to the top of 
the atmosphere averaged over the whole globe.

## Physics
- Start with random temperature distribution over all layers
- Calculate change of temperature due to radiative transfer based on [Schwarzschild's equation](https://en.wikipedia.org/wiki/Schwarzschild%27s_equation_for_radiative_transfer) for each time step and each layer
- Update temperatures until equilibrium is reached

## Installation
- This is easiest to use in a Linux environment
- `git clone https://github.com/HeckeSiegel/Climate_Model.git`
- In the time_loop.c file you can decide which kind of atmosphere you want to simulate and change different parameters like number of layers,
greenhouse gas concentrations...
- Run the model from the command line with `time_loop`, this will automatically print the results to the bash shell
