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
greenhouse gas concentrations etc...
- Run the model from the command line with `time_loop`, this will automatically print the results to the bash shell

## Example
- The most important insight of this project was to see how different concentrations of CO2 change the temperature of earth's surface
- 400 ppm is our current concentration, if we would double this to 800 ppm the surface temperature would rise around 2 Kelvin
- (If there was no CO2 at all we would live on a snowball)
![CO2](https://github.com/HeckeSiegel/Climate_Model/blob/master/lbl.png)
