# SECAR_PID code

### User inputs:
1. The reaction informations: isotopes, beam energy, target configuration.
2. Stopper foil thickness.
3. End-detectors' configuration: MCP, IC, DSSD thicknesses and gas pressure.
4. Choose the stopping power calculation toold: VICAR or SRIM. If the user has selected SRIM, the code will guide user how to generate the required files using SRIM and how to name them. If VICAR (http://peiluan-tai.com/programs/stp_lib.html) is selected, there is a C++ code an executable that will generate the files automatically. 

### The code:
1. Generates isotropic reaction angles in COM and translates them to lab frame.
2. Calculates the reaction kinematics for the above generated reaction angles.
3. Calculates energy losses of recoil an the leaky beam through target and stopper foil.
4. Calculates the energy deposited in different sections of the end detectors and the the E-delE particle ID plots.