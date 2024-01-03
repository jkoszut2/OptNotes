# OptNotes
Notes on optimization and interior-point methods with example code.

## Example Code
### Overview
Example code is included which solves a conic program using an interior-point method. The example Python code is a simple custom solver based on the paper by Lievenberghe[^1]. The example C code is essentially the ECOS solver[^2] with some simplifications.

### Requirements
The following packages are required to run the example Python code: `numpy` and `scipy`.
All dependencies for the example C code are included.

### Simulation Instructions
To run the example code, use the command line to navigate to the location of the `main.py` or `main.c` file. For Python, run `python main.py`. For C, following the notes in the `notes_on_compiling.txt` file to compile and then run the code.

## References
[^1]: Lieven Vandenberghe. “The CVXOPT linear and quadratic cone program solvers” (2010). 
[^2]: Domahidi, Alexander, Eric King-Wah Chu and Stephen P. Boyd. “ECOS: An SOCP solver for embedded systems.” 2013 European Control Conference (ECC) (2013): 3071-3076.