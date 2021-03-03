# NepoSolver
## Intro
I've made this in order to make my life easier studying Electrical Power Systems.
The name is based on my professor, Leonardo Nepomuceno, and I've used this set of tools to study and
solve his exams.

## nepo_tools.py
This is NepoSolver's main library. All of the tools used are here. You can take a look inside to
see what each function does. Hopefully I have documented them all.
To import the file from the working directory, usually it is common to do:

    import lib.nepo_tools as nepo

    np.set_printoptions(precision=4)
    np.set_printoptions(suppress=True)

The printout settings round the results to 4 decimal places and suppress scientific notation.
## NepoSolver_Octave_MatLab.m
This is what I originally created before I decided to use Python for good. I do not plan to update this file.
## Requirements: 
* Python 3.5+
* Numpy
## To be done:
* easier way to add arrays (through string perhaps, octave like etc.)
* Add language change (between Portuguese and English)
* ~~non symmetrical faults~~ 
* ~~Conversion of matrix to polar form~~
* Add how to use to readme
