# Humbert surfaces

This repository contains code related to the paper [On the geometry of the Humbert surface of square discriminant](https://arxiv.org/) by Sam Frengley. You are free to use this code in your work, please cite the above paper. Let me know if it is useful!

## Requirements
You should have Python 3.8 or newer installed on your system (try `python3 --version`), certainly 3.6 and below will break because of f-strings (and I have mainly been testing on 3.12).

## Dependencies
The python code here relies on matplotlib, sympy, and numpy. See the `.toml` file for the exact versions. If you wish to output the images of intersections of Hirzebruch--Zagier divisors in a neighbourhood of the cusps, I have configured this to work with pdflatex (tested on version 3.141592653-2.6-1.40.26).

## Directory structure
The directory `humbert` contains the python code related to this project. This includes functionality to evaluate the numerical invariants of the surfaces $W_{N,r}$ from the paper. Most of these results are simple evaluations of combinatorially defined invariants. The non-trivial functionality which we provide is the "drawing" capabilities. In particular, we provide functions which use matplotlib to illustrate intersection behaviours of Hirzebruch--Zagier divisors in a neighbourhood of the cusps on $\widetilde{Z}_{N,r}$.

Less important is the `magma-files` directory. This contains some claims which are made in the paper which are verified with computer algebra (please see the specific documentation therein).

## The python implementation

### Installation
A naive installation of the python package is as follows (better computer scientists than I probably have their own preferences):

1. Clone the repository and navigate into the directory:
   ``` shell
   git clone https://github.com/samfrengley/humbert.git
   cd humbert
   ```

2. Check the directory, `ls` should show the `pyproject.toml` file.

3. [optional/recommended] Create a [virtual environment](https://docs.python.org/3/library/venv.html) with 
   ``` shell
   python3 -m venv .venv
   source .venv/bin/activate
   ```

4. Install the `humbert` package and its dependencies:
    ``` shell
    .venv/bin/python3 -m pip install .
    ```

At this point when you run `.venv/bin/python3` you should be able to `import humbert`. You can test that the installation has proceeded as expected by running
``` shell
.venv/bin/python3 tests/install-test.py
```

### Functionality

#### Numerical invariants
As stated above, the primary role of the code is to compute the numerical invariants of the surface which are stated in the paper. An example with $\widetilde{Z}_{N,r}$:
``` python
from humbert.znr import Ztil                                       # import the znr sub-module
Z = Ztil(16,1)                                                     # the surface ~Z(16,1)
print(Z)                                                           # display what we know
```

And an example with $W^{\textsf{small}}_{N,r}$
``` python
from humbert.wnr import WNr_o
W = WNr_o(19, 1)                                                   # W^small_(19,1)
print(W)                                                           # display what we know
```

More examples of this can be found in the documentation in the relevant source code files in `./humbert/` and in the `./scripts/` directory where the code is used to (1) verify some numerical claims in the paper and (2) generate the tables in Appendix B.

#### Graphical display of intersection behaviour
Something which is not used in the paper (though it was useful for finding fibrations) is the ability to generate figures which depict intersection behaviour of Hirzebruch--Zagier divisors and resolutions of cusp singularities. For example compare the output of the following with the image in Figure 16.
``` python
from humbert.znr import Ztil
Z = Ztil(21, 5)
_ = Z.sketch_cusps(disp=True)                                      # Doesn't make pdf
_ = Z.sketch_cusps(tex=True, compile_tex=False, open_pdf=False)    # Makes tex in ./figs/21-5/ but doesn't compile a pdf
_ = Z.sketch_cusps(tex=True, compile_tex=True, open_pdf=True)      # Compiles pdf in ./figs/21-5/
```

##### Non-features
This above snippet illustrates a few to-dos. In particular, it would be nice if the code:

1. worked for modular curves which 'meet the ends' $\widetilde{C}_{\infty,i}$ (for example in the $\widetilde{Z}_{21,5}$ example above the curve $\widetilde{F}_{41}$ is not displayed) -- this shouldn't be too hard,
2. was able to compute the CM intersection points of the Hirzebruch--Zagier divisors and illustrate them (using Lemma 4.11, or translating this to BQFs, cf. Remark 6.11),
3. included the curves $\widetilde{F}_{m \circ g}^+$,
4. computed/drew an intersection matrix/graph for a bunch of known curves (this would give a "more automated" way of finding fibrations).

If this would be useful to your research please let me know. Feel free to create a fork, feature, or pull request. The most ambitious (and probably the "most honest") thing to do would be to implement the code so that:

5. it worked for any positive discriminant (i.e., $D \equiv 0,1 \mod{4}$) and for any level structure. 

## The magma files
The directory `magma-files` contains some simple [magma](http://magma.maths.usyd.edu.au/magma/) scripts which are designed to verify some tedious claims in Appendix A. Files are named accordingly.

## Known issues
The functionality to output a `.pdf` file of the intersection illustrations is almost certainly highly dependent on your LaTex installation. If anyone comes along and actually knows what they're doing, I'd be very grateful of a pull request! In particular, this is almost certainly broken on Windows machines (though I have not tried it to see -- please log an issue if there is one).
