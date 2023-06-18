# coreEntropy_neighborIFS
The up-to-date and documented code for my research on the connection between the Mandelbrot set and the Thurston set.

One endgoal of this project is to create a RESTful API from which one can obtain all the necessary information regarding a rational angle in the context of the Mandelbrot and Thurston set.  


[THEORY SLIDES](https://slides.com/silvas/qsconjcoreentropy)

## Folder's description
 
- `data`: raw data for the project. 
- `docs`: documentation. (_not yet implemented; there are comments in the docsting of each function, though._)
- `log` : it contains the logging configuration file and eventually the logs 
- `results`: results, including checkpoints, as well as figures and tables. 
- `scripts`: scripts - Python and bash alike - as well as .ipynb notebooks.
- `src`: reusable Python modules for the project. 
- `tests`: tests for the code.

## Some worked out examples (jupyter notebook)
More to come...
- [Generating the Neighbor Graph](/scripts/Generating%20the%20Neighbor%20Graph.ipynb)
- [Parameter Spaces Using the Green Function](/scripts/Parameter%20Spaces%20Using%20The%20Green%20Function.ipynb)
- [The Core Entropy plot](/scripts/The%20Core%20Entropy%20Plot.ipynb)


## Install requirements
> This package uses python version 3.8 (or above, but I haven't tested them yet)

### with pip
1. Create a virtual environment called `mandel_thurston`  (you can call it whatever you want)
    ```shell
    $ python -m venv mandel_thurston
    ``` 
2. Activate it
    ```shell
    $ source mandel_thurston/Scripts/activate
    ``` 
    if you are using PowerShell, activate it with
    ```powershell
    > mandel_thurston/Scripts/Activate.ps1
    ``` 
3. Install the required packages, make sure you are in the root directory of this repo
    ```shell
    $(mandel_thurston) pip install -r requirements.txt
    $(mandel_thurston) pip install -e . # install this package
    ```

### with conda
1. Create a virtual environment and install dependencies
    ```shell
    $ conda env create -f mandel_thurston_env.yml
    ```
2. Activate it
    ```shell
    $ conda activate mandel_thurston
    ```
3. Install this package, make sure you are in the root directory of this repo
    ```shell
    $(mandel_thurston) pip install -e .
    ```


## Run tests
Install the dev requirements
```shell
$(mandel_thurston) pip install -r requirements_dev.txt
```
and then run all tests
```shell
$(mandel_thurston) pytest -vv --cov=src
``` 
- the flag `-vv` is for verbose
- the flag `--cov` is for the coverage report

or run a single test
```shell
$(mandel_thurston) pytest -vv tests/test_functions.py -k "test_core_entropy_exact"
``` 
- the flag `-k` is used to specify the test inside the file 

## TODOs :
- [x] write python code to calculate the Core Entropy 
- [x] fix eigenvalue problem in the core_entropy
- [x] refactor the KneadToRat function
- [x] create python modules from the code in the jupyter notebooks
- [ ] make some example on how to use this code
- [ ] write down the theory
- [ ] figure out another way to draw graphs other than Gephi
- [x] write a recursive function that erases vertices with out degree equal to 0.
- [ ] import Mathematica notebooks (?)
- [ ] write function that obtains the landing parameter in the Mandelbrot set from angle (check Wolf Jung's work)
