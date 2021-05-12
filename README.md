# Utopia Models: Kronecker Network Generator

This repository contains [Utopia] models by Thomas Gaskin

#### Table of Contents
* [Installation](#installation)
* [Model Documentation](#model-documentation)
* [Quickstart](#quickstart)
* [Information for Developers](#information-for-developers)
* [Dependencies](#dependencies)

---

## Installation
**Note:** If not mentioned explicitly, all instructions and considerations from the [main repository][Utopia] still apply here. ☝️

The following instructions assume that you built Utopia in a development environment, as indicated in the framework repository's [README](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia#utopia).

### Step-by-step Instructions
These instructions are intended for **[homebrew]-based** Ubuntu or macOS setups with Utopia already _built_ by your user.

:warning: Be aware that this repository has some additional [dependencies](#dependencies) which might not be satisfiable easily on Ubuntu _without_ [homebrew].
If you don't use homebrew, you will have to fulfil the required dependency versions in some other way.


#### 0 — Setup Utopia
Follow the setup procedure of the Utopia framework repository as described in its [README](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia#utopia).


#### 1 — Clone this repository
Enter the `Utopia` directory into which you previously cloned the main repository.

With access to the Utopia GitLab group, you can clone the repository to that directory using the following command:

```bash
git clone <CLONE-URL>
```

Inside your top level `Utopia` directory, there should now be (at least) two repositories:
* `utopia`: the framework repository
* `KronGen`: *this* repository with your own models


#### 2 — Install dependencies
Install the third-party dependencies using [homebrew].

If you use some other package manager, make sure to fulfil the [dependencies requirements](#dependencies).
Additionally, consider the notes given in the section on [installing dependencies in the Utopia framework repository][Utopia_deps].

```bash
brew update
brew install <additional-dependencies-here>
```

If you want, install the optional dependencies:

```bash
brew install fftw doxygen
```


#### 3 — Configure and build
Enter the repository and create your desired build directory:

```bash
cd models
mkdir build
```

Now, enter the build directory and invoke CMake:

```bash
cd build
cmake ..
```

The terminal output will show the configuration steps, which includes the installation of further Python dependencies and the creation of a virtual environment.

After this, you can build a specific or all Utopia models using:

```bash
make KronGen
make -j4 all    # builds all models, using 4 CPUs
```


#### 4 — Run a model :tada:
You should now be able to run Utopia models from both the framework and this repository via the Utopia CLI.

Upon configuration, CMake creates symlinks to the Python virtual environment of the Utopia repository.
You can enter it by navigating to either of the build directories and calling

```bash
source ./activate
```

You can now execute the Utopia frontend as you are used to, only now the models of this repository should be available.


### Troubleshooting
* If the `cmake ..` command fails because it cannot locate the Utopia main
    repository, something went wrong with the CMake User Package Registry.
    You can always specify the location of the package _build_ directory via
    the CMake variable `Utopia_DIR` and circumvent this CMake feature:

    ```bash
    cmake -DUtopia_DIR=<path/to/>utopia/build ..
    ```

* If your model requires a certain branch of the Utopia framework, make sure
    to re-configure and re-build the framework repository after switching
    the branch.



## Model Documentation
Unlike the [main Utopia documentation](https://hermes.iup.uni-heidelberg.de/utopia_doc/latest/html/), the models included in this repository come with their own documentation which has to be built locally.
It is *not* available online.

To build these docs locally, navigate to the `build` directory and execute

```bash
make doc
```

The [Sphinx](http://www.sphinx-doc.org/en/master/index.html)-built user documentation will then be located at `build/doc/html/`, and the C++ [doxygen](http://www.stack.nl/~dimitri/doxygen/)-documentation can be found at `build/doc/doxygen/html/`.
Open the respective `index.html` files to browse the documentation.



## Quickstart
### How to run a model?
The Utopia command line interface (CLI) is, by default, only available in a Python virtual environment, in which `utopya` (the Utopia frontend) and its dependencies are installed.
This virtual environment is located in the **main** (`utopia`) repository's build directory.
However, symlinks to the `activate` and `run-in-utopia-env` scripts are provided within the build directory of this project.
You can enter it by specifying the correct path:

```bash
source <path/to/{utopia,models}>/build/activate
```

Now, your shell should be prefixed with `(utopia-env)`.
All the following should take place inside this virtual environment.

As you have already done with the `dummy` model, the basic command to run a model named `KronGen` is:

```bash
utopia run KronGen
```

You can list all models registered in the frontend with

```bash
utopia model ls
```



## Information for Developers
### New to Utopia? How do I implement a model?
Please refer to the [documentation of the framework repository](https://hermes.iup.uni-heidelberg.de/utopia_doc/latest/html/).


### Testing
Not all test groups of the main project are available in this repository.

| Identifier        | Test description   |
| ----------------- | ------------------ |
| `model_<name>`    | The C++ model tests of model with name `<name>` |
| `models`†         | The C++ and Python tests for _all_ models |
| `models_python`†‡ | All python model tests (from `python/model_tests`) |
| `all`             | All of the above. (Go make yourself a hot beverage, when invoking this.) |

_Note:_
* Identifiers marked with `†` require all models to be built (by running `make all`).
* Identifiers marked with `‡` do _not_ have a corresponding `build_tests_*` target.
* Notice that you cannot execute tests for models of the main project in this repository.

#### Evaluating Test Code Coverage
Code coverage is useful information when writing and evaluating tests.
The coverage percentage of the C++ code is reported via the GitLab CI pipeline.
Check the [`README.md` in the main repository](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia#c-code-coverage) for information on how to compile your code with coverage report flags and how to retrieve the coverage information.


## Dependencies
| Software                  | Version  | Comments |
| ------------------------- | -------- | -------- |
| [FFTW] _(optional)_       | >= 3.3.8 | for fast Fourier transformations |
| [doxygen] _(optional)_    |          | for building the documentation |

These dependencies come on top of [Utopia's dependencies][Utopia_deps].
For some dependencies, a more recent version number is required than by the Utopia framework; the numbers given here take precedence and are enforced.
Note that the homebrew-based CI/CD testing image might use more recent versions, so it's actually best to have those installed.

:warning: There seem to be some weird errors with boost v1.74, observed only in the pipeline, thus the testing image uses v1.71.


### Additional Python Dependencies
If your models require additional Python packages (e.g., for plotting and testing), they can be added to the `python/requirements.txt` file.
All packages listed there will be automatically installed into the Python virtual environment as part of the CMake configuration.

Note that if such an installation fails, it will _not_ lead to a failing CMake configuration.


[Utopia]: https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia
[Utopia_deps]: https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia#dependencies

[homebrew]: https://brew.sh
[armadillo]: http://arma.sourceforge.net/
[fmt]: https://fmt.dev/latest/index.html
[spdlog]: https://github.com/gabime/spdlog
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[graphviz]: https://www.graphviz.org/

[FFTW]: http://www.fftw.org/
[range-v3]: https://ericniebler.github.io/range-v3/
[doxygen]: http://www.doxygen.nl
