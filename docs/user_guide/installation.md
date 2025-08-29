# Installation

## Pixi (recommended)

We recommend working with coincident in an isolated environment, with specific
versions of dependencies managed by [pixi.sh](https://pixi.sh/latest/):

```bash
gh repo clone uw-cryo/coincident
cd coincident
```

```
export PIXI_FROZEN=true
pixi install
pixi shell # type "exit" to deactivate environment
```

### Executing Jupyter Notebooks

You can run and modify coincident example notebooks (in the `docs/examples/`
folder) using your editor of choice. VSCode automatically detects pixi
environments, but JupyterLab requires the
[pixi-kernel extension](https://github.com/renan-r-santos/pixi-kernel).

#### Local workstations

If working on your personal computer it can be convenient to install a global
jupyterlab executable that can be used with any pixi environment:

```
pixi global install --environment jupyterlab jupyterlab jupyterlab-myst pixi-kernel
```

#### Remote workstations

If you're using a remote JupyterHub, like
[CryoCloud JupyterHub](https://hub.cryointhecloud.com/), use the following
command to create an environment for running the `coincident` example notebooks
(you only need to run this once):

```bash
pip install pixi-kernel --user
```

You will then be able to select the `Pixi (Python)` kernel from the JupyterLab
interface. Note that you may need to run `pixi shell` in a terminal before the
environment is available in JupyterLab.

```{note}
Pixi environments are scoped to folders, so you will only be able to use this environment for notebooks under the `coincident/` repository folder.
```

## pip install (not recommended)

The latest release of `coincident` can be installed from PyPi, but this is not
recommended because exact dependency versions are not guaranteed.

```bash
pip install coincident
```

### conda

We do yet have a conda-forge package for coincident. In the meantime, you can
try installing coincident with pip _into_ an existing conda environment, but
dependencies are not guaranteed to work:

```bash
conda activate myenv
pip install coincident
```

Or, if possible, install into a clean conda environment:

```bash
mamba create -n coincident python=3.12 pip
conda activate coincident
pip install coincident
```

## Authentication

Some datasets require authentication to _search_ (Maxar) others only require
authentication to _download_ data (NASA). `coincident` assumes you have the
following Environment Variables defined:

```bash
export EARTHDATA_USERNAME=aaaaa
export EARTHDATA_PASSWORD=bbbbb
export MAXAR_API_KEY=ccccc
```

Sign up for credentials at the following webpages:

- <https://urs.earthdata.nasa.gov>
- <https://developers.maxar.com/docs/authentication/guides/api-key>
- <https://planetarycomputer.developer.azure-api.net>
