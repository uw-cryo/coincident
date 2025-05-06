# Installation

We recommend working with coincident in an isolated environment, managed by
[pixi.sh](https://pixi.sh/latest/):

## Pixi

Alternatively, you can install the latest development version using the
[GitHub CLI](https://cli.github.com) and :

```bash
gh repo clone uw-cryo/coincident
cd coincident
pixi install
pixi shell
```

### Pixi+JupyterLab

## pip

The latest release of `coincident` can be installed from PyPi:

```bash
pip install coincident
```

## conda

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
