# Installation

The latest release of `coincident` can be installed from PyPi:

```bash
pip install coincident
```

Alternatively, you can install the latest development version using the
[GitHub CLI](https://cli.github.com) and [pixi.sh](https://pixi.sh/latest/):

```bash
gh repo clone uw-cryo/coincident
cd coincident
pixi install
pixi shell
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
