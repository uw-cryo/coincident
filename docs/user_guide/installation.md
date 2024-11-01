# Installation

If you'd like to install into an existing environment we recommend installing
coincident directly from GitHub:

```bash
pip install git+https://github.com/uw-cryo/coincident.git@main
```

Alternatively, you can install a fresh locked environment using the
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
export PC_SDK_SUBSCRIPTION_KEY=ddddd
```

Sign up for credentials at the following webpages:

- https://urs.earthdata.nasa.gov
- https://developers.maxar.com/docs/authentication/guides/api-key
- https://planetarycomputer.developer.azure-api.net
