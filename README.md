# coincident

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]
[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/uw-cryo/coincident)

<!-- SPHINX-START -->

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/uw-cryo/coincident/workflows/CI/badge.svg
[actions-link]:             https://github.com/uw-cryo/coincident/actions
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/coincident
[conda-link]:               https://github.com/conda-forge/coincident-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/uw-cryo/coincident/discussions
[pypi-link]:                https://pypi.org/project/coincident/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/coincident
[pypi-version]:             https://img.shields.io/pypi/v/coincident
[rtd-badge]:                https://readthedocs.org/projects/coincident/badge/?version=latest
[rtd-link]:                 https://coincident.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->

Identify, curate, process and distribute coincident datasets contained within
select archives across a range of representative terrain and landcover types.
Such datasets are intended to be used by the NASA STV community for
calibration/validation, fusion algorithm development, and discipline-specific
scientific analysis.

**This tool is under active development, there are no stable releases yet!**
https://science.nasa.gov/earth-science/decadal-surveys/decadal-stv/coincident-datasets/

## Development

Use [pixi](https://pixi.sh) for environment management

```bash
git clone https://github.com/uw-cryo/coincident.git
cd coincident
git checkout -b newfeature
pixi shell --environment dev # type `exit` to deactivate
pre-commit install

# Or run pre-configured environments and commands
pixi run test
pixi run lint
pixi run docs
```

## Authentication

Some datasets require authentication to _search_ (Maxar) others only require
authentication to _download_ data (NASA). `coincident` assumes you have the
following Environment Variables defined:

```bash
export EARTHDATA_USERNAME=xxxxx
export EARTHDATA_PASSWORD=yyyyy
export MAXAR_API_KEY=zzzzz
```

### Acknowledgements

- Python packaging template provided by
  https://github.com/scientific-python/cookie
