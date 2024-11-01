# Contributing Guide

Contributions are welcome!

We recommend use [pixi](https://pixi.sh) for environment management

```bash
git clone https://github.com/uw-cryo/coincident.git
cd coincident
git checkout -b newfeature
pixi shell --environment dev # type `exit` to deactivate
pre-commit install

# Or run pre-configured commands:
pixi run networktest # or 'test'
pixi run precommit # also runs automatically upon commits
pixi run lint
pixi run docs
```
