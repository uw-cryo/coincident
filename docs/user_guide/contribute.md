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

## Without internet access

In situations where you don't have reliable internet it can be convenient to not
run any checks for package version updates with pixi commands
(`--no-lockfile-update`)

```bash
pixi shell -e dev --no-lockfile-update
```

You can also run only the tests that don't require searching via external APIs:

```bash
pixi run test
```
