# Contributing Guide

Contributions are welcome!

Please feel free to open [issues](https://github.com/uw-cryo/coincident/issues)
in this repository if you come across bugs.

If you'd like to make modifications to the code, please open an issue first to
explain to motivate the changes.

## Setup development environment

We recommend using [pixi](https://pixi.sh) for environment management.

Also use a git branch workflow where "newfeature" is short name related to your
changes.

```bash
git clone https://github.com/uw-cryo/coincident.git
cd coincident
git checkout -b newfeature
pixi shell --environment dev # type `exit` to deactivate
pre-commit install
```

### Using a locked environment

If you are not modifying dependencies, we recommend using the latest locked
environment with pixi. This can be set with an environment variable:

```
export PIXI_FROZEN=true
```

That setting can also be convenient for developing without internet access, as
pixi likes to always check the environment lock-file is up-to-date with
dependency specifications in `pyproject.toml`. Alternatively, you can pass
`--no-lockfile-update` to pixi commands.

### Running tests and checks

```
pixi run test
pixi run networktest # tests requiring access to external APIs
pixi run precommit # also runs automatically upon commits
pixi run lint
pixi run docs
```

### Updating dependencies

If you're new to using pixi for environment management, first have a look at
[their docs](https://pixi.sh/latest/getting_started/).

For dependencies, we use a `latest-up`
[pinning strategy](https://pixi.sh/latest/reference/pixi_configuration/#pinning-strategy),
which is not the pixi default of using `semver` and placing upper bounds on
package versions. You can set this by adding `pinning-strategy = "latest-up"` in
`./.pixi/config.toml` in the coincident folder.

To update the environment, `unset PIXI_FROZEN` if you have that enabled, and add
or upgrade packages via `pixi add [package]` or `pixi upgrade [package]`, or
edit the dependencies in the pyproject.toml file.

## Open a pull request

After you're satisfied with your changes and you've run tests locally, please
open a pull request!
