# ARCHIVED: SEE https://github.com/Lab41/hypothesis-bio FOR UPDATED REPO

[![Build Status][status]][status_link]
[![codecov][codecov]][codecov_link]

[status]: https://travis-ci.com/luizirber/hypothesis-bio.svg?branch=master
[status_link]: https://travis-ci.com/luizirber/hypothesis-bio
[codecov]: https://codecov.io/gh/luizirber/hypothesis-bio/branch/master/graph/badge.svg
[codecov_link]: https://codecov.io/gh/luizirber/hypothesis-bio

# Hypothesis-bio

This module provides a [Hypothesis][Hypothesis] strategy for generating biological data formats.
This can be used to efficiently and thoroughly test your code.

## Installation

This module can be installed via `pip`:
```
pip install hypothesis-bio
```

## Development

Tests can be executed using tox:
`tox -e py37` runs all tests in a Python 3.7 environment.

There are also [pre-commit][pre-commit] files available, so you can run
`pre-commit install` to add them to your repo.
(These are enforced by CI)

[pre-commit]: https://pre-commit.com/

## See also

This module is based on [Hypothesis-networkx][nx], [hypothesis-csv][csv] and
[hypothesis-jsonschema][jsonschema].

[nx]: https://github.com/pckroon/hypothesis-networkx/
[csv]: https://github.com/chobeat/hypothesis-csv
[jsonschema]: https://github.com/Zac-HD/hypothesis-jsonschema
[Hypothesis]: https://hypothesis.readthedocs.io/en/latest/index.html
