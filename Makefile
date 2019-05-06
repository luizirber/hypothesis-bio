all: test

test:
	tox -e py37

build:
	flit build

format:
	tox -e black
