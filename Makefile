# Builds a cache of binaries which can just be copied for CI
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

PYTHON ?= python3
VENV ?= venv

venv: ${VENV}/bin/activate
IN_VENV=. ./${VENV}/bin/activate

$(VENV)/bin/activate:
	test -d $(VENV) || $(PYTHON) -m venv $(VENV) --prompt "duplex_tools"
	${IN_VENV} && pip install pip --upgrade


.PHONY: install
install: venv
	${IN_VENV} && pip install -r requirements.txt
	${IN_VENV} && python setup.py install


.PHONY: develop
develop: venv
	${IN_VENV} && pip install -r requirements.txt
	${IN_VENV} && python setup.py develop


.PHONY: flake8
flake8: venv
	${IN_VENV} && pip install flake8 flake8-rst-docstrings flake8-docstrings flake8-import-order
	${IN_VENV} && flake8 duplex_tools --statistics --import-order-style google --application-import-names duplex_tools


.PHONY: test
test: develop
	${IN_VENV} && pip install pytest pytest-cov hypothesis
	${IN_VENV} && pytest


.PHONY: clean
clean:
	rm -rf dist build duplex_tools.egg-info pypi_build venv

.PHONY: build_env
build_env: pypi_build/bin/activate
IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || $(PYTHON) -m venv pypi_build --prompt "pypi"
	${IN_BUILD} && pip install pip --upgrade
	${IN_BUILD} && pip install --upgrade pip setuptools twine wheel readme_renderer[md]


.PHONY: sdist
sdist: pypi_build/bin/activate
	${IN_BUILD} && python setup.py sdist
