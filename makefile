# Run mypy, ruff, and pytest

.PHONY: all mypy ruff test

PACKAGE_NAME = mento

all: mypy ruff test

mypy:
	@echo "Running mypy..."
	mypy .

ruff:
	@echo "Running ruff..."
	ruff .

test:
	@echo "Running pytest..."
	pytest tests/
