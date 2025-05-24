# Run mypy, ruff, and pytest

.PHONY: all mypy ruff test

PACKAGE_NAME = mento

all: mypy ruff test

mypy:
	@echo "Running mypy..."
	mypy . --exclude build

ruff:
	@echo "Running ruff..."
	ruff check . --fix --exclude build

test:
	@echo "Running pytest..."
	pytest tests/
