# KVP: A multiscale kurtosis approach for seismic phase picking

(in preparation)

To install:

The C sources must be compiled first.

Linux:

```bash

make
pip install .

```

MacOS:

```bash

make mac
pip install .

```

Windows: (requires gcc installed)

```bash

make
pip install .

```

Platform wheels will be provided soon on PyPI.

The examples require both KVP and ObsPy installed, note that the later is not a required dependency of the former.

GitHub pages documentation will be provided soon.

# Repository items

```bash
.
├── kvp                     # Package directory
├── make.bat                # Windows "makefile"
├── Makefile                # Normal makefile
├── pyproject.toml          # pyproject file
├── setup.py                # setup file (enforces platform when building the package)
├── LICENSE
└── README.md
```
