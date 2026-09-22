"""IPython's ``display`` and ``Markdown``, imported when something is shown.

The views render Markdown in a notebook, and every element imports the views;
importing IPython at the top of them made it part of ``import mento``. It is a
large import for a calculation that never shows anything, and an environment
that runs mento headless -- a script, a server, Pyodide -- need not have it at
all. So the two names live here as thin stand-ins: they reach for IPython when
called, and fall back to plain text on ``print`` when it is not installed.
"""

from __future__ import annotations

from typing import Any


class _PlainMarkdown:
    """What ``Markdown`` returns without IPython: the text, and nothing else.

    It keeps the ``data`` attribute of IPython's class, so code that reads the
    source back works against either.
    """

    def __init__(self, data: str) -> None:
        self.data = data

    def __str__(self) -> str:
        return self.data


def Markdown(data: str) -> Any:
    """An ``IPython.display.Markdown`` of ``data``, or its plain-text stand-in."""
    try:
        from IPython.display import Markdown as _Markdown
    except ImportError:
        return _PlainMarkdown(data)
    return _Markdown(data)  # type: ignore[no-untyped-call]


def display(obj: Any) -> None:
    """``IPython.display.display``, or ``print`` when IPython is not installed."""
    try:
        from IPython.display import display as _display
    except ImportError:
        print(obj)
        return
    _display(obj)  # type: ignore[no-untyped-call]
