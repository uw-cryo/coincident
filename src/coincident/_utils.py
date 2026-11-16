from __future__ import annotations

import os
from collections.abc import Callable
from functools import wraps
from importlib.util import find_spec

# NOTE: not really sure about type hinting here!
# https://stackoverflow.com/questions/65621789/mypy-untyped-decorator-makes-function-my-method-untyped
from typing import Any, TypeVar, cast

F = TypeVar("F", bound=Callable[..., Any])


def depends_on_optional(module_name: str) -> Callable[[F], F]:
    """borrowed from https://stackoverflow.com/questions/27361427/how-to-properly-deal-with-optional-features-in-python
    Usage: add decorator above function that uses optional dependency
    Example: @depends_on_optional("matplotlib")
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):  # type: ignore[no-untyped-def]
            spec = find_spec(module_name)
            if spec is None:
                message = (
                    f"Optional dependency {module_name} not found ({func.__name__})."
                )
                raise ImportError(message)
            return func(*args, **kwargs)

        return cast(F, wrapper)

    return decorator


# NOTE: Vantor was formerly Maxar, and the API still expects a 'MAXAR-API-KEY' header
VANTOR_API_KEY_HEADER = "MAXAR-API-KEY"


def vantor_api_key() -> str:
    """Vantor API key from VANTOR_API_KEY (or legacy MAXAR_API_KEY) environment variable"""
    return os.environ.get("VANTOR_API_KEY") or os.environ.get("MAXAR_API_KEY") or ""
