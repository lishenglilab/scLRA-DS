"""Shared OpenAI-compatible client configuration.

The API key is read from an environment variable by default so that it is not
embedded in Nextflow parameters, generated command scripts, or run logs.
"""

from __future__ import annotations

import argparse
import os
import warnings
from typing import Any
from urllib.parse import urlparse

from openai import OpenAI


DEFAULT_API_KEY_ENV = "LLM_API_KEY"
DEFAULT_MODEL = "deepseek-reasoner"  # Preserved for backward compatibility.


def add_llm_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_model: str = DEFAULT_MODEL,
) -> None:
    """Add the common endpoint, model, and credential arguments."""
    parser.add_argument(
        "--base_url",
        default=os.getenv("LLM_BASE_URL", ""),
        help="Base URL of an OpenAI-compatible API (env: LLM_BASE_URL)",
    )
    parser.add_argument(
        "--model",
        default=os.getenv("LLM_MODEL", default_model),
        help="Model name exposed by the endpoint (env: LLM_MODEL)",
    )
    parser.add_argument(
        "--api_key_env",
        default=DEFAULT_API_KEY_ENV,
        help=f"Environment variable containing the API key [default: {DEFAULT_API_KEY_ENV}]",
    )
    parser.add_argument(
        "--api_key",
        default=None,
        help=argparse.SUPPRESS,
    )


def resolve_api_key(args: argparse.Namespace) -> str:
    """Resolve a key from the environment, with a deprecated CLI fallback."""
    if getattr(args, "api_key", None):
        warnings.warn(
            "--api_key is deprecated because command-line secrets can leak; "
            f"set {args.api_key_env} instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return args.api_key

    api_key = os.getenv(args.api_key_env, "").strip()
    if not api_key:
        raise ValueError(
            f"API key is missing; set the {args.api_key_env} environment variable"
        )
    return api_key


def client_from_args(args: argparse.Namespace) -> OpenAI:
    """Build an OpenAI SDK client for any compatible endpoint."""
    base_url = str(getattr(args, "base_url", "")).strip().rstrip("/")
    if not base_url:
        raise ValueError("LLM base URL is empty; provide --base_url or LLM_BASE_URL")
    parsed_url = urlparse(base_url)
    if parsed_url.scheme not in {"http", "https"} or not parsed_url.netloc:
        raise ValueError("LLM base URL must be an absolute http(s) URL")

    model = str(getattr(args, "model", "")).strip()
    if not model:
        raise ValueError("LLM model is empty; provide --model or LLM_MODEL")
    args.model = model

    return OpenAI(api_key=resolve_api_key(args), base_url=base_url)


def response_text(response: Any) -> str:
    """Extract a non-empty text response with a provider-neutral error."""
    try:
        content = response.choices[0].message.content
    except (AttributeError, IndexError, TypeError) as exc:
        raise RuntimeError("LLM response did not contain a chat-completion choice") from exc

    if not isinstance(content, str) or not content.strip():
        raise RuntimeError("LLM response did not contain non-empty text")
    return content


def response_total_tokens(response: Any) -> int | str:
    """Return token usage when the compatible provider reports it."""
    usage = getattr(response, "usage", None)
    total_tokens = getattr(usage, "total_tokens", None)
    return total_tokens if total_tokens is not None else "unknown"
