import argparse
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib


def parse_config_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to a TOML config file with script-specific settings.",
    )
    return parser.parse_args()


def get_section(config, section_name):
    section = config
    for part in section_name.split("."):
        section = section.get(part, {})
    return section


def config_key_to_constant(key):
    return key.upper()


def apply_config(section_name, namespace):
    args = parse_config_args()
    if args.config is None:
        return

    with args.config.open("rb") as config_file:
        config = tomllib.load(config_file)

    section = get_section(config, section_name)
    for key, value in section.items():
        constant_name = config_key_to_constant(key)
        if constant_name in namespace:
            namespace[constant_name] = value
