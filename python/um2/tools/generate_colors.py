import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

COLOR_HPP = ROOT / "include" / "um2" / "common" / "color.hpp"
C_HEADER = ROOT / "include" / "um2c_colors.h"
C_SOURCE = ROOT / "src" / "um2c_colors.cpp"
PYTHON_COLORS = ROOT / "python" / "um2" / "common_colors.py"

COLOR_PATTERN = re.compile(
    r"inline\s+constexpr\s+Color\s+([A-Za-z_][A-Za-z0-9_]*)\s*"
    r"\{\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\}\s*;"
)


def parse_colors(text: str):
    colors = []
    for match in COLOR_PATTERN.finditer(text):
        name, r, g, b, a = match.groups()
        colors.append(
            {
                "name": name,
                "upper_name": name.upper(),
                "c_name": f"UM2_COLOR_{name.upper()}",
                "r": int(r),
                "g": int(g),
                "b": int(b),
                "a": int(a),
            }
        )
    return colors


def generate_c_header(colors):
    lines = []
    lines.append("/* Generated from include/um2/common/color.hpp. Do not edit by hand. */")
    lines.append("")
    lines.append("#pragma once")
    lines.append("")
    lines.append('#include "um2c.h"')
    lines.append("")
    lines.append("#ifdef __cplusplus")
    lines.append('extern "C" {')
    lines.append("#endif")
    lines.append("")
    for c in colors:
        lines.append(f"extern UM2Color const {c['c_name']};")
    lines.append("")
    lines.append("#ifdef __cplusplus")
    lines.append("}")
    lines.append("#endif")
    lines.append("")
    return "\n".join(lines)


def generate_c_source(colors):
    lines = []
    lines.append("/* Generated from include/um2/common/color.hpp. Do not edit by hand. */")
    lines.append("")
    lines.append('#include "um2c_colors.h"')
    lines.append("#include <um2/common/color.hpp>")
    lines.append("")
    lines.append('extern "C" {')
    lines.append("")
    for c in colors:
        lines.append(f"UM2Color const {c['c_name']} = {{")
        lines.append(
            f"    um2::{c['name']}.r(), um2::{c['name']}.g(), "
            f"um2::{c['name']}.b(), um2::{c['name']}.a()"
        )
        lines.append("};")
        lines.append("")
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def generate_python_file(colors):
    lines = []
    lines.append('"""Generated from include/um2/common/color.hpp. Do not edit by hand."""')
    lines.append("")
    lines.append("from um2 import Color")
    lines.append("")

    for c in colors:
        lines.append(
            f"{c['c_name']} = Color({c['r']}, {c['g']}, {c['b']}, {c['a']})"
        )
    lines.append("")

    lines.append("COMMON_COLORS = {")
    for c in colors:
        lines.append(f'    "{c["name"]}": {c["c_name"]},')
    lines.append("}")
    lines.append("")

    lines.append("__all__ = [")
    for c in colors:
        lines.append(f'    "{c["c_name"]}",')
    lines.append('    "COMMON_COLORS",')
    lines.append("]")
    lines.append("")

    return "\n".join(lines)


def generate_common_colors():
    text = COLOR_HPP.read_text(encoding="utf-8")
    colors = parse_colors(text)

    if not colors:
        raise RuntimeError(f"No colors found in {COLOR_HPP}")

    C_HEADER.write_text(generate_c_header(colors), encoding="utf-8")
    C_SOURCE.write_text(generate_c_source(colors), encoding="utf-8")
    PYTHON_COLORS.write_text(generate_python_file(colors), encoding="utf-8")

    print(f"Generated {len(colors)} colors.")
    print(f"  {C_HEADER}")
    print(f"  {C_SOURCE}")
    print(f"  {PYTHON_COLORS}")


def main():
    generate_common_colors()


if __name__ == "__main__":
    main()