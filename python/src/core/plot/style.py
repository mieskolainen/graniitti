# Shared lightweight Matplotlib presentation helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

from __future__ import annotations


# Validate one explicit plot brand
def normalize_plot_brand(value: str) -> str:
    brand = str(value).strip()
    if not brand:
        raise ValueError("Plot brand must not be empty")
    if len(brand) > 64 or any(ord(character) < 32 for character in brand):
        raise ValueError("Plot brand must be a single printable line of at most 64 characters")
    return brand


# Convert one simulator driver identifier into its display brand
def driver_plot_brand(simdriver: str | None) -> str | None:
    if simdriver is None:
        return None
    name = str(simdriver).strip().replace("\\", "/").rsplit("/", 1)[-1]
    if name.lower().endswith(".py"):
        name = name[:-3]
    if "." in name:
        name = name.rsplit(".", 1)[-1]
    lowered = name.lower()
    for suffix in ("_driver", "driver"):
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return normalize_plot_brand(name.upper()) if name else None


# Infer the producer brand from parameter naming when metadata is unavailable
def parameter_plot_brand(param_names) -> str:
    names = [str(name) for name in param_names]
    pandora_prefixes = ("PXML_", "PBOOL_", "PWRAP_")
    if names and all(name.startswith(pandora_prefixes) for name in names):
        return "PANDORA"
    return "GRANIITTI"


# Resolve CLI, saved and driver-derived plot branding in precedence order
def resolve_plot_brand(
    *, override: str | None = None, saved_brand: str | None = None, simdriver: str | None = None, param_names=()
) -> str:
    if override is not None:
        return normalize_plot_brand(override)
    if saved_brand is not None:
        return normalize_plot_brand(saved_brand)
    driver_brand = driver_plot_brand(simdriver)
    if driver_brand is not None:
        return driver_brand
    return parameter_plot_brand(param_names)


# Add collision-free producer and icetune branding to one figure
def add_icetune_branding(
    fig,
    *,
    brand: str = "GRANIITTI",
    x: float = 0.01,
    y: float = 0.995,
    gap_points: float = 5.0,
    brand_fontsize: float = 13.0,
    icetune_fontsize: float = 10.0,
):
    primary = fig.text(x, y, normalize_plot_brand(brand), fontsize=brand_fontsize, fontweight="bold", va="top")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    primary_box = primary.get_window_extent(renderer=renderer)
    primary_right = fig.transFigure.inverted().transform((primary_box.x1, primary_box.y0))[0]
    gap_fraction = float(gap_points) / (72.0 * float(fig.get_figwidth()))
    secondary = fig.text(
        primary_right + gap_fraction, y, "icetune", fontsize=icetune_fontsize, style="italic", va="top"
    )
    return primary, secondary
