# Tests for standalone iceweb gallery generation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


# Load the extensionless iceweb entry point for direct helper testing
def load_iceweb_module():
    from importlib import import_module

    return import_module('core.iceweb')


# Verify galleries use the shared template and safely encode dynamic content
def test_create_html_template_encodes_content(tmp_path):
    iceweb = load_iceweb_module()
    output = tmp_path / "html" / "main.html"
    image = output.parent / "plots" / "mass spectrum.png"
    output.parent.mkdir()

    iceweb.create_html(
        {"analysis & data/sub<set": [("M<&>", str(image), "2026-01-02 03:04:05")]},
        str(output),
        image_size=320,
        base_folder_name="unused",
    )

    document = output.read_text(encoding="utf-8")
    assert "__ICEWEB_" not in document
    assert "minmax(320px, 1fr)" in document
    assert "analysis &amp; data" in document
    assert "sub&lt;set" in document
    assert "M&lt;&amp;&gt;" in document
    assert "plots/mass%20spectrum.png" in document
    assert "2026-01-02 03:04:05" in document
