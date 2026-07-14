from __future__ import annotations

import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import anchorscope


class ReleaseMetadataTests(unittest.TestCase):
    def test_version_is_consistent(self):
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        project_section = re.search(
            r"(?ms)^\[project\]\s*(.*?)(?=^\[|\Z)", pyproject
        )
        self.assertIsNotNone(project_section, "missing [project] section")
        version_match = re.search(
            r'(?m)^version\s*=\s*"([^"]+)"\s*$', project_section.group(1)
        )
        self.assertIsNotNone(version_match, "missing project version")
        version = version_match.group(1)
        self.assertEqual(version, anchorscope.__version__)
        self.assertIn(f"version: {version}", (ROOT / "CITATION.cff").read_text(encoding="utf-8"))

    def test_json_schema_and_v2_example_are_valid_json(self):
        json.loads((ROOT / "src" / "anchorscope" / "config.schema.json").read_text(encoding="utf-8"))
        json.loads((ROOT / "examples2" / "config_v2.json").read_text(encoding="utf-8"))

    def test_publication_documents_exist(self):
        for relative in (
            "LICENSE",
            "CITATION.cff",
            "CHANGELOG.md",
            "docs/METHODS.md",
            "docs/PUBLICATION_CHECKLIST.md",
        ):
            self.assertTrue((ROOT / relative).is_file(), relative)


if __name__ == "__main__":
    unittest.main()
