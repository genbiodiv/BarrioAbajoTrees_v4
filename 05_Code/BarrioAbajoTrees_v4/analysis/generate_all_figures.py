#!/usr/bin/env python3
"""Public entry point for the complete Barrio Abajo analysis.

Run from the repository root with::

    python analysis/generate_all_figures.py

The command reads the curated supplementary database, regenerates all eight
figures and result tables, and stops if any census reconciliation check fails.
"""

from thematic_analysis import main


if __name__ == "__main__":
    main()
