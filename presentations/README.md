# Counterflow HX research review

Download `Counterflow_HX_Code_and_Results_Review.pptx` for a 30-slide editable
PowerPoint with source notes, charts, code changes and all-case appendix tables.

The deck uses the model-results snapshot `da28165` (October 13, 2026). It
separates historical fits, corrected uncalibrated comparisons, leave-one-case-out
smooth predictions and fitted dimpled transfer. It retains the heat-flux mismatch
and does not claim experimental uncertainty was reduced or physical roughness
validated. All 40 dimpled and 38 smooth cases appear in the appendices.

Rebuild from the repository root with:

```bash
python -m pip install python-pptx pandas numpy matplotlib Pillow
python presentations/build_review.py
```

The slide-generation dependencies are separate from solver dependencies. Charts
are saved as standalone PNG files in `assets/`. Slide notes contain source
references. No solver reruns are required to rebuild this results snapshot.
