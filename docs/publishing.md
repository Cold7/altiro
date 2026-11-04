# Publishing and maintaining orthoresponse

Notes for the maintainer.

## 1. First upload to GitHub

1. Create an **empty** repository named `orthoresponse` on GitHub (no README,
   no license, no .gitignore: they already exist here).
2. Commit and push:

   ```bash
   git remote add origin https://github.com/Cold7/orthoresponse.git
   git push -u origin main
   ```

3. Push the release tag so users can pin it:

   ```bash
   git push origin v0.1.0
   ```

   Then on GitHub: *Releases -> Draft a new release -> choose v0.1.0*.

The workflows in `.github/workflows/` run the Python tests and `R CMD check`
on every push.

## 2. Installing a specific version

```bash
pip install "orthoresponse[all] @ git+https://github.com/Cold7/orthoresponse.git@v0.1.0#subdirectory=python"
```

```r
remotes::install_github("Cold7/orthoresponse@v0.1.0", subdir = "r-package")
```

## 3. Optional: PyPI

```bash
cd python
python -m pip install build twine
python -m build                     # creates dist/*.whl and dist/*.tar.gz
python -m twine upload dist/*       # needs a PyPI account and API token
```

After that, `pip install "orthoresponse[all]"` works without the GitHub URL.
Test first on TestPyPI (`twine upload --repository testpypi dist/*`).

## 4. Optional: R-universe

R-universe builds binary packages from GitHub. Create a repository
`Cold7.r-universe.dev` (named exactly `Cold7.r-universe.dev` on
GitHub) with a `packages.json`:

```json
[{"package": "orthoresponse", "url": "https://github.com/Cold7/orthoresponse", "subdir": "r-package"}]
```

and install the R-universe GitHub app. Users can then run
`install.packages("orthoresponse", repos = c("https://Cold7.r-universe.dev", "https://cloud.r-project.org", BiocManager::repositories()))`.

## 5. Making a new version

1. Update the version in `python/pyproject.toml`,
   `python/src/orthoresponse/__init__.py`, `r-package/DESCRIPTION`,
   `CITATION.cff`, and add an entry to `r-package/NEWS.md`.
2. If the classifier changed, regenerate the shared fixture:
   `python tools/make_fixture.py` (writes to both packages).
3. Run the tests:

   ```bash
   pytest python/tests
   Rscript -e 'roxygen2::roxygenise("r-package")'
   R CMD build r-package && R CMD check --as-cran orthoresponse_*.tar.gz
   ```

4. Commit, tag (`vX.Y.Z`) and push.

## 6. Archiving with a DOI (optional)

Connect the repository to Zenodo (https://zenodo.org, *GitHub* tab), then
create a GitHub release: Zenodo archives it and assigns a DOI you can add to
`CITATION.cff` and to publications.
