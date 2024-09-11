---
name: New source
about: Adding a new type of source and corresponding functionality
title: 'New source: Gateway cloning'
labels: new-source
assignees: ''

---

You can do some of these tasks in parallel, specially if you need to wait for approval, but it's good to keep in mind the intended order:

* [ ] Make a branch in [ShareYourCloning_LinkML](https://github.com/genestorian/ShareYourCloning_LinkML) where you implement the new source (see extra docs in there) on how to do this.
* [ ] Add your branch as a dependency (via commit id) using poetry in this repository, to update `shareyourcloning-linkml` dependency in `pyproject.toml`:
    ```
    poetry add git+hhttps://github.com/genestorian/ShareYourCloning_LinkML#<commit-id>
    ```
* [ ] Implement the new source in a branch in this repository. You will need to add a new endpoint in `main.py` and a new class in `pydantic_models.py` that will handle the new source. You can use the existing sources as a template. Note that you can add extra validation or methods (see examples as well).
* [ ] Write tests for the new source in `tests/test_endpoints.py`. You can use the existing tests as a template.
* [ ] Once the tests pass, merge the PR of the new source in [ShareYourCloning_LinkML] and make a release of the `shareyourcloning-linkml` package.
* [ ] Update the dependency in the branch of this repository repository to the new version.
* [ ] Implement the frontend functionality in a branch, following [the frontend wiki](https://github.com/manulera/ShareYourCloning_frontend/wiki/Checklist-%E2%80%90--adding-a-source).
* [ ] Make tests for the frontend functionality.
* [ ] Merge the PR in this repository into master.
* [ ]  Update the backend submodule in the frontend repository to the latest master.
* [ ]  Merge the PR in the frontend repository into master.
* [ ] You are done!