# Post-processing with Qorona

[Qorona](https://github.com/RayanDhib/Qorona) renders a COCONUT solution into eclipse-like
synthetic imagery: the line-of-sight-integrated perpendicular squashing factor Q⊥, which
lights up the loops, streamers, and current sheets seen at a total solar eclipse, plus
synthetic white-light / polarized-brightness channels. The images are made for direct
morphological comparison with eclipse and coronagraph observations.

It reads COCONUT output directly, both `.CFmesh` restarts and Tecplot `.plt` files, so the
results of the testcases in this directory (e.g. `Eclipse`, `Map`) can be rendered without
any conversion step.

<p align="center">
  <img src="https://raw.githubusercontent.com/RayanDhib/Qorona/main/docs/assets/eclipse.png" width="440" alt="Synthetic eclipse render of a COCONUT corona">
</p>

```bash
pip install qorona

# inspect the solution (mesh, variables, boundaries)
qorona info corona_restart.CFmesh

# build the viewpoint-independent Q-perp volume once (minutes)
qorona build corona_restart.CFmesh -o corona.qor

# render any number of viewpoints off it (seconds each)
qorona render corona.qor -o eclipse.png --fov 8
```

- Documentation: https://rayandhib.github.io/Qorona/
- Quickstart on an example COCONUT corona: https://rayandhib.github.io/Qorona/getting-started/first-eclipse-image/
- Running on a cluster: https://rayandhib.github.io/Qorona/getting-started/hpc/

If Qorona contributes to a publication, please cite it (see
[Citing](https://github.com/RayanDhib/Qorona#citing)).
