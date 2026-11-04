# Running on Perlmutter and Sherlock

`pixi run pipeline` (and the other `pipeline-*` tasks) chooses a profile for the host it runs
on, unless `-p` names one:

| Host | Default profile | Executor | Containers |
| --- | --- | --- | --- |
| Perlmutter (`NERSC_HOST=perlmutter`, or Slurm with `podman-hpc`) | `nersc` | Slurm | podman-hpc |
| Sherlock (Slurm with `apptainer`) | `sherlock` | Slurm | Apptainer |
| Anything else | `local,docker` | local | Docker |

`--queue` chooses where the jobs go, and means what each site calls it:

| Profile | `--queue` is | Default |
| --- | --- | --- |
| `nersc` | a QOS; Perlmutter picks the partition | `shared` |
| `sherlock` | a partition | `normal` |
| `slurm` | a partition | `regular` |

Another Slurm cluster also gets `local,docker`, with a note on stderr; choose its profile with
`-p`. Run the command on a login node from a checkout that compute nodes can read. The
Nextflow manager stays on that login node and submits each task as a Slurm job, so keep it
running for the whole run, for example under `tmux` or `nohup`.

Each run lives in `WORKSPACE/RUN_NAME/`. On both clusters the default workspace is
`$SCRATCH/workspace`, which the site purges: copy the results you keep out of `output/`.

References, and the indexes runs publish beside them, default to project storage, since home
folders are too small: `$CFS/PROJECT/$USER/references` with the `nersc` profile, PROJECT being
your default Slurm account (`sacctmgr show user $USER format=defaultaccount`), and
`$OAK/$USER/genesis/references` with `sherlock`. Elsewhere they default to the checkout's
`references/`. Pass `-R`/`--references` to use another folder; either way, `scripts/run-pipeline.sh`
creates it if missing. Put each sheet's FASTAs and annotations there before a run. The templates
in `conf/sites/` show these settings; copy one outside the checkout and pass it with
`-c /path/site.config`.

## Perlmutter (NERSC)

```bash
cd ~/genesis-preprocessing
pixi run pipeline /path/samples.tsv
```

**Images.** Compute nodes only run images that podman-hpc has pulled and squashed into
`$SCRATCH/storage`; a task that pulls for itself runs out of memory. So, for the `nersc`
profile, `scripts/run-pipeline.sh` first runs `podman-hpc pull` for every image pinned in
`conf/genesis-images.lock`, on the login node. That takes about 3 minutes the first time and
about a minute when the images are already there, and completes any pull that was interrupted
earlier. A launch therefore needs Docker Hub to be reachable. The pins carry no registry and
Perlmutter's Podman only searches `registry.suse.com`, so `podman.registry` is `docker.io`.

**Podman.** Nextflow runs containers with `podman`, but on Perlmutter that must be podman-hpc.
`scripts/run-pipeline.sh` puts `scripts/podman-hpc-bin/` first on `PATH`, whose `podman` runs
`podman-hpc`; Slurm passes `PATH` on to the jobs. Containers run as you
(`--userns=keep-id`), so task files stay readable.

**Account and QOS.** Jobs are charged to your default NERSC project. For another project, or
another QOS for one run, set `SBATCH_ACCOUNT` or `SBATCH_QOS` in your shell; `sbatch` reads them
and they override the job scripts' settings. Avoid setting `process.clusterOptions`: it holds
the profile's QOS and memory settings. Jobs use the `shared` QOS (`--queue`, which is a QOS
here), which is charged only for the cores and memory used and allows at most half a node for 48 hours, so the
profile caps each task at 128 CPUs, 247 GB and 48 hours. Without a QOS, NERSC would use `debug`,
which is limited to 30 minutes.

**Memory.** podman-hpc's image mounts count against a job's memory, so each job asks Slurm for
4 GB more than its task; the container itself stays limited to the task's memory.

**Scheduling.** As [NERSC recommends](https://docs.nersc.gov/jobs/workflow/nextflow/), the
profile checks Slurm every 5 minutes and submits at most 20 jobs a minute, 15 at a time. Each
pipeline stage therefore waits up to 5 minutes, which matters only for short tasks; the `test`
profile (`-p nersc,test`) checks every minute and caps tasks at 1 CPU, 2 GB and 5 minutes. A
login node suits campaigns of a few days; NERSC's `workflow` QOS is meant for longer ones.

**Checking a setup.** A stub run starts every task in its container without running the tools:

```bash
pixi run pipeline -p nersc,test -R /path/references /path/samples.tsv -stub-run
```

## Sherlock (Stanford)

```bash
pixi run pipeline /path/samples.tsv
```

The `sherlock` profile submits to the `normal` partition (`--queue`, a partition here, selects
another) and
runs the images with Apptainer. Nextflow pulls each image once, on the login node, into
`$SCRATCH/workspace/apptainer_cache`, where compute nodes read it. Set the Slurm account, if you
need one, with `SBATCH_ACCOUNT`. References default to `$OAK/$USER/genesis/references`, which
needs `$OAK` set (your group has Oak storage); otherwise pass `-R`/`--references`, e.g. a folder on
`$GROUP_HOME`. Keep results on Oak or `$GROUP_HOME`: `$HOME` is small and `$SCRATCH` is purged.
