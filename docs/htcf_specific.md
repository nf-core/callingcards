
# Running the pipeline on HTCF

You will need the following system dependencies. These can be installed with
`spack`, and therefore only need to be done one time by one person in the lab.

__NOTE__: You may need to update the available versions on spack, in particular
for nextflow, which requires Nextflow `>=21.10.3`. [See this guide](./updating_modules_on_spack.md)

```

$ interactive

$ spack install git

$ spack install openjdk

# note the version requirement -- you may need to add a newer version manually.
# this is something that may need to be done once per lab (not once per user),
# just like any other spack install
$ spack install nextflow

# it is critical that you set `suid=False` for singularity to work on the
# cluster (with any program)
$ spack install singularityce suid=False

```

Follow the instructions in the README quickstart to create a directory to
store the pipeline output, and git pull the repo.

Next, copy and paste the script below into a file called, for example, `run_nf.sh`

```bash
#!/usr/bin/env bash

#SBATCH --mem-per-cpu=10G
#SBATCH -J cc_nf_test.out
#SBATCH -o cc_nf_test.out

# load system dependencies -- on HTCF, we use spack
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow)

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir singularity
mkdir local_tmp

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

# if the repo is in this directory, then this relative path callingcards-mammals
# will work. Otherwise, replace callingcards-mammals with the path, relative or
# absolute, to the correct main.nf
# NOTE: you can also try the test_human profile
nextflow run callingcards-mammals/main.nf -profile test_human,singularity,htcf -resume
```

Launch this with the command

```bash
$ sbatch run_nf.sh

```
At this point, you can check progress with `squeue -u $USER`. You can also
track progress by looking at the nextflow process log, which in this case
will be called `cc_nf_test.out`:

```bash
$ tail -150 cc_nf_test.out
```
Note that it sometimes takes HTCF a long time to start scheduling processes. If
there is no error message in `cc_nf_test.out`, then just keep waiting.

### Running the pipeline on your laptop

See the [README]("../README.md") Quick Start section for instructions on installing
nextflow and one of docker, singularity or conda (conda should be an absolute last choice).
Then follow the instructions in the HTCF tutorial above to create a directory
in which you will clone the callingcards pipeline repository, and launch the pipeline.
Next, copy and paste the script below into a script named something like `run_nf.sh`

```bash

#!/usr/bin/bash

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)
mkdir singularity

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

nextflow run nf-core-callingcards-mammals -profile test_human,singularity,htcf -resume

```
You'll run this with the command `.run_nf.sh` and the process logger will
print to `stdout`.
