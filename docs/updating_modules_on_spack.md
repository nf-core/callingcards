# Adding new versions to spack configuration files

Spack is good because it allows users to update software versions directly, as
opposed to having to ask the sys admins to do it for us. It is instructive to note
that most of the `module` scripts available on the old partition of HTCF are
very out of date as a result.

That shouldn't be a problem on the new partition because you, the user, are
empowered to keep your lab software up to date now.

__NOTE__: you should do all of this in an interactive session

If you notice that a version you would like to use is unavailable, do the following:

1. Ensure that you are in an interactive session
2. Open the software spack configuration file with the command:

      `spack edit <software_name>`

   For example:

      `spack edit nextflow`

This will open the python script which spack uses to maintain the given software.
It will look something like this:

```
# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

 class Nextflow(Package):
     """Data-driven computational pipelines."""

     homepage = "https://www.nextflow.io"
     url      = "https://github.com/nextflow-io/nextflow/releases/download/v21.04.3/nextflow"

     maintainers = ['dialvarezs']

     version('21.10.6', sha256='104c0352c592924233ea7897cbfb2ece41795be348f97d6dfbc8d66e6271e4ad', expand=False)
     version('21.04.3', sha256='80c7ecd94b55da8eb0e17040dbd0c43ee80e252cd999374e16c00d54d3d3abf3', expand=False)
     version('20.10.0', sha256='54f76c83cbabe8ec68d6a878dcf921e647284499f4ae917356e594d873cb78dd', expand=False)
     version('20.07.1', sha256='de4db5747a801af645d9b021c7b36f4a25c3ce1a8fda7705a5f37e8f9357443a', expand=False)
     version('20.04.1', sha256='b46833ad75b9b7db72668235b53d5c295a9ab02b50d36506bbbe53f383239bde', expand=False)
     version('20.01.0', sha256='fe1900284fd658c0781e6d8048839541afe5818d0b53f6ee8ae81f59d47ad662', expand=False)
     version('19.10.0', sha256='45497eb4bea62dd5477ebe75a6dabfd6905554c46321ca40aec6edfec61c59f4', expand=False)
     version('19.07.0', sha256='e6e7ba4770cd6230bd5410a6fd8c071d6c6dde7a7765880ecabc820b84d38fe5', expand=False)
     version('19.04.1', sha256='21318d8b64095a548f6baf0ef2811f33452e4f9f8a502a46a0aab7815ee34c69', expand=False)
     version('0.25.6', sha256='9498806596c96ba87396194fa6f1d7d1cdb739990f83e7e89d1d055366c5a943', expand=False, deprecated=True)
     version('0.24.1', sha256='0bfde5335b385e3cff99bf4aab619e583de5dc0849767240f675037a2e7c1d83', expand=False, deprecated=True)
     version('0.23.3', sha256='ffe1c314962ff97ebf47b0567883e152522acfbf6fd5800200b1a7a0ca2896d2', expand=False, deprecated=True)
     version('0.21.0', sha256='076089079479da0d91fe1ad7aad06816164ecbcf17f73c55e795b1db8462b28d', expand=False, deprecated=True)
     version('0.20.1', sha256='02635f3371f76a10e12f7366508c90bacf532ab7c23ae03c895317a150a39bd4', expand=False, deprecated=True)
     version('0.17.3', sha256='05563ee1474fbef22f65fa3080792dcb08d218dd1b1561c517ebff4346559dbe', expand=False, deprecated=True)

     depends_on('java')

     def install(self, spec, prefix):
         mkdirp(prefix.bin)

          install(self.stage.archive_file, join_path(prefix.bin, "nextflow"))
          set_executable(join_path(prefix.bin, "nextflow"))
```

Add a new version by adding a line to the configuration file like so:


```
    ...
    version('22.04.5', url='https://github.com/nextflow-io/nextflow/releases/download/v22.04.5/nextflow')  <- THE NEW LINE!
    version('21.10.6', sha256='104c0352c592924233ea7897cbfb2ece41795be348f97d6dfbc8d66e6271e4ad', expand=False)
    version('21.04.3', sha256='80c7ecd94b55da8eb0e17040dbd0c43ee80e252cd999374e16c00d54d3d3abf3', expand=False)
    version('20.10.0', sha256='54f76c83cbabe8ec68d6a878dcf921e647284499f4ae917356e594d873cb78dd', expand=False)
    version('20.07.1', sha256='de4db5747a801af645d9b021c7b36f4a25c3ce1a8fda7705a5f37e8f9357443a', expand=False)
    version('21.10.6', sha256='104c0352c592924233ea7897cbfb2ece41795be348f97d6dfbc8d66e6271e4ad', expand=False)
    ...

```
3. Save and close the file

4. Then, get the sha256 for this new version:

      `spack checksum <software> <version>`

   For example:

      `spack checksum nextflow 2022.05.4`

The return will look something like this:

```
==> Found 1 version of nextflow:

  22.04.5  https://github.com/nextflow-io/nextflow/releases/download/v22.04.5/nextflow

==> Fetching https://github.com/nextflow-io/nextflow/releases/download/v22.04.5/nextflow

    version('22.04.5', sha256='b9155a27e11eef920739ce10db5e1c624951aa8300e2b75d4e43e8a287d566a6')
```

5. spack will return the sha256, which you then need to add into the configuration file, so once
again, do:  

      `spack edit <software_name>`  

And add the line to the version, which you just added. In the rstudio example, this would look like:

```
    ...
    version('22.04.5',
             sha256='b9155a27e11eef920739ce10db5e1c624951aa8300e2b75d4e43e8a287d566a6',
             url='https://github.com/nextflow-io/nextflow/releases/download/v22.04.5/nextflow')  <- THE NEW LINE!
    version('21.10.6', sha256='104c0352c592924233ea7897cbfb2ece41795be348f97d6dfbc8d66e6271e4ad', expand=False)
    version('21.04.3', sha256='80c7ecd94b55da8eb0e17040dbd0c43ee80e252cd999374e16c00d54d3d3abf3', expand=False)
    version('20.10.0', sha256='54f76c83cbabe8ec68d6a878dcf921e647284499f4ae917356e594d873cb78dd', expand=False)
    version('20.07.1', sha256='de4db5747a801af645d9b021c7b36f4a25c3ce1a8fda7705a5f37e8f9357443a', expand=False)
    version('21.10.6', sha256='104c0352c592924233ea7897cbfb2ece41795be348f97d6dfbc8d66e6271e4ad', expand=False)
    ...
```

5. Then, you can install the new version with:

      `spack install <software_name>`

   or more explicitly:

      `spack install <software_name@version_name>`

    In the rstudio example, it would look like so:

    `spack install nextflow@22.04.5`
