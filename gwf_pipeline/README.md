# Low frequency mutation in cancer cells

The workflow is created using gwf:
http://gwf.readthedocs.io/en/latest/

gwf can be installed using conda (see link above) or alternatively from source:
https://github.com/gwforg/gwf

Once gwf is installed you should run the following command in this directory to
tell it that jobs should be submitted using slurm:

gwf config set backend slurm

To submit all jobs just type:

gwf run

To check status of jobs type:

gwf status

To check error message of finished job type:

gwf -e logs {job_name}
