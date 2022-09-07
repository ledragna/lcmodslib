Install
=======

From the root folder run:

.. code-block:: bash

    pip install -e . --user

The library will be installed in local and two entry point will be placed in `$USER/.local/bin`: (check the PATH)
    * lmmkinp -- to generate Gaussian input files
    * lmprcdt -- processes the Gaussian fchk and print the local mode results

In the folder utilities examples of bash script to run the jobs are given, after running the jobs the script formats the chk and trim the resulting fchk in order to reduce the disk usage.

.. note:: 
    At the moment the number of energy and property points is expected to be the same.
