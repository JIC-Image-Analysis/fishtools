fishtools
=========

A Python package for analysing microscopy data for FISH.

Installation
------------

Instructions for installation in Linux or Linux environment in Windows (WSL)

step 1 - Clone repository
----------------

Clone or download the code repository.

step 2 - Install bftools
----------------

Next install bioformats tools:

.. code-block:: bash

    bash install-bftools.sh
    export PATH=$PATH:~/packages/bftools
    bfconvert -version

The last command checks that the install has worked by printing the bioformats
version.

If there is a problem running this, you may need additional software. Try installing the additional needed software if missing:

.. code-block:: bash

    sudo apt update
    sudo apt upgrade
    sudo apt install unzip
    sudo apt install openjdk-11-jdk

(sudo means you have to give the password)
Now try the commands to install bioformats tools again.

step 3 - Create Python virtualenv
------------------------


Test you have working python

.. code-block:: bash

    python --version

You might need to write "python3" instead of "python".
This should print out the python version.
Make sure the python version was 3.

Create the Python virtual environment:

.. code-block:: bash

    python -m venv venv
    source venv/bin/activate

If this does not work, you may have to install some software:

.. code-block:: bash

    sudo apt install python-venv

Install requirements:

.. code-block:: bash

    pip install -r requirements.txt

This might take a while.

step 4 - Test image conversion
---------------------

To test the image conversion, first create a folder called ``local-data``. Then
copy some test .czi files into that folder. Then, to test the conversion:

.. code-block:: bash

    mkdir scratch
    python scripts/quick_convert_to_ids.py local-data scratch myids

step  - OPTIONAL
---------------------

For Windows WSL
In case you need to use a folder not in your home directory you need to mount it so that your computer recognises it.
The first time you need to do this (e.g. for a hypothetical drive Z)

.. code-block:: bash

    mkdir /mnt/z

Each time you restart WSL you may have to do this:

.. code-block:: bash

    sudo mount -t drvfs Z: /mnt/z

To go into this folder:

.. code-block:: bash

    cd /mnt/z/
