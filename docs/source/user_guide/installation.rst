============
Installation
============

.. Pre-requisites should be automagically handled so I'm not sure they need to be listed here?
.. 
.. Pre-requisites
.. --------------

.. Espresso requires Python 3.6+, and the following dependencies:

.. - numpy>=1.18
.. - scipy>=1.0.0
.. - matplotlib>=3.1
.. - tqdm>=4.0

**Step 1**: (*Optional*) Set up a virtual environment.

We strongly recommend installing Espresso within a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`_. This ensures that Espresso can install the various modules that it needs without the risk of breaking anything else on your system. There are a number of tools that can facilitate this, including `venv`, `virtualenv`, `conda` and `mamba`.

.. tabbed:: venv

  Ensure you have `python>=3.6`. Then, you can create a new virtual environment by running the command

  .. code-block:: console

    $ python -m venv <path-to-new-env>/<env-name>

  where :code:`<path-to-new-env>` is your prefered location for storing information about this environment, and :code:`<env-name>` is your preferred name for the virtual environmment. For example,

  .. code-block:: console

    $ python -m venv ~/my_envs/espresso 

  will create a virtual environment named :code:`espresso` and store everything within a sub-directory of your home-space named :code:`my_envs`.

  To 'activate' or 'switch on' the virtual environment, run the command
  
  .. code-block:: console

    $ source <path-to-new-env>/<env-name>/bin/activate

  At this point you effectively have a 'clean' Python installation. You can now install and use Espresso, following the instructions at step 2. When you are finished, you can run the command
  
  .. code-block:: console

    $ deactivate

  and your system will return to its default state. If you want to use Espresso again, simply re-run the 'activate' step above; you do not need to repeat the installation process. Alternatively, you can remove Espresso and the virtual environment from your system by running

  .. code-block:: console

    $ rm -rf <path-to-new-env>/<env-name>

.. tabbed:: virtualenv

  You can create a new virtual environment (using Python version 3.10) by running the command

  .. code-block:: console

    $ virtualenv <path-to-new-env>/<env-name> -p=3.10
  
  where :code:`<path-to-new-env>` is your prefered location for storing information about this environment, and :code:`<env-name>` is your preferred name for the virtual environmment. For example,

  .. code-block:: console

    $ virtualenv ~/my_envs/espresso -p=3.10

  will create a virtual environment named :code:`espresso` and store everything within a sub-directory of your home-space named :code:`my_envs`.

  To 'activate' or 'switch on' the virtual environment, run the command

  .. code-block:: console

    $ source <path-to-new-env>/<env-name>/bin/activate

  At this point you effectively have a 'clean' Python installation. You can now install and use Espresso, following the instructions at step 2. When you are finished, you can run the command

  .. code-block:: console

    $ deactivate

  and your system will return to its default state. If you want to use Espresso again, simply re-run the 'activate' step above; you do not need to repeat the installation process. Alternatively, you can remove Espresso and the virtual environment from your system by running

  .. code-block:: console

    $ rm -rf <path-to-new-env>/<env-name>

.. tabbed:: conda / mamba

  You can create a new virtual environment (using Python version 3.10) by running the command

  .. code-block:: console

    $ conda create -n <env-name> python=3.10

  where :code:`<env-name>` is your preferred name for the virtual environmment. For example,

  .. code-block:: console

    $ conda create -n espresso python=3.10

  will create a virtual environment named :code:`espresso`.
  
  To 'activate' or 'switch on' the virtual environment, run the command

  .. code-block:: console

    $ conda activate <env-name>

  At this point you effectively have a 'clean' Python installation. You can now install and use Espresso, following the instructions at step 2. When you are finished, you can run the command
  
  .. code-block:: console

    $ conda deactivate

  and your system will return to its default state. If you want to use Espresso again, simply re-run the 'activate' step above; you do not need to repeat the installation process. Alternatively, you can remove Espresso and the virtual environment from your system by running
  
  .. code-block:: console

    $ conda env remove -n <env-name>


**Step 2**: Install Espresso

.. tabbed:: pip

  Espresso is available on `PyPI <https://pypi.org/project/geo-espresso/>`_, so for most users installation is as simple as:

  .. code-block:: console

    $ pip install geo-espresso

.. tabbed:: From source

  You can build Espresso from source. You are most likely to want to do this if you want to work in 'developer mode', and make changes to Espresso's source code.

  .. code-block:: console

    $ git clone https://github.com/inlab-geo/espresso.git
    $ cd espresso
    $ pip install -e .

  The :code:`-e` flag ensures that the module is installed in editable mode; you can omit this if you do not intend to make any changes.

If all has gone well, you should now be able to successfully :code:`import espresso` within your Python interpreter or script.