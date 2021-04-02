# JAX Synteny Browser
An interactive web-based conserved synteny browser application, The Jackson Laboratory (JAX) Synteny Browser. The browser 
allows researchers to highlight or selectively display genome features in the reference and/or the comparison genomes 
based on the biological attributes of the features. The current implementation for the browser supports the reference 
genomes of the laboratory mouse and human.

There is a live, working version of Synteny Browser available at: [syntenybrowser.jax.org](http://syntenybrowser.jax.org/) 

User documentation can be found [here](http://syntenybrowser.jax.org/static/docs/SB-UserManual_v1.pdf).

# Loading a Synteny Browser Database
### Prerequisites
Before starting the setup process, you'll need:

1. A bash terminal (Mac OS X & Linux will have this included) or a way of running shell scripts
2. A version of Python installed on your machine
3. The `pip` Python library installed on your machine (your version of Python might have Pip included by default, but if
not, you'll have to install it manually)

Assuming you have all of the above items, you'll need a Python library called Virtualenv which will allow you to run
Python scripts in isolated environments that can have their own dependencies, versions, and permissions without messing
with those belonging to your machine. To set up Virtualenv open up a bash and navigate to the root `syntenybrowser/`
directory (all of the following commands are run from this directory unless otherwise noted) and install Virtualenv:

    pip install virtualenv
    
*Note: if you're running Python 3, you may have to run `pip install virtualenv`.*


### Getting a Database
To load your own database, you'll need a virtual environment that runs in Python3.7:

    python -m venv venv-db

Once created, activate the virtual environment:

    . venv-db/bin/activate

Install necessary packages:

    pip install -r requirements.txt

Run the database creation script with the required parameter:

    ./create_database.sh synteny.db

This will take several minutes and when it's finished, it will yield a file named 'synteny.db' in the root `synteny-database` directory. Shut down the `venv-db` virtual environment:

    deactivate

