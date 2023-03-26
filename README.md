# CFHackathon2023

This is the code base for the 2023 FAME CF Hackathon, held at Flinders, Victoria Square from March 20-24th. 

There are two ways you can join this repo. Ask Rob to add you as a collaborator and then you can pull/commit/push as usual, or your can fork this repository to your own account and then use pull-requests to add your code. It's up to you!


# Installation

1. Clone this repository.
2. Create a virtual environment

```bash
virtualenv venv
```

3. Activate the virtual environment

```bash
source venv/bin/activate
```

4. Check that you are using the right virtual environment:

```bash
which pip
```

5. Install the requirements:

```bash
pip install -r requirements.txt
```

6. You might need to add this directory to the `PYTHONPATH` (because Rob is too lazy to add a `setup.py` file)

```
export PYTHONPATH=$PYTHONPATH:$PWD
```

You are up and running!

# Collaboration rules!

- Play nice with everyone's code
- If you want to use a python library, please add it (and if you want, the minimum version) to requirements.txt. That means everyone can easily get it installed in their own directories
- If you want to add `tests/` that would be awesome. Then people should run the tests before pushing their code back

