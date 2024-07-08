import os
import subprocess
import sys

def run_command(command):
    """Run a system command and print its output."""
    print(f"Running command: {command}")
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        sys.exit(result.returncode)

def check_conda():
    """Check if Conda is installed."""
    try:
        subprocess.run(["conda", "--version"], check=True)
    except subprocess.CalledProcessError:
        print("Conda is not installed. Please install Conda first.")
        sys.exit(1)

def create_environment(env_file):
    """Create a Conda environment from an environment.yml file."""
    run_command(f"conda env create -f {env_file}")

if __name__ == "__main__":
    check_conda()
    env_file = "environment.yml"
    if os.path.exists(env_file):
        create_environment(env_file)
    else:
        print(f"{env_file} not found. Please provide the correct path to the environment.yml file.")

