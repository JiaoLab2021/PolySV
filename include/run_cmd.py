#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import subprocess
import tempfile
import os


# Run
def run(command, envPath):
    """
    Executes a command in the shell and captures the standard output, standard error, and log.

    Args:
        command (str): The command to be executed.
        envPath (dict): The environment variables to be set for the command.

    Returns:
        tuple: A tuple containing the standard output (str), standard error (str), and log (str).
    """
    # Create output file
    stdout_file = tempfile.NamedTemporaryFile(delete=False)
    stderr_file = tempfile.NamedTemporaryFile(delete=False)
    # Submit task
    proc = subprocess.Popen(command, shell=True, stdout=stdout_file, stderr=stderr_file, env=envPath)
    # Wait for the command to complete
    proc.wait()
    # Read the output file and close the file handle
    with open(stdout_file.name, 'rb') as f:
        stdout_data = f.read()
    os.unlink(stdout_file.name)
    with open(stderr_file.name, 'rb') as f:
        stderr_data = f.read()
    os.unlink(stderr_file.name)

    exit_state = proc.returncode
    stdout = stdout_data.strip().decode('utf-8')
    stderr = stderr_data.strip().decode('utf-8')
    log = ""

    # standard output and error output
    if exit_state != 0 or "FileNotFoundError" in stderr or \
            "command not found" in stderr or \
            "error" in stderr or \
            "Error" in stderr:
        log = f'Error occurred with status code: {exit_state}. Error message: {stderr}.'

    return stdout, stderr, log